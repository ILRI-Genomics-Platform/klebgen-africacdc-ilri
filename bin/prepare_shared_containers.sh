#!/usr/bin/env bash
#
# prepare_shared_containers.sh - fill and lock the shared, read-only Singularity/Apptainer
# image folder used with `-profile singularity -c conf/slurm.config` (or conf/local.config).
# Run once by a trainer (not participants).
#
#   bin/prepare_shared_containers.sh [--dir DIR] [--group GROUP] [--check] [--no-lock]
#
#   --dir DIR      image folder (default: /var/scratch/global/jjuma/neonatal-amr/shared/containers). Use the same value as
#                  shared_container_dir in your site config (conf/slurm.config or conf/local.config).
#   --group GROUP  make images readable only by this group (e.g. the course's Unix group):
#                  folders 0550, files 0440. Without it, everyone can read: folders 0555, files 0444.
#                  Nobody, including the owner, keeps write permission.
#   --check        only report which images are present or missing (exit 1 if any are missing);
#                  changes nothing.
#   --no-lock      leave permissions unchanged (for repeated test runs).
#
# What it does
#   1. Asks Nextflow which container every process uses under -profile singularity with
#      conf/singularity_images.config (nextflow inspect), so the list always matches the
#      pipeline code and configuration. Optional steps (Prokka, Bakta, --phylo_method ska2) are included.
#   2. nf-core module images (https://depot.galaxyproject.org/..., Seqera community images):
#      downloads each into DIR/nf-core/ under the file name Nextflow looks for.
#   3. Local module images (DIR/local/*.sif): builds each from its module's environment.yml,
#      and the course-only tool images (DIR/course/<tool>.sif: seqkit, FastQC, CheckM2, Bandage, Bakta)
#      from course/envs/<tool>/environment.yml, for the hands-on sessions (course/tools.sh).
#      with `singularity build`. Building needs root or --fakeroot; override the command with
#      BUILD_CMD, e.g.  BUILD_CMD="sudo singularity build"  or  BUILD_CMD="apptainer build --fakeroot".
#   4. Writes DIR/MANIFEST.tsv (image, source, size, sha256, processes) as a record for the course.
#   5. Removes write permission everywhere under DIR.
#
# Needs: nextflow (>= 25.04), python3, curl, and singularity or apptainer (for building and for
# any docker:// images). Network access is needed only here, never during the course.
#
set -euo pipefail

DIR=/var/scratch/global/jjuma/neonatal-amr/shared/containers
GROUP=""
CHECK_ONLY=0
LOCK=1
BUILD_CMD=${BUILD_CMD:-"singularity build --fakeroot"}
MICROMAMBA_IMAGE=${MICROMAMBA_IMAGE:-"docker://mambaorg/micromamba:1.5.10"}

while [ $# -gt 0 ]; do
    case "$1" in
        --dir)     DIR="$2"; shift 2 ;;
        --group)   GROUP="$2"; shift 2 ;;
        --check)   CHECK_ONLY=1; shift ;;
        --no-lock) LOCK=0; shift ;;
        -h|--help) awk 'NR == 1 {next} /^#/ {sub(/^# ?/, ""); print; next} {exit}' "$0"; exit 0 ;;
        *) echo "Unknown option: $1 (see --help)" >&2; exit 2 ;;
    esac
done

PIPELINE_DIR=$(cd "$(dirname "$0")/.." && pwd)
WORK=$(mktemp -d "${TMPDIR:-/tmp}/prepare_containers.XXXXXX")
trap 'rm -rf "$WORK"' EXIT

for tool in nextflow python3; do
    command -v "$tool" >/dev/null || { echo "ERROR: '$tool' not found on PATH" >&2; exit 1; }
done

DIR=$(python3 -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "$DIR")
echo "Image folder: $DIR"

# ---------------------------------------------------------------------------------------------
# 1. Which image does each process use?
# ---------------------------------------------------------------------------------------------
echo "Resolving containers with nextflow inspect (profile singularity + conf/singularity_images.config)..."
# Twice: once with the default phylogeny method (Snippy) and once with --phylo_method ska2, with the optional
# annotation steps switched on, so that every image any run can need is prepared.
for method in snippy ska2; do
    ( cd "$WORK" && nextflow -q inspect "$PIPELINE_DIR/main.nf" \
          -profile test,singularity -c "$PIPELINE_DIR/conf/singularity_images.config" \
          --shared_container_dir "$DIR" --run_prokka true --run_bakta true --bakta_db "$WORK" \
          --phylo_method "$method" -format json ) > "$WORK/inspect_$method.json"
done
python3 - "$WORK/inspect_snippy.json" "$WORK/inspect_ska2.json" > "$WORK/inspect.json" <<'MERGE'
import json, sys
procs = {}
for path in sys.argv[1:]:
    for p in json.load(open(path))["processes"]:
        procs.setdefault(p["name"], p)
print(json.dumps({"processes": list(procs.values())}))
MERGE

# Plan: one line per image -> kind <tab> source <tab> target file <tab> processes
python3 - "$WORK/inspect.json" "$DIR" "$PIPELINE_DIR" > "$WORK/plan.tsv" <<'PY'
import json, os, sys

inspect, dest, pipeline = sys.argv[1], sys.argv[2], sys.argv[3]
procs = json.load(open(inspect))["processes"]

def simple_name(url):
    """File name Nextflow uses for an image in singularity.libraryDir / cacheDir."""
    p = url.find("://")
    name = url[p + 3:] if p != -1 else url
    ext = ".img"
    if ".sif:" in name:
        ext, name = ".sif", name.replace(".sif:", "-")
    elif name.endswith(".sif"):
        ext, name = ".sif", name[:-4]
    return name.replace(":", "-").replace("/", "-") + ext

images = {}
for p in procs:
    c = (p.get("container") or "").strip()
    if not c:
        sys.exit(f"ERROR: process {p['name']} has no container (check conf/singularity_images.config)")
    if c.startswith(dest + "/local/"):
        module = os.path.join(pipeline, "modules", "local", *p["name"].lower().split("_"))
        env = os.path.join(module, "environment.yml")
        if not os.path.exists(env):
            sys.exit(f"ERROR: no environment.yml for {p['name']} (looked in {module})")
        key = ("build", env, c)
    elif c.startswith(("http://", "https://")):
        key = ("download", c, os.path.join(dest, "nf-core", simple_name(c)))
    else:
        src = c if "://" in c else "docker://quay.io/" + c   # nextflow.config: singularity.registry = 'quay.io'
        key = ("pull", src, os.path.join(dest, "nf-core", simple_name(src)))
    images.setdefault(key, []).append(p["name"])

# Course-only tools (hands-on sessions, not the pipeline): course/envs/<tool>/environment.yml
# -> DIR/course/<tool>.sif, listed in MANIFEST.tsv as process COURSE_<TOOL> (course/tools.sh looks them up)
course_envs = os.path.join(pipeline, "course", "envs")
if os.path.isdir(course_envs):
    for name in sorted(os.listdir(course_envs)):
        env = os.path.join(course_envs, name, "environment.yml")
        if os.path.exists(env):
            images.setdefault(("build", env, os.path.join(dest, "course", name + ".sif")), []).append("COURSE_" + name.upper())

# Several local modules may share one image (e.g. KAPTIVE uses KLEBORATE's): build it once,
# after checking that their environment.yml files ask for the same packages
merged = {}
for (kind, src, target), names in images.items():
    if target in merged:
        prev_kind, prev_src, prev_names = merged[target]
        if kind == "build" and open(prev_src).read().split("dependencies:")[-1] != open(src).read().split("dependencies:")[-1]:
            sys.exit(f"ERROR: {target} is used by {prev_names + names} but their environment.yml files differ")
        merged[target] = (prev_kind, prev_src, prev_names + names)
    else:
        merged[target] = (kind, src, names)

for target, (kind, src, names) in sorted(merged.items()):
    print(kind, src, target, ",".join(sorted(names)), sep="\t")
PY

n_total=$(wc -l < "$WORK/plan.tsv")
echo "$n_total images needed for $(python3 -c 'import json,sys;print(len(json.load(open(sys.argv[1]))["processes"]))' "$WORK/inspect.json") processes."

# ---------------------------------------------------------------------------------------------
# --check: report and stop
# ---------------------------------------------------------------------------------------------
report() {
    local missing=0
    printf "%-9s %-8s %s\n" "STATUS" "KIND" "IMAGE"
    while IFS=$'\t' read -r kind src target procs; do
        if [ -s "$target" ]; then status=present; else status=MISSING; missing=$((missing + 1)); fi
        printf "%-9s %-8s %s\n" "$status" "$kind" "${target#$DIR/}"
    done < "$WORK/plan.tsv"
    echo "$((n_total - missing)) of $n_total images present in $DIR"
    return $missing
}

if [ "$CHECK_ONLY" -eq 1 ]; then
    report && exit 0 || exit 1
fi

# ---------------------------------------------------------------------------------------------
# 2-3. Download and build what is missing
# ---------------------------------------------------------------------------------------------
owner_of() { stat -c %U "$1" 2>/dev/null || stat -f %Su "$1" 2>/dev/null || echo unknown; }

# Nearest existing parent of a path
existing_parent() {
    local p="$1"
    while [ ! -e "$p" ]; do p=$(dirname "$p"); done
    echo "$p"
}

if [ ! -e "$DIR" ]; then
    if ! mkdir -p "$DIR" 2>/dev/null; then
        parent=$(existing_parent "$DIR")
        echo "ERROR: the image folder $DIR does not exist, and you cannot create it:" >&2
        echo "       you have no write permission in $parent (owner: $(owner_of "$parent"))." >&2
        echo "       Either ask the owner of $parent to create $DIR for you, or choose a folder you can" >&2
        echo "       create, e.g.  --dir /var/scratch/global/$USER/neonatal-amr/shared/containers" >&2
        echo "       and set the same path as shared_container_dir in your site config." >&2
        exit 1
    fi
    echo "Created $DIR"
elif [ ! -d "$DIR" ]; then
    echo "ERROR: $DIR exists but is not a folder." >&2
    exit 1
fi

# The folders may already be locked from an earlier run: the owner gets write access back
chmod u+w "$DIR" 2>/dev/null || true
mkdir -p "$DIR/nf-core" "$DIR/local" "$DIR/course" 2>/dev/null || true
chmod u+w "$DIR/nf-core" "$DIR/local" "$DIR/course" 2>/dev/null || true
for d in "$DIR" "$DIR/nf-core" "$DIR/local" "$DIR/course"; do
    if [ ! -d "$d" ] || [ ! -w "$d" ]; then
        owner=$(owner_of "$d")
        if [ "$owner" != "$(id -un)" ]; then
            echo "ERROR: you cannot write to $d: it belongs to '$owner' and you are '$(id -un)'." >&2
            echo "       Run this script as '$owner', or use --dir with a folder of your own." >&2
        else
            echo "ERROR: you own $d but cannot make it writable (read-only file system or immutable flag?)." >&2
        fi
        echo "       Use --check to see which images are present or missing (it changes nothing)." >&2
        exit 1
    fi
done

singularity_cmd() {
    if command -v singularity >/dev/null; then echo singularity
    elif command -v apptainer >/dev/null; then echo apptainer
    else echo "ERROR: neither singularity nor apptainer is on PATH" >&2; return 1; fi
}

while IFS=$'\t' read -r kind src target procs; do
    if [ -s "$target" ]; then
        echo "[skip]     ${target#$DIR/} (already present)"
        continue
    fi
    tmp="$target.partial"
    case "$kind" in
        download)
            echo "[download] ${target#$DIR/}  <- $src"
            curl -fL --retry 3 --retry-delay 10 -o "$tmp" "$src"
            ;;
        pull)
            echo "[pull]     ${target#$DIR/}  <- $src"
            "$(singularity_cmd)" pull --force "$tmp" "$src"
            ;;
        build)
            echo "[build]    ${target#$DIR/}  <- ${src#$PIPELINE_DIR/}  ($procs)"
            def="$WORK/$(basename "$target" .sif).def"
            cat > "$def" <<DEF
Bootstrap: docker
From: ${MICROMAMBA_IMAGE#docker://}

%files
    $src /opt/environment.yml

%post
    export MAMBA_ROOT_PREFIX=/opt/conda
    micromamba install -y -n base -f /opt/environment.yml
    micromamba install -y -n base -c conda-forge procps-ng
    micromamba clean -a -y

%environment
    export MAMBA_ROOT_PREFIX=/opt/conda
    export PATH=/opt/conda/bin:\$PATH

%labels
    neonatal-amr-nf.module ${src#$PIPELINE_DIR/}
    neonatal-amr-nf.processes $procs
DEF
            $BUILD_CMD "$tmp" "$def"
            ;;
    esac
    mv "$tmp" "$target"
done < "$WORK/plan.tsv"

# ---------------------------------------------------------------------------------------------
# 4. Manifest
# ---------------------------------------------------------------------------------------------
manifest="$DIR/MANIFEST.tsv"
chmod u+w "$manifest" 2>/dev/null || true
{
    printf "image\tsource\tsize_bytes\tsha256\tprocesses\n"
    while IFS=$'\t' read -r kind src target procs; do
        [ "$kind" = build ] && src="built from ${src#$PIPELINE_DIR/}"
        printf "%s\t%s\t%s\t%s\t%s\n" "${target#$DIR/}" "$src" "$(wc -c < "$target" | tr -d ' ')" \
            "$( (sha256sum "$target" 2>/dev/null || shasum -a 256 "$target") | cut -d' ' -f1)" "$procs"
    done < "$WORK/plan.tsv"
} > "$manifest"
echo "Manifest written: $manifest"

report || { echo "ERROR: some images are still missing" >&2; exit 1; }

# ---------------------------------------------------------------------------------------------
# 5. Lock: read-only for everyone (or one group); no write permission for anyone
# ---------------------------------------------------------------------------------------------
if [ "$LOCK" -eq 1 ]; then
    if [ -n "$GROUP" ]; then
        chgrp -R "$GROUP" "$DIR"
        find "$DIR" -type d -exec chmod 0550 {} +
        find "$DIR" -type f -exec chmod 0440 {} +
        echo "Locked: readable by owner and group '$GROUP' only; nobody can write."
    else
        find "$DIR" -type d -exec chmod 0555 {} +
        find "$DIR" -type f -exec chmod 0444 {} +
        echo "Locked: readable by everyone; nobody can write."
    fi
    echo "To add or update images later, run this script again as the owner of $DIR."
fi

# ---------------------------------------------------------------------------------------------
# 6. Participants must be able to pass through every folder above $DIR (execute permission).
#    Folders you own get 'x' for others: they can pass through but not list or change them.
# ---------------------------------------------------------------------------------------------
blocked=""
p=$(dirname "$DIR")
while [ "$p" != "/" ]; do
    perms=$(stat -c %A "$p" 2>/dev/null || stat -f %Sp "$p")
    if [ "${perms:9:1}" != "x" ] && [ "${perms:9:1}" != "t" ]; then
        if [ "$(owner_of "$p")" = "$(id -un)" ] && [ "$LOCK" -eq 1 ]; then
            chmod o+x "$p" && echo "Allowed others to pass through (not list) $p"
        else
            blocked="$blocked $p"
        fi
    fi
    p=$(dirname "$p")
done
if [ -n "$blocked" ]; then
    echo "WARNING: participants cannot pass through:$blocked" >&2
    echo "         They will not reach the images until the owner runs  chmod o+x  on these folders." >&2
fi
