#!/usr/bin/env bash
#
# config.sh - settings and helper functions shared by the step-by-step scripts.
# Every step script starts with:   source "$(dirname "$0")/config.sh"
#
# Edit the SETTINGS block (or set any of these variables before running a step, e.g.
#     SAMPLESHEET=/path/sheet.csv bash step_by_step/01_fastp.sh all
# ). Nothing below the SETTINGS block needs changing.

# ============================================================================== SETTINGS
PIPELINE_DIR="${PIPELINE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"

# Samples: the same samplesheet the pipeline uses (reads: fastq_1/fastq_2[/long_fastq]; assemblies: fasta)
SAMPLESHEET="${SAMPLESHEET:-${PIPELINE_DIR}/assets/samplesheet_test_5samples.csv}"

# Where the results go. The folder layout is the same as the pipeline's --outdir, so you can compare them.
OUTDIR="${OUTDIR:-/var/scratch/global/${USER}/neonatal-amr/step_by_step}"

# Shared, read-only resources prepared by the trainer (the same ones the pipeline uses)
CONTAINERS="${CONTAINERS:-/var/scratch/global/jjuma/neonatal-amr/shared/containers}"
REFERENCES="${REFERENCES:-/var/scratch/global/jjuma/neonatal-amr/shared/references}"
KRAKEN2_DB="${KRAKEN2_DB:-/export/data/bio/kraken2/db/standard_20241228}"
AMRFINDER_DB="${AMRFINDER_DB:-}"        # .tar.gz or folder; empty = download it in step 08 (needs internet)
REFERENCE_SHEET="${REFERENCE_SHEET:-${PIPELINE_DIR}/conf/phylo_references.csv}"
MIN_GROUP="${MIN_GROUP:-2}"             # phylogeny: smallest group (species + ST) that is analysed

# How the tools are run
#   singularity : from the pipeline's shared images (same versions as the pipeline; nothing to install)
#   module      : from whatever is on your PATH, e.g. after `module load fastp/1.3.6 ...` (see README)
TOOL_MODE="${TOOL_MODE:-singularity}"
BINDS="${BINDS:-/var/scratch,/export/data,${HOME}}"   # folders the containers must see (only existing ones are used)

# Resources: the Slurm allocation when there is one (1 CPU in a course session; run_all.sh as a batch job gets 8),
# otherwise $THREADS (1 in the course's ~/.bashrc)
THREADS="${SLURM_CPUS_PER_TASK:-${THREADS:-1}}"
MEM_GB="${MEM_GB:-$(( ${SLURM_MEM_PER_NODE:-16384} / 1024 ))}"
# ======================================================================================

set -euo pipefail

# --------------------------------------------------------------------------- messages
step()  { printf '\n\033[1;34m==== %s\033[0m\n' "$*"; }
info()  { printf '     %s\n' "$*"; }
warn()  { printf '\033[1;33mWARNING:\033[0m %s\n' "$*" >&2; }
fail()  { printf '\033[1;31mERROR:\033[0m %s\n' "$*" >&2; exit 1; }

# ------------------------------------------------------------------- container images
# The image of a pipeline process (FASTP, SHOVILL, ...), looked up in the shared folder's MANIFEST.tsv
# (written by bin/prepare_shared_containers.sh). Override one with IMG_<PROCESS>=/path/image.sif.
# Processes that share another process's image are looked up under that process too, so an older
# MANIFEST.tsv (written before these processes were added) still works.
declare -A SAME_IMAGE_AS=(
    [KAPTIVE]=KLEBORATE                       # Kaptive is installed with Kleborate
    [SNIPPY_CONTIGS]=SNIPPY_RUN               # all Snippy steps use the Snippy 4.6.0 image
    [SNIPPY_CLEAN]=SNIPPY_RUN
    [SNIPPY_CORE]=SNIPPY_RUN
    [AMRFINDERPLUS_UPDATE]=AMRFINDERPLUS_RUN
)
manifest_image() {   # manifest_image PROCESS -> image path relative to CONTAINERS, or nothing
    awk -F'\t' -v p="$1" 'NR > 1 {
            n = split($5, procs, ",")
            for (i = 1; i <= n; i++) { k = split(procs[i], parts, ":"); if (parts[k] == p) { print $1; exit } }
        }' "${CONTAINERS}/MANIFEST.tsv"
}
image_for() {
    local process=$1 override="IMG_$1" image
    if [ -n "${!override:-}" ]; then echo "${!override}"; return; fi
    [ -r "${CONTAINERS}/MANIFEST.tsv" ] || fail "no ${CONTAINERS}/MANIFEST.tsv; set CONTAINERS, or use TOOL_MODE=module"
    image=$(manifest_image "$process")
    if [ -z "$image" ] && [ -n "${SAME_IMAGE_AS[$process]:-}" ]; then
        image=$(manifest_image "${SAME_IMAGE_AS[$process]}")
    fi
    [ -n "$image" ] || fail "no image for process ${process} in ${CONTAINERS}/MANIFEST.tsv (set IMG_${process}=/path/image.sif)"
    echo "${CONTAINERS}/${image}"
}

SINGULARITY_BIN=$(command -v singularity || command -v apptainer || true)
bind_opts() {
    local opts="" d
    IFS=',' read -ra dirs <<< "${BINDS}"
    for d in "${dirs[@]}"; do [ -d "$d" ] && opts="${opts} -B ${d}"; done
    echo "$opts"
}

# run PROCESS command [args...]
#   Prints the command, then runs it: inside PROCESS's image (TOOL_MODE=singularity) or directly
#   (TOOL_MODE=module). PROCESS is the pipeline's process name, so you can find the same step in the pipeline.
in_tool_env() {
    local process=$1; shift
    if [ "${TOOL_MODE}" = singularity ]; then
        [ -n "${SINGULARITY_BIN}" ] || fail "singularity/apptainer not found; load its module or use TOOL_MODE=module"
        local image; image=$(image_for "$process") || exit 1
        [ -r "$image" ] || fail "image not found: $image"
        # shellcheck disable=SC2046
        "${SINGULARITY_BIN}" exec $(bind_opts) "$image" "$@"
    else
        "$@"
    fi
}
run() {
    local process=$1; shift
    printf '\033[0;32m$ %s\033[0m\n' "$*"
    in_tool_env "$process" "$@"
}

# run_sh PROCESS 'shell command line'   - for commands with pipes, loops or redirections
run_sh() {
    printf '\033[0;32m$ %s\033[0m\n' "$(printf '%s' "$2" | sed -E 's/([^ ])  +/\1 /g; s/^\n//')"
    in_tool_env "$1" bash -c "set -euo pipefail; $2"
}

# --------------------------------------------------------------------------- samplesheet
[ -r "${SAMPLESHEET}" ] || fail "samplesheet not found: ${SAMPLESHEET}"
SHEET_DIR=$(cd "$(dirname "${SAMPLESHEET}")" && pwd)

# field SAMPLE COLUMN  -> value of COLUMN for SAMPLE ('' if the column or value is missing)
field() {
    awk -F',' -v s="$1" -v c="$2" 'NR == 1 { for (i = 1; i <= NF; i++) col[$i] = i; next }
                                   $col["sample"] == s { if (c in col) print $col[c]; exit }' "${SAMPLESHEET}" | tr -d '\r'
}
# path SAMPLE COLUMN -> like field, with relative paths made absolute (relative to the samplesheet)
path() {
    local p; p=$(field "$1" "$2")
    [ -z "$p" ] || [[ "$p" = /* ]] || p="${SHEET_DIR}/${p}"
    echo "$p"
}
all_samples() { awk -F',' 'NR == 1 { for (i = 1; i <= NF; i++) col[$i] = i; next } { print $col["sample"] }' "${SAMPLESHEET}" | tr -d '\r'; }
is_assembly_sheet() { head -n 1 "${SAMPLESHEET}" | tr ',' '\n' | grep -x fasta > /dev/null; }

# Samples to process: the script's arguments ('all' = every sample in the samplesheet; none = the first one)
samples_from_args() {
    if [ $# -eq 0 ]; then all_samples | awk 'NR == 1'
    elif [ "$1" = all ]; then all_samples
    else printf '%s\n' "$@"
    fi
}

# ------------------------------------------------------------------- results of earlier steps
trimmed_r1()   { echo "${OUTDIR}/01_qc/fastp/$1_R1.fastp.fastq.gz"; }
trimmed_r2()   { echo "${OUTDIR}/01_qc/fastp/$1_R2.fastp.fastq.gz"; }
filtered_ont() { echo "${OUTDIR}/01_qc/filtlong/$1.filtlong.fastq.gz"; }

# best_assembly SAMPLE -> the assembly the pipeline sends downstream:
#   the samplesheet's fasta (assembly input), else the hybrid assembly if there is one, else Shovill's contigs
best_assembly() {
    if is_assembly_sheet; then path "$1" fasta; return; fi
    local hybrid="${OUTDIR}/02_assembly/hybrid/$1.hybrid.fasta" short="${OUTDIR}/02_assembly/shovill/$1.contigs.fa"
    if [ -s "$hybrid" ]; then echo "$hybrid"; elif [ -s "$short" ]; then echo "$short"; else echo ""; fi
}
need_file() { [ -s "$1" ] || fail "$1 not found or empty: ${2:-run the previous step first}"; }

mkdir -p "${OUTDIR}"
