# tools.sh - the course's software, from the shared Singularity images or from modules.
#
# Loaded from ~/.bashrc (Day 1, section 0.1):
#     source $SHARED/neonatal-amr-nf/course/tools.sh
#
# Then, before using a group of tools:
#     use_tools qc          # seqkit, FastQC, fastp, NanoPlot, Filtlong, Kraken2, Bracken, MultiQC
#     use_tools assembly    # Shovill, Autocycler, Flye, Raven, BWA, Polypolish, Pypolca, QUAST, CheckM2, Bandage
#     use_tools medaka      # Medaka
#     use_tools typing      # Kleborate, Kaptive, ECTyper, EzClermont
#     use_tools amr         # AMRFinderPlus
#     use_tools plasmids    # MOB-suite, Mash, Bakta (bakta, bakta_db)
#     use_tools phylo       # Snippy, SKA2, Gubbins, snp-sites, snp-dists, IQ-TREE, seqkit
#     use_tools r           # R (module) with your own package library
#     use_tools nextflow    # Nextflow (module)
#     use_tools all         # every group above at once
#     list_tools            # every group, its tools and where each comes from
#
# Forgot `use_tools`, or logged in again after a disconnect? A course tool typed on its own still runs:
# when bash cannot find a command, command_not_found_handle (below) runs it from its image (or loads its
# module) for you. `use_tools` is still useful: it checks the images up front and makes `type <tool>` work.
# Switch the automatic loading off with: export TOOLS_AUTO=0
#
# Two ways to get the tools (TOOL_MODE, set before the `source` line, or in ~/.bashrc):
#   singularity (default) : each command runs inside the pipeline's own image, so you use exactly the
#                           versions the pipeline uses. You type the commands as usual: `fastp ...` runs
#                           `singularity exec <fastp image> fastp ...` for you.
#   module                : `use_tools qc` runs `module load` for the group's modules (see GROUP_MODULES
#                           below; your trainer sets the names used on your HPC).
#
# This file is read by an interactive shell: it never exits it and never sets -e.

# ------------------------------------------------------------------------------ settings
TOOL_MODE="${TOOL_MODE:-singularity}"
CONTAINERS="${CONTAINERS:-${SHARED:-/var/scratch/global/jjuma/neonatal-amr/shared}/containers}"
TOOL_BINDS="${TOOL_BINDS:-/var/scratch,/export/data,${HOME}}"   # folders the images must see (existing ones only)
R_MODULE="${R_MODULE:-R/4.5}"
NEXTFLOW_MODULE="${NEXTFLOW_MODULE:-nextflow/25.04}"
SINGULARITY_MODULE="${SINGULARITY_MODULE:-}"                    # e.g. singularity/3.11, if `singularity` is not on PATH

# Groups and their commands (the same groups as the course's software table)
declare -gA GROUP_TOOLS=(
    [qc]="seqkit fastqc fastp NanoPlot filtlong kraken2 bracken multiqc"
    [assembly]="shovill autocycler flye raven bwa polypolish pypolca quast.py checkm2 Bandage seqkit"
    [medaka]="medaka_consensus medaka"
    [typing]="kleborate kaptive ectyper ezclermont"
    [amr]="amrfinder amrfinder_update"
    [plasmids]="mob_recon mob_typer mash bakta bakta_db"
    [phylo]="snippy snippy-core snippy-clean_full_aln ska run_gubbins.py snp-sites snp-dists iqtree seqkit"
)

# Which image each command runs in: the pipeline process that uses it (looked up in MANIFEST.tsv),
# or COURSE_<TOOL> for tools only the course uses (built from course/envs/<tool>/environment.yml)
declare -gA TOOL_IMAGE=(
    [fastp]=FASTP [NanoPlot]=NANOPLOT [filtlong]=FILTLONG [kraken2]=KRAKEN2_KRAKEN2 [bracken]=BRACKEN_BRACKEN
    [multiqc]=MULTIQC [shovill]=SHOVILL [autocycler]=AUTOCYCLER_ASSEMBLE [flye]=AUTOCYCLER_ASSEMBLE
    [raven]=AUTOCYCLER_ASSEMBLE [bwa]=POLYPOLISH_SHORTREAD [polypolish]=POLYPOLISH_SHORTREAD
    [pypolca]=PYPOLCA_RUN [quast.py]=QUAST [medaka_consensus]=MEDAKA [medaka]=MEDAKA
    [kleborate]=KLEBORATE [kaptive]=KLEBORATE [ectyper]=ECTYPER [ezclermont]=EZCLERMONT
    [amrfinder]=AMRFINDERPLUS_RUN [amrfinder_update]=AMRFINDERPLUS_RUN
    [mob_recon]=MOBSUITE_RECON [mob_typer]=MOBSUITE_RECON [mash]=MOBSUITE_RECON
    [snippy]=SNIPPY_RUN [snippy-core]=SNIPPY_RUN [snippy-clean_full_aln]=SNIPPY_RUN [snp-sites]=SNIPPY_RUN
    [run_gubbins.py]=GUBBINS [snp-dists]=SNPDISTS [iqtree]=IQTREE [ska]=SKA2 [bakta]=BAKTA_BAKTA [bakta_db]=BAKTA_BAKTA
    [seqkit]=COURSE_SEQKIT [fastqc]=COURSE_FASTQC [checkm2]=COURSE_CHECKM2 [Bandage]=COURSE_BANDAGE
)

# TOOL_MODE=module: the modules each group loads. Trainer: change these to the names on your HPC
# (`module avail`), keeping the pipeline's versions where there is one.
declare -gA GROUP_MODULES=(
    [qc]="seqkit/2.8.2 fastqc/0.12.1 fastp/1.3.6 nanoplot/1.47.0 filtlong/0.2.1 kraken2/2.1.6 bracken/3.1 multiqc/1.35"
    [assembly]="shovill/1.1.0 autocycler/0.6.2 flye/2.9.5 raven/1.6.1 bwa/0.7.19 polypolish/0.6.0 pypolca/0.4.0 quast/5.3.0 checkm2/1.0.2 bandage/0.8.1 seqkit/2.8.2"
    [medaka]="medaka/2.2.2"
    [typing]="kleborate/3.2.4 ectyper/1.0.0 ezclermont/1.0.0"
    [amr]="amrfinderplus/4.2.7"
    [plasmids]="mob_suite/3.1.9 mash/2.3 bakta/1.12.0"
    [phylo]="snippy/4.6.0 ska2/0.5.1 gubbins/3.3.5 snp-sites/2.5.1 snp-dists/0.8.2 iqtree/3.1.3 seqkit/2.8.2"
)

# Processes that share another process's image (an older MANIFEST.tsv may list only the first)
declare -gA _SAME_IMAGE_AS=( [KAPTIVE]=KLEBORATE [SNIPPY_CORE]=SNIPPY_RUN [AMRFINDERPLUS_UPDATE]=AMRFINDERPLUS_RUN )

export TOOL_MODE CONTAINERS TOOL_BINDS
# This file's path, so that a child shell (e.g. `bash my_script.sh`) can load the tables above again
_TOOLS_SH=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")
export _TOOLS_SH

# ------------------------------------------------------------------------------ helpers
_tools_msg() { printf '%s\n' "$*" >&2; }

_tools_module() {   # module load, quietly tolerating module scripts that use unset variables
    if ! type module >/dev/null 2>&1 && [ -f /etc/profile.d/modules.sh ]; then source /etc/profile.d/modules.sh; fi
    type module >/dev/null 2>&1 || { _tools_msg "ERROR: the 'module' command is not available"; return 1; }
    module load "$@"
}

_tools_singularity() {
    # The path is checked on every machine: a path found on the login node (exported into an interactive
    # session) may not exist on the compute node.
    if [ -n "${_TOOLS_SIF_BIN:-}" ] && [ -x "${_TOOLS_SIF_BIN}" ]; then return 0; fi
    _TOOLS_SIF_BIN=$(command -v singularity || command -v apptainer || true)
    if [ -z "${_TOOLS_SIF_BIN}" ]; then
        local m
        for m in ${SINGULARITY_MODULE} singularity apptainer; do
            _tools_module "$m" >/dev/null 2>&1 || continue
            _TOOLS_SIF_BIN=$(command -v singularity || command -v apptainer || true)
            [ -n "${_TOOLS_SIF_BIN}" ] && break
        done
    fi
    if [ -z "${_TOOLS_SIF_BIN}" ]; then
        _tools_msg "ERROR: singularity/apptainer not found on $(hostname): load its module (set SINGULARITY_MODULE) or use TOOL_MODE=module"
        return 1
    fi
    export _TOOLS_SIF_BIN
}

_tools_image() {    # _tools_image PROCESS -> image path (cached per shell)
    local process=$1 rel="" key
    key="_TOOLS_IMG_${process}"
    if [ -n "${!key:-}" ]; then echo "${!key}"; return 0; fi
    [ -r "${CONTAINERS}/MANIFEST.tsv" ] || { _tools_msg "ERROR: ${CONTAINERS}/MANIFEST.tsv not found (is SHARED set? ask your trainer)"; return 1; }
    for p in "$process" "${_SAME_IMAGE_AS[$process]:-}"; do
        [ -n "$p" ] || continue
        rel=$(awk -F'\t' -v p="$p" 'NR > 1 { n = split($5, procs, ",")
                 for (i = 1; i <= n; i++) { k = split(procs[i], parts, ":"); if (parts[k] == p) { print $1; exit } } }' \
                 "${CONTAINERS}/MANIFEST.tsv")
        [ -n "$rel" ] && break
    done
    [ -n "$rel" ] || { _tools_msg "ERROR: no image for ${process} in ${CONTAINERS}/MANIFEST.tsv (ask your trainer)"; return 1; }
    [ -r "${CONTAINERS}/${rel}" ] || { _tools_msg "ERROR: image not readable: ${CONTAINERS}/${rel}"; return 1; }
    printf -v "$key" '%s' "${CONTAINERS}/${rel}"; export "${key?}"
    echo "${CONTAINERS}/${rel}"
}

# _tools_exec PROCESS COMMAND [args...] : run COMMAND inside PROCESS's image
_tools_exec() {
    local process=$1 tool=$2 image binds="" d
    local -a dirs
    shift 2
    _tools_singularity || return 127
    image=$(_tools_image "$process") || return 127
    IFS=',' read -ra dirs <<< "${TOOL_BINDS}"
    for d in "${dirs[@]}"; do [ -d "$d" ] && binds="${binds} -B ${d}"; done
    (
        # The wrapper functions (fastp, medaka, ...) are exported so that they also work in child shells.
        # Keep them out of the container: there, a tool's own bash scripts (medaka_consensus calls 'medaka')
        # would run the wrapper instead of the tool. Only this subshell is changed, not your shell.
        while read -r _ _ f; do export -fn "$f" 2>/dev/null; done < <(declare -Fx)
        if [ "$tool" = Bandage ]; then        # Bandage draws images without a display
            export SINGULARITYENV_QT_QPA_PLATFORM=offscreen APPTAINERENV_QT_QPA_PLATFORM=offscreen
            # Qt needs a private runtime folder. The login node's /run/user/<uid> is copied into an interactive
            # session but does not exist on the compute node, so use one in /tmp (shared with the container).
            if [ ! -d "${XDG_RUNTIME_DIR:-/nonexistent}" ] || [ ! -w "${XDG_RUNTIME_DIR:-/nonexistent}" ]; then
                export XDG_RUNTIME_DIR="${TMPDIR:-/tmp}/runtime-${USER:-$(id -un)}"
                mkdir -p -m 700 "$XDG_RUNTIME_DIR"
            fi
            export SINGULARITYENV_XDG_RUNTIME_DIR="$XDG_RUNTIME_DIR" APPTAINERENV_XDG_RUNTIME_DIR="$XDG_RUNTIME_DIR"
        fi
        # The images have no en_US.UTF-8 locale: without this, Perl tools (Prokka, Snippy, ...) warn
        # "Setting locale failed". (A "complete: unsupported builtin" line from some images is harmless.)
        export LC_ALL=C
        exec "${_TOOLS_SIF_BIN}" exec ${binds} "$image" "$tool" "$@"
    )
}
export -f _tools_msg _tools_module _tools_singularity _tools_image _tools_exec 2>/dev/null

_tools_group_of() {  # _tools_group_of TOOL -> the first group that lists it
    local g
    for g in qc assembly medaka typing amr plasmids phylo; do
        case " ${GROUP_TOOLS[$g]} " in *" $1 "*) echo "$g"; return 0 ;; esac
    done
    return 1
}
export -f _tools_group_of 2>/dev/null

# ------------------------------------------------------------------------------ what you type
use_tools() {
    local group=${1:-} t
    if [ -z "$group" ]; then _tools_msg "usage: use_tools qc|assembly|medaka|typing|amr|plasmids|phylo|all|r|nextflow"; return 1; fi
    case "$group" in
        r|R)
            _tools_module "${R_MODULE}" || return 1
            # settings from another R break its start-up ('methods' errors): use only your own library
            unset R_HOME R_LIBS R_LIBS_SITE R_ENVIRON R_ENVIRON_USER R_PROFILE R_PROFILE_USER
            _tools_msg "R: $(command -v R) ($R_MODULE), with packages from your own library in ~/R/"
            return 0 ;;
        nextflow)
            _tools_module "${NEXTFLOW_MODULE}" || return 1
            _tools_msg "Nextflow: $(command -v nextflow) ($NEXTFLOW_MODULE)"
            return 0 ;;
    esac
    if [ "$group" = all ]; then
        local g status=0
        for g in qc assembly medaka typing amr plasmids phylo; do use_tools "$g" || status=1; done
        return $status
    fi
    [ -n "${GROUP_TOOLS[$group]:-}" ] || { _tools_msg "unknown group '$group': use one of ${!GROUP_TOOLS[*]} all r nextflow"; return 1; }

    if [ "${TOOL_MODE}" = module ]; then
        _tools_module ${GROUP_MODULES[$group]} || return 1
        _tools_msg "Loaded modules for '$group': ${GROUP_MODULES[$group]}"
    else
        _tools_singularity || return 1
        local ready="" missing=""
        for t in ${GROUP_TOOLS[$group]}; do
            if _tools_image "${TOOL_IMAGE[$t]}" >/dev/null 2>&1; then
                eval "${t}() { _tools_exec ${TOOL_IMAGE[$t]} ${t} \"\$@\"; }"
                export -f "$t" 2>/dev/null
                ready="$ready $t"
            else
                missing="$missing $t"
            fi
        done
        _tools_msg "Tools for '$group' (from images in ${CONTAINERS}):${ready:- none}"
        if [ -n "$missing" ]; then
            _tools_msg "WARNING: no image for:$missing (ask your trainer; \`list_tools\` shows what is missing)"
            return 1
        fi
    fi
}

# Which image (or module) gives each tool
list_tools() {
    local g t
    for g in qc assembly medaka typing amr plasmids phylo; do
        echo "== $g"
        if [ "${TOOL_MODE}" = module ]; then
            echo "   modules: ${GROUP_MODULES[$g]}"
        else
            for t in ${GROUP_TOOLS[$g]}; do
                printf '   %-22s %s\n' "$t" "$(_tools_image "${TOOL_IMAGE[$t]}" 2>&1 | sed "s#^${CONTAINERS}/##")"
            done
        fi
    done
    echo "== r          module ${R_MODULE}"
    echo "== nextflow   module ${NEXTFLOW_MODULE}"
}

# ------------------------------------------------------------------------------ automatic loading
# bash calls command_not_found_handle when a command is not found. For a course tool, run it from its image
# (TOOL_MODE=singularity) or load its group's modules first (TOOL_MODE=module), so that the course's commands
# work without `use_tools`, e.g. in a new session after a disconnect. Anything else is passed on to the
# system's own handler, if there is one.
# (Only a handler that is not this one is kept, so sourcing this file again does not chain it to itself.)
if declare -F command_not_found_handle >/dev/null 2>&1 && ! declare -F _tools_prev_cnf >/dev/null 2>&1 \
        && ! declare -f command_not_found_handle | grep -q _TOOLS_SH; then
    eval "$(declare -f command_not_found_handle | sed '1s/^command_not_found_handle/_tools_prev_cnf/')"
fi
command_not_found_handle() {
    local cmd=$1 g
    # In a child shell the tables are not inherited (bash does not export arrays): read this file again
    if ! declare -p TOOL_IMAGE >/dev/null 2>&1 && [ -r "${_TOOLS_SH:-}" ]; then
        source "${_TOOLS_SH}" >/dev/null 2>&1
    fi
    # nextflow, R and Rscript come from modules in both modes
    if [ "${TOOLS_AUTO:-1}" = 1 ] && declare -F use_tools >/dev/null 2>&1; then
        case "$cmd" in
            nextflow)  use_tools nextflow >/dev/null 2>&1 && hash -r && type -P "$cmd" >/dev/null && { "$@"; return $?; } ;;
            R|Rscript) use_tools r        >/dev/null 2>&1 && hash -r && type -P "$cmd" >/dev/null && { "$@"; return $?; } ;;
        esac
    fi
    if [ "${TOOLS_AUTO:-1}" = 1 ] && declare -p TOOL_IMAGE >/dev/null 2>&1 && [ -n "${TOOL_IMAGE[$cmd]:-}" ]; then
        if [ "${TOOL_MODE}" = module ]; then
            g=$(_tools_group_of "$cmd") && _tools_module ${GROUP_MODULES[$g]} >/dev/null 2>&1 \
                && hash -r && type -P "$cmd" >/dev/null && { "$@"; return $?; }
        else
            _tools_exec "${TOOL_IMAGE[$cmd]}" "$@"; return $?
        fi
    fi
    if declare -F _tools_prev_cnf >/dev/null 2>&1; then _tools_prev_cnf "$@"; return $?; fi
    printf 'bash: %s: command not found\n' "$cmd" >&2
    return 127
}
export -f command_not_found_handle 2>/dev/null
