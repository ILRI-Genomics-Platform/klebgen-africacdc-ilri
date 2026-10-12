#!/usr/bin/env bash
#
# run_R.sh - run one of the course's R scripts with a clean, known R.
#
#   bin/run_R.sh make_figures.R --results results --metadata training_data/metadata_compact.tsv
#   bin/run_R.sh run_bactdating.R --gubbins results/06_phylogeny/gubbins/Kpneumoniae_ST39 --metadata ...
#   bin/run_R.sh transmission_pairs.R --dist results/06_phylogeny --metadata ...
#   bin/run_R.sh install_r_packages.R                 # once, on a node with internet access
#
# Why: errors such as "package 'methods' ... was not found" or "'methods' does not have a namespace" mean that
# Rscript started with another R's settings: a conda R next to the R module, or R_HOME / R_LIBS / ~/.Rprofile /
# ~/.Renviron left over from another R version. This wrapper loads one R (module R/4.5 by default; set R_MODULE
# to change, or R_MODULE=none to use the Rscript on PATH), clears those settings, and runs Rscript --vanilla.
# The scripts then use only your own library, ~/R/<platform>-library/<version> (bin/r_setup.R).

set -euo pipefail

if [ $# -lt 1 ] || [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
    sed -n '3,15p' "$0" | sed 's/^# \{0,1\}//'
    exit 0
fi

BIN_DIR=$(cd "$(dirname "$0")" && pwd)
SCRIPT="$1"; shift
[ -f "$SCRIPT" ] || SCRIPT="$BIN_DIR/$SCRIPT"
[ -f "$SCRIPT" ] || { echo "ERROR: R script not found: $(basename "$SCRIPT") (looked in . and $BIN_DIR)" >&2; exit 1; }

# A conda environment's R would shadow the module's R
if [ -n "${CONDA_PREFIX:-}" ] && [ -x "$CONDA_PREFIX/bin/Rscript" ] && [ "${R_MODULE:-R/4.5}" != "none" ]; then
    echo "NOTE: ignoring the R in the active conda environment ($CONDA_PREFIX); using module ${R_MODULE:-R/4.5}" >&2
    PATH=$(echo "$PATH" | tr ':' '\n' | grep -vx "$CONDA_PREFIX/bin" | paste -sd: -)
fi

R_MODULE=${R_MODULE:-R/4.5}
if [ "$R_MODULE" != "none" ]; then
    if ! command -v module >/dev/null 2>&1 && [ -f /etc/profile.d/modules.sh ]; then
        # shellcheck disable=SC1091
        source /etc/profile.d/modules.sh
    fi
    set +u
    module load "$R_MODULE" || { echo "ERROR: 'module load $R_MODULE' failed; check 'module avail R' and set R_MODULE" >&2; exit 1; }
    set -u
fi

# Settings from another R installation, or other libraries, break startup or mix package versions
unset R_HOME R_LIBS R_LIBS_USER R_LIBS_SITE R_ENVIRON R_ENVIRON_USER R_PROFILE R_PROFILE_USER

RSCRIPT=$(command -v Rscript) || { echo "ERROR: Rscript not found on PATH" >&2; exit 1; }
echo "Using $RSCRIPT ($("$RSCRIPT" --version 2>&1 | head -n 1))" >&2

# --vanilla: no ~/.Rprofile or ~/.Renviron (a common source of the 'methods' error)
exec "$RSCRIPT" --vanilla "$SCRIPT" "$@"
