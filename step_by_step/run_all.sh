#!/usr/bin/env bash
#
# run_all.sh - every step, one after the other, for every sample in the samplesheet.
# This is what the Nextflow pipeline automates (in parallel, with resume and retries).
#
#   mkdir -p logs && sbatch step_by_step/run_all.sh          # as a Slurm job
#   bash step_by_step/run_all.sh                             # in an interactive session
#   SKIP="05 10 10b" bash step_by_step/run_all.sh            # leave out steps (05 hybrid is the slowest)
#   SKIP=05 PHYLO_METHOD=both bash step_by_step/run_all.sh   # with Prokka and Bakta; Snippy and SKA2 in step 11
#
#SBATCH --job-name=amr-steps
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/%u_steps_%j.txt
#SBATCH -w compute05

set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
[ -f "$HERE/config.sh" ] || HERE="${SLURM_SUBMIT_DIR:-$PWD}/step_by_step"   # sbatch runs a copy of the script
SKIP=" ${SKIP:-10 10b} "      # Prokka (10) and Bakta (10b) are optional, as in the pipeline

source "$HERE/config.sh"
reads=$(is_assembly_sheet && echo no || echo yes)

# 00_, 01_, ... 10_, 10b_, 11_, 12_ in that order
for script in $(printf '%s\n' "$HERE"/[0-9][0-9]*_*.sh | LC_ALL=C sort); do
    n=$(basename "$script" | cut -d_ -f1)
    [[ "$SKIP" == *" $n "* ]] && { info "step $n skipped"; continue; }
    [ "$reads" = no ] && [[ " 01 02 03 04 05 " == *" $n "* ]] && continue    # assemblies: no read steps
    case "$n" in
        00|11|12) bash "$script" ;;         # whole-run steps
        *)        bash "$script" all ;;     # per-sample steps
    esac
done
step "Done: ${OUTDIR}"
