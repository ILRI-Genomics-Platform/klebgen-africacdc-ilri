#!/bin/bash
#
# plot_results.sbatch - after the pipeline: transmission pairs, figures and (optionally) BactDating.
# Edit the I/O block, then from the pipeline folder either
#     bash bin/plot_results.sbatch                      # on the head node, where R/4.5 is installed
# or, if the compute nodes have the same R module,
#     mkdir -p logs && sbatch bin/plot_results.sbatch
# Install the R packages first, once, on the login node (it has internet; compute nodes may not):
#     bin/run_R.sh install_r_packages.R

#SBATCH -w compute06
#SBATCH --job-name=plot-neonatal-amr
#SBATCH --output=logs/%u_plots_%j.out
#SBATCH --error=logs/%u_plots_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --partition=batch

set -euo pipefail

# load the required module(s); packages come from your own ~/R/x86_64-pc-linux-gnu-library/<version>
R_MODULE="${R_MODULE:-R/4.5}"
module load "${R_MODULE}"
export NEONATAL_R_NO_INSTALL=1        # do not try to download packages during the run (install them first)

# I/O
# For genomes given as assemblies, use the metadata written by bin/make_assembly_samplesheet.py
# (metadata_assemblies.tsv: sample, collection_year, collection_date, ward, ...).
PIPELINE_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
RESULTS="/var/scratch/global/${USER}/neonatal-amr/outdir"
METADATA="${PIPELINE_DIR}/training_data/metadata_compact.tsv"
AST="${PIPELINE_DIR}/training_data/ast_compact.tsv"   # leave empty ("") if there is no AST table
OUTDIR="${RESULTS}/plots"
RUN_BACTDATING=false                  # true: one BactDating run per Gubbins group (slow: minutes to hours)
MRCA=""                               # e.g. "CAA65Z,CAA99X": date these isolates' common ancestor

R="${PIPELINE_DIR}/bin/run_R.sh"
export R_MODULE=none                  # run_R.sh: the R module is already loaded above
mkdir -p "${OUTDIR}"

echo "== Transmission pairs (surveillance question 5)"
"$R" transmission_pairs.R \
    --dist "${RESULTS}/06_phylogeny" \
    --metadata "${METADATA}" \
    --out "${OUTDIR}/transmission_pairs.tsv"

echo "== Figures"
"$R" make_figures.R \
    --results "${RESULTS}" \
    --metadata "${METADATA}" \
    ${AST:+--ast "${AST}"} \
    --outdir "${OUTDIR}/figures"

if [ "${RUN_BACTDATING}" = true ]; then
    for prefix in "${RESULTS}"/06_phylogeny/gubbins/*.final_tree.tre; do
        [ -e "$prefix" ] || continue
        group=$(basename "$prefix" .final_tree.tre)
        echo "== BactDating: ${group}"
        "$R" run_bactdating.R \
            --gubbins "${RESULTS}/06_phylogeny/gubbins/${group}" \
            --metadata "${METADATA}" \
            --outdir "${OUTDIR}/bactdating_${group}" \
            ${MRCA:+--mrca "${MRCA}"} || echo "BactDating failed for ${group}; continuing"
    done
fi

echo "Done: ${OUTDIR}"
