#!/bin/bash
#
# download_terra_outputs.sbatch - download Terra outputs (assemblies, Kleborate, MLST) in one parallel pass.
# Edit the settings block, then from the pipeline folder:
#     mkdir -p logs && sbatch bin/download_terra_outputs.sbatch
# or override at submission:
#     sbatch --export=ALL,TERRA_TABLE=/path/table.tsv,OUTDIR=/path/out,BY_COUNTRY=1 bin/download_terra_outputs.sbatch
#
# Downloads are network-bound: more CPUs only help the one-file-at-a-time fallback (--jobs).
# Needs the Google Cloud CLI, logged in once on a login node:  gcloud auth login
# Rerunning is safe: files already downloaded are skipped.
#
# Then build the pipeline samplesheet from the assemblies (course Day 5):
#     python3 bin/make_assembly_samplesheet.py --metadata <metadata.tsv> --contigs-dir $OUTDIR/assemblies ...
#   (with BY_COUNTRY=1: --contigs-dir $OUTDIR/<Country>/assemblies)

#SBATCH --job-name=terra-download
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=4G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%u_terra_download_%j.txt
#SBATCH --error=logs/%u_terra_download_error_%j.txt
#SBATCH --mail-type=END,FAIL
##SBATCH --mail-user=you@example.org     # remove one '#' and set your address to get e-mails
##SBATCH -w compute05                   # optional: pin to a node with good internet access

set -euo pipefail

# ------------------------------------------------------------------ settings
PIPELINE_DIR="${PIPELINE_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
TERRA_TABLE="${TERRA_TABLE:-${HOME}/trainings/africacdc-pgi-klebgen/metadata/KlebGenNov2025.tsv}"
OUTDIR="${OUTDIR:-/var/scratch/global/${USER}/trainings/africacdc-pgi-klebgen/results/africa}"
FILES="${FILES:-assembly_fasta kleborate_klebsiella_output_file kleborate_ecoli_output_file ts_mlst_results}"
BY_COUNTRY="${BY_COUNTRY:-0}"                   # 1: <OUTDIR>/<Country>/<type>/ instead of <OUTDIR>/<type>/
COUNTRIES="${COUNTRIES:-}"                      # e.g. "Kenya Ethiopia"; empty = all
BILLING_PROJECT="${BILLING_PROJECT:-}"          # only for requester-pays buckets
GCLOUD_ENV="${GCLOUD_ENV:-${HOME}/miniconda3/envs/gcloud-env}"   # conda env with the Google Cloud CLI

# ------------------------------------------------------------------ environment
if ! command -v gcloud >/dev/null 2>&1; then
    set +u
    source "${HOME}/miniconda3/bin/activate"
    conda activate "${GCLOUD_ENV}"
    set -u
fi
command -v gcloud >/dev/null || { echo "[ERROR] gcloud not found (set GCLOUD_ENV)" >&2; exit 1; }
gcloud auth list --filter=status:ACTIVE --format="value(account)" | grep -q . \
    || { echo "[ERROR] gcloud is not logged in: run 'gcloud auth login' on a login node" >&2; exit 1; }
[[ -f "${TERRA_TABLE}" ]] || { echo "[ERROR] Terra table not found: ${TERRA_TABLE}" >&2; exit 1; }

args=(--terra-output "${TERRA_TABLE}" --outdir "${OUTDIR}" --files ${FILES} --jobs "${SLURM_CPUS_PER_TASK:-8}")
[[ "${BY_COUNTRY}" == 1 ]] && args+=(--by-country)
[[ -n "${COUNTRIES}" ]] && args+=(--countries ${COUNTRIES})
[[ -n "${BILLING_PROJECT}" ]] && args+=(--billing-project "${BILLING_PROJECT}")

echo "Start: $(date)   Node: $(hostname)   Table: ${TERRA_TABLE}"
python3 "${PIPELINE_DIR}/bin/download_terra_outputs.py" "${args[@]}"
echo "End:   $(date)"
