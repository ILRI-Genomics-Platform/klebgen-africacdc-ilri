#!/bin/bash
#
# submit_neonatal_amr_assemblies.sbatch - run neonatal-amr-nf on ASSEMBLIES (contigs) on the ILRI HPC.
# For sequencing reads (FASTQ), use submit_neonatal_amr_reads.sbatch.
#
# The samplesheets come from bin/make_assembly_samplesheet.py --split-by-species (columns sample,species,fasta,phylo_group),
# one per species; each species is a separate run. Read QC, Kraken2 and assembly are skipped.
#
# This job runs only the Nextflow "head" process: it uses little CPU and memory but must stay alive for the
# whole run. Nextflow submits every analysis step as its own Slurm job (conf/slurm.config chooses the nodes).
#
# From the pipeline folder:
#     mkdir -p logs                          # Slurm does not create the log folder; the job fails without it
#     sbatch submit_neonatal_amr_assemblies.sbatch                                  # K. pneumoniae (default)
#     sbatch --export=ALL,SPECIES=ecoli submit_neonatal_amr_assemblies.sbatch       # E. coli
#
# Samplesheets in another folder:
#     sbatch --export=ALL,SPECIES=ecoli,ASSEMBLY_DIR=/path/to/samplesheets submit_neonatal_amr_assemblies.sbatch
#
# Follow it:   squeue -u $USER       tail -f logs/<user>_assemblies_<jobid>.txt
# Stop it:     scancel <jobid>       (the analysis jobs stop too; submit again and -resume continues)

#SBATCH --job-name=amr-assemblies
#SBATCH --partition=batch
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=3-00:00:00
#SBATCH --output=logs/%u_assemblies_%j.txt
#SBATCH --error=logs/%u_assemblies_error_%j.txt
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=J.Juma@cgiar.org
#SBATCH -w compute05                   # optional: pin the head job to one node

set -euo pipefail

echo "================================================="
echo "Node Hostname : $(hostname)"
echo "Job ID        : ${SLURM_JOB_ID:-none}"
echo "Start Time    : $(date)"
echo "================================================="

#-------------------------------------------------------------------------------
# Paths (each can be overridden with sbatch --export=ALL,NAME=value)
#-------------------------------------------------------------------------------
PIPELINE_DIR="${PIPELINE_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"       # folder with main.nf (where you ran sbatch)
BASE_DIR="${BASE_DIR:-/var/scratch/global/${USER}/neonatal-amr}"

# Folder where bin/make_assembly_samplesheet.py --split-by-species wrote
#   samplesheet_assemblies_kpneumoniae.csv / samplesheet_assemblies_ecoli.csv   (edit to where yours are)
ASSEMBLY_DIR="${ASSEMBLY_DIR:-${HOME}/trainings/africacdc-pgi-klebgen/results/Ethiopia/assemblies/Klebsiella_pneumoniae/}"
SPECIES="${SPECIES:-kpneumoniae}"                                 # kpneumoniae | ecoli

RUN_NAME="${RUN_NAME:-ethiopia_${SPECIES}}"
INPUT="${INPUT:-${HOME}/trainings/africacdc-pgi-klebgen/results/Ethiopia/samplesheet_assemblies_${SPECIES}.csv}"
OUTDIR="${OUTDIR:-${BASE_DIR}/outdir_${RUN_NAME}}"

WORKDIR="${WORKDIR:-${BASE_DIR}/work}"
MIN_GROUP="${MIN_GROUP:-3}"                                       # --phylo_min_group_size

# Nextflow keeps its run history (needed by -resume) in the folder it starts from.
# A fixed launch folder per run name means a resubmission always resumes the same run.
LAUNCH_DIR="${LAUNCH_DIR:-${BASE_DIR}/launch/${RUN_NAME}}"

#-------------------------------------------------------------------------------
# Environment
#-------------------------------------------------------------------------------
set +u                                        # module scripts can use unset variables
module load nextflow/25.04
set -u

export NXF_OPTS="-Xms512M -Xmx4G"
export NXF_ANSI_LOG=false                     # plain progress lines in the log file
# Singularity's own temporary files on scratch, not in your home quota.
# (The pipeline's images come from the shared read-only folder; nothing is downloaded.)
export SINGULARITY_CACHEDIR="${BASE_DIR}/.singularity_cache"
export APPTAINER_CACHEDIR="${SINGULARITY_CACHEDIR}"
export SINGULARITY_TMPDIR="${BASE_DIR}/.singularity_tmp"
export APPTAINER_TMPDIR="${SINGULARITY_TMPDIR}"

#-------------------------------------------------------------------------------
# Checks before anything is submitted
#-------------------------------------------------------------------------------
fail() { echo "[ERROR] $*" >&2; exit 1; }

[[ -f "${PIPELINE_DIR}/main.nf" ]]           || fail "main.nf not found in ${PIPELINE_DIR}; run sbatch from the pipeline folder or set PIPELINE_DIR"
[[ -f "${PIPELINE_DIR}/conf/slurm.config" ]] || fail "conf/slurm.config not found in ${PIPELINE_DIR}"
[[ -f "${INPUT}" ]]                          || fail "samplesheet not found: ${INPUT} (set ASSEMBLY_DIR or INPUT)"
head -n 1 "${INPUT}" | grep -q 'fasta' \
    || fail "${INPUT} has no fasta column: it is not an assembly samplesheet (for reads use submit_neonatal_amr_reads.sbatch)"

# Every contigs file in the samplesheet must exist and be readable from the compute nodes
INPUT_DIR=$(cd "$(dirname "${INPUT}")" && pwd)
missing=0
while IFS= read -r f; do
    [[ -z "$f" ]] && continue
    [[ "$f" = /* ]] || f="${INPUT_DIR}/${f}"
    [[ -r "$f" ]] || { echo "[ERROR] contigs file missing or unreadable: $f" >&2; missing=$((missing + 1)); }
done < <(awk -F',' 'NR==1 { for (i = 1; i <= NF; i++) col[$i] = i; next } { print $col["fasta"] }' "${INPUT}")
(( missing == 0 )) || fail "${missing} contigs file(s) missing"

INPUT="${INPUT_DIR}/$(basename "${INPUT}")"
mkdir -p "${LAUNCH_DIR}" "${OUTDIR}" "${WORKDIR}" "${SINGULARITY_CACHEDIR}" "${SINGULARITY_TMPDIR}"

echo " Input      : assemblies (${SPECIES})"
echo " Pipeline   : ${PIPELINE_DIR}"
echo " Samplesheet: ${INPUT} ($(($(wc -l < "${INPUT}") - 1)) samples)"
echo " Outdir     : ${OUTDIR}"
echo " Work dir   : ${WORKDIR}"
echo " Launch dir : ${LAUNCH_DIR}"
echo "================================================="

#-------------------------------------------------------------------------------
# Run
#-------------------------------------------------------------------------------
cd "${LAUNCH_DIR}"

status=0
nextflow run "${PIPELINE_DIR}/main.nf" \
    -profile singularity \
    -c "${PIPELINE_DIR}/conf/slurm.config" \
    --input "${INPUT}" \
    --input_type assembly \
    --outdir "${OUTDIR}" \
    -work-dir "${WORKDIR}" \
    --phylo_min_group_size "${MIN_GROUP}" \
    -resume \
    || status=$?

#-------------------------------------------------------------------------------
# Summary
#-------------------------------------------------------------------------------
echo "================================================="
echo " Nextflow exit status: ${status}"
echo " End Time            : $(date)"
echo " Nextflow log        : ${LAUNCH_DIR}/.nextflow.log"

trace=$(ls -t "${OUTDIR}"/pipeline_info/execution_trace_*.txt 2>/dev/null | head -n 1 || true)
if [[ -n "${trace}" ]]; then
    echo " Trace               : ${trace}"
    echo "-------------------------------------------------"
    echo " Tasks by status:"
    awk -F'\t' 'NR==1 { for (i = 1; i <= NF; i++) if ($i == "status") s = i; next } { n[$s]++ }
                END { for (k in n) printf "   %-10s %d\n", k, n[k] }' "${trace}"
    failed=$(awk -F'\t' -v wd="${WORKDIR}" 'NR==1 { for (i = 1; i <= NF; i++) { if ($i == "status") s = i; if ($i == "name") m = i; if ($i == "hash") h = i }; next }
                         $s == "FAILED" { print "   " $m "  ->  " wd "/" $h "*" }' "${trace}")
    if [[ -n "${failed}" ]]; then
        echo " Failed tasks (see .command.err in each work folder):"
        echo "${failed}"
    fi
fi
echo "================================================="

exit "${status}"

# After the run: figures and transmission pairs for this species. In bin/plot_results.sbatch set
#   RESULTS="/var/scratch/global/${USER}/neonatal-amr/outdir_ethiopia_kpneumoniae"   (or _ecoli)
#   METADATA="<ASSEMBLY_DIR>/metadata_assemblies_kpneumoniae.tsv"                    (or _ecoli)
#   AST=""                                                                           (no AST table)
# then:  bash bin/plot_results.sbatch
