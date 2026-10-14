#!/bin/bash
#
# submit_neonatal_amr_reads.sbatch - run neonatal-amr-nf on SEQUENCING READS (FASTQ) on the ILRI HPC.
# For genomes that are already assembled (contigs), use submit_neonatal_amr_assemblies.sbatch.
#
# The samplesheet has the columns  sample,fastq_1,fastq_2[,long_fastq,...]  (Illumina, optionally with ONT).
#
# This job runs only the Nextflow "head" process: it uses little CPU and memory but must stay alive for the
# whole run. Nextflow submits every analysis step as its own Slurm job (conf/slurm.config chooses the nodes).
#
# From the pipeline folder:
#     mkdir -p logs                          # Slurm does not create the log folder; the job fails without it
#     sbatch submit_neonatal_amr_reads.sbatch                         # the 29-isolate training set (compact)
#
# Databases: the course's shared ones in $SHARED/databases (Kraken2 Standard-8, AMRFinderPlus), for every run.
# Another one for a run:  --export=ALL,KRAKEN2_DB=/export/data/bio/kraken2/db/standard_20241228   (full Standard)
#
# Another samplesheet without editing this file:
#     sbatch --export=ALL,INPUT=/path/samplesheet.csv,RUN_NAME=myrun submit_neonatal_amr_reads.sbatch
#     (results go to <BASE_DIR>/outdir_<RUN_NAME>)
#
# The trainer's pre-computed runs for the course (course/TRAINER_GUIDE.md, section 8.1), from $SHARED/neonatal-amr-nf:
#     sbatch --export=ALL,PRESET=training submit_neonatal_amr_reads.sbatch   # 29 isolates -> $SHARED/precomputed/pipeline_results
#     sbatch --export=ALL,PRESET=capstone submit_neonatal_amr_reads.sbatch   # capstone    -> $SHARED/capstone/precomputed
#
# Follow it:   squeue -u $USER       tail -f logs/<user>_reads_<jobid>.txt
# Stop it:     scancel <jobid>       (the analysis jobs stop too; submit again and -resume continues)

#SBATCH --job-name=amr-reads
#SBATCH --partition=batch
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=3-00:00:00
#SBATCH --output=logs/%u_reads_%j.txt
#SBATCH --error=logs/%u_reads_error_%j.txt
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

# Pre-computed course runs (PRESET=training or capstone) set the samplesheet, results, work folder and databases
SHARED="${SHARED:-/var/scratch/global/jjuma/neonatal-amr/shared}"
PRESET="${PRESET:-}"
case "${PRESET}" in
    training)
        RUN_NAME="${RUN_NAME:-precomputed_training}"
        INPUT="${INPUT:-${SHARED}/samplesheet.csv}"
        OUTDIR="${OUTDIR:-${SHARED}/precomputed/pipeline_results}"
        WORKDIR="${WORKDIR:-${BASE_DIR}/work_precomputed}" ;;
    capstone)
        RUN_NAME="${RUN_NAME:-precomputed_capstone}"
        INPUT="${INPUT:-${SHARED}/capstone/samplesheet.csv}"
        OUTDIR="${OUTDIR:-${SHARED}/capstone/precomputed}"
        WORKDIR="${WORKDIR:-${BASE_DIR}/work_capstone}" ;;
    "") ;;
    *) echo "[ERROR] PRESET must be 'training' or 'capstone' (got '${PRESET}')" >&2; exit 1 ;;
esac

RUN_NAME="${RUN_NAME:-compact}"
INPUT="${INPUT:-${PIPELINE_DIR}/assets/samplesheet_training_set_compact.csv}"
# 5-sample test instead:
#   sbatch --export=ALL,RUN_NAME=test,INPUT=${PWD}/assets/samplesheet_test_5samples.csv submit_neonatal_amr_reads.sbatch
OUTDIR="${OUTDIR:-${BASE_DIR}/outdir_${RUN_NAME}}"

WORKDIR="${WORKDIR:-${BASE_DIR}/work}"
KRAKEN2_DB="${KRAKEN2_DB:-${SHARED}/databases/kraken2}"            # shared Kraken2 Standard-8 (course/TRAINER_GUIDE.md, 5.1)
AMRFINDER_DB="${AMRFINDER_DB:-${SHARED}/databases/amrfinderplus/amrfinderdb.tar.gz}"   # shared AMRFinderPlus database (course/TRAINER_GUIDE.md, 5.1)
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
[[ -f "${INPUT}" ]]                          || fail "samplesheet not found: ${INPUT}"
head -n 1 "${INPUT}" | grep -q 'fastq_1' \
    || fail "${INPUT} has no fastq_1 column: it is not a reads samplesheet (for contigs use submit_neonatal_amr_assemblies.sbatch)"
[[ -d "${KRAKEN2_DB}" ]]                     || fail "Kraken2 database folder not found: ${KRAKEN2_DB}"
[[ -e "${AMRFINDER_DB}" ]] || fail "AMRFinderPlus database not found: ${AMRFINDER_DB} (course/TRAINER_GUIDE.md, section 5.1)"
KRAKEN2_DB="${KRAKEN2_DB%/}"                 # no trailing slash in the Singularity bind

# Every read file in the samplesheet must exist and be readable from the compute nodes
INPUT_DIR=$(cd "$(dirname "${INPUT}")" && pwd)
missing=0
while IFS= read -r f; do
    [[ -z "$f" ]] && continue
    [[ "$f" = /* ]] || f="${INPUT_DIR}/${f}"
    [[ -r "$f" ]] || { echo "[ERROR] read file missing or unreadable: $f" >&2; missing=$((missing + 1)); }
done < <(awk -F',' 'NR==1 { for (i = 1; i <= NF; i++) col[$i] = i; next }
                     { print $col["fastq_1"]; print $col["fastq_2"]; if ("long_fastq" in col) print $col["long_fastq"] }' "${INPUT}")
(( missing == 0 )) || fail "${missing} read file(s) missing"

INPUT="${INPUT_DIR}/$(basename "${INPUT}")"
mkdir -p "${LAUNCH_DIR}" "${OUTDIR}" "${WORKDIR}" "${SINGULARITY_CACHEDIR}" "${SINGULARITY_TMPDIR}"

echo " Input      : reads"
echo " Pipeline   : ${PIPELINE_DIR}"
echo " Samplesheet: ${INPUT} ($(($(wc -l < "${INPUT}") - 1)) samples)"
echo " Outdir     : ${OUTDIR}"
echo " Work dir   : ${WORKDIR}"
echo " Launch dir : ${LAUNCH_DIR}"
echo " Kraken2 DB : ${KRAKEN2_DB}"
echo " AMRFinder  : ${AMRFINDER_DB}"
[[ -n "${PRESET}" ]] && echo " Preset     : ${PRESET}"
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
    --input_type reads \
    --outdir "${OUTDIR}" \
    -work-dir "${WORKDIR}" \
    --kraken2_db "${KRAKEN2_DB}" \
    --amrfinderplus_db "${AMRFINDER_DB}" \
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

# After the run: figures and transmission pairs. In bin/plot_results.sbatch set
#   RESULTS="/var/scratch/global/${USER}/neonatal-amr/outdir_compact"
#   METADATA="${PIPELINE_DIR}/training_data/metadata_compact.tsv"
#   AST="${PIPELINE_DIR}/training_data/ast_compact.tsv"
# then:  bash bin/plot_results.sbatch
