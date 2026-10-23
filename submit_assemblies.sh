#!/bin/bash
#
# submit_assemblies.sh - one country, end to end: download the Terra assemblies, filter to neonatal blood
# isolates, build the samplesheet, and run neonatal-amr-nf on them (ILRI HPC).
#
# Everything is derived from COUNTRY, so the same script serves every country. Nothing is hard-coded.
#
#   STEPS picks what runs (default: all four, in this order)
#     download     bin/download_terra_outputs.py  -> $COUNTRY_DIR/assemblies/<sample>_contigs.fasta
#     filter       bin/filter_terra_table.py      -> $FILTERED   (blood, <= MAX_AGE_DAYS days, this species)
#     samplesheet  bin/make_assembly_samplesheet.py -> $INPUT
#     run          nextflow run main.nf
#
# IMPORTANT: 'download' needs internet and a logged-in gcloud ('gcloud auth login'), which compute nodes
# usually do not have. Run the preparation on the LOGIN node, then submit the pipeline:
#
#     mkdir -p logs
#     COUNTRY=Ethiopia STEPS="download filter samplesheet" bash submit_assemblies.sh      # login node
#     sbatch --export=ALL,COUNTRY=Ethiopia,STEPS=run submit_assemblies.sh                 # compute node
#
# Or, where compute nodes do have internet and gcloud credentials, all four in one job:
#
#     sbatch --export=ALL,COUNTRY=Ethiopia submit_assemblies.sh
#
# More examples:
#     sbatch --export=ALL,COUNTRY=Kenya,SPECIES=ecoli submit_assemblies.sh
#     sbatch --export=ALL,COUNTRY=Malawi,MAX_AGE_DAYS=30 submit_assemblies.sh
#     sbatch --export=ALL,COUNTRY=Zimbabwe,TERRA_TSV=/path/to/Zimbabwe_sample.tsv submit_assemblies.sh
#     COUNTRY=Ethiopia STEPS=samplesheet bash submit_assemblies.sh        # rebuild one samplesheet only
#
# Folder layout under DATA_DIR (set DATA_DIR=$SHARED/countries to share it with the participants):
#     $DATA_DIR/terra/<COUNTRY>.tsv                       the Terra table you exported (input)
#     $DATA_DIR/<COUNTRY>/assemblies/<sample>_contigs.fasta
#     $DATA_DIR/<COUNTRY>/terra_downloads_<COUNTRY>.tsv    download log
#     $DATA_DIR/<COUNTRY>/<COUNTRY>_le<N>d.tsv             filtered metadata
#     $DATA_DIR/<COUNTRY>/<COUNTRY>_le<N>d_<species>_samplesheet.csv
#     $DATA_DIR/<COUNTRY>/<COUNTRY>_le<N>d_<species>_metadata.tsv   for bin/plot_results.sbatch
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

fail() { echo "[ERROR] $*" >&2; exit 1; }
step() { printf '\n==== %s  (%s)\n' "$*" "$(date '+%F %T')"; }
want() { [[ " ${STEPS} " == *" $1 "* ]]; }

echo "================================================="
echo "Node Hostname : $(hostname)"
echo "Job ID        : ${SLURM_JOB_ID:-none}"
echo "Start Time    : $(date)"
echo "================================================="

#-------------------------------------------------------------------------------
# Settings (each can be overridden with sbatch --export=ALL,NAME=value, or NAME=value in the environment)
#-------------------------------------------------------------------------------
COUNTRY="${COUNTRY:-Ethiopia}"                 # the one setting that has to be right
SPECIES="${SPECIES:-kpneumoniae}"              # kpneumoniae | ecoli
MAX_AGE_DAYS="${MAX_AGE_DAYS:-28}"             # neonates: age <= this many days
STEPS="${STEPS:-download filter samplesheet run}"

PIPELINE_DIR="${PIPELINE_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"       # folder with main.nf (where you ran sbatch)
USER="${USER:-$(id -un)}"                      # set -u: not every shell exports USER
BASE_DIR="${BASE_DIR:-/var/scratch/global/${USER}/neonatal-amr}"
SHARED="${SHARED:-/var/scratch/global/jjuma/neonatal-amr/shared}"   # course shared folder (databases)

# Where the country data lives. $SHARED/countries makes it readable by the participants.
DATA_DIR="${DATA_DIR:-${HOME}/trainings/klebgen-africacdc-ilri/Terra_outputs}"
COUNTRY_DIR="${COUNTRY_DIR:-${DATA_DIR}/${COUNTRY}}"
TERRA_TSV="${TERRA_TSV:-${DATA_DIR}/terra/${COUNTRY}.tsv}"          # the Terra export for this country
ASSEMBLY_DIR="${ASSEMBLY_DIR:-${COUNTRY_DIR}/assemblies}"

TAG="${TAG:-${COUNTRY}_le${MAX_AGE_DAYS}d_${SPECIES}}"
FILTERED="${FILTERED:-${COUNTRY_DIR}/${COUNTRY}_le${MAX_AGE_DAYS}d.tsv}"
INPUT="${INPUT:-${COUNTRY_DIR}/${TAG}_samplesheet.csv}"
METADATA_OUT="${METADATA_OUT:-${COUNTRY_DIR}/${TAG}_metadata.tsv}"

# Terra column names (these differ between workspaces: check the header of your own table)
COUNTRY_COL="${COUNTRY_COL:-Country}"
ORGANISM_COL="${ORGANISM_COL:-gambit_predicted_taxon}"
SPECIMEN_COL="${SPECIMEN_COL:-Specimen_type}"
AGE_COL="${AGE_COL:-Age}"
GROUP_COL="${GROUP_COL:-kleborate_klebsiella_mlst_sequence_type}"
COUNTRY_MAP="${COUNTRY_MAP:-}"                 # optional TSV 'raw<TAB>country' for GAZ-only country values

case "${SPECIES}" in
    kpneumoniae) ORGANISM="${ORGANISM:-Klebsiella pneumoniae}" ;;
    ecoli)       ORGANISM="${ORGANISM:-Escherichia coli}" ;;
    *) fail "SPECIES must be kpneumoniae or ecoli (got '${SPECIES}')" ;;
esac

RUN_NAME="${RUN_NAME:-${TAG}}"
OUTDIR="${OUTDIR:-${BASE_DIR}/outdir_${RUN_NAME}}"
WORKDIR="${WORKDIR:-${BASE_DIR}/work}"
MIN_GROUP="${MIN_GROUP:-3}"                                       # --phylo_min_group_size
PHYLO_METHOD="${PHYLO_METHOD:-ska2}"                              # --phylo_method: snippy, ska2 or both
# --phylo_grouping: species   = all isolates of the species in one group, whatever the ST (one tree, one SNP matrix)
#                   kleborate = one group per species + ST (one tree and SNP matrix per ST)
#                   samplesheet = the samplesheet's phylo_group column
# Isolates without an ST from Kleborate are left out in species and kleborate modes.
PHYLO_GROUPING="${PHYLO_GROUPING:-species}"
AMRFINDER_DB="${AMRFINDER_DB:-${SHARED}/databases/amrfinderplus/amrfinderdb.tar.gz}"   # course/TRAINER_GUIDE.md, 5.1

# Nextflow keeps its run history (needed by -resume) in the folder it starts from.
# A fixed launch folder per run name means a resubmission always resumes the same run.
LAUNCH_DIR="${LAUNCH_DIR:-${BASE_DIR}/launch/${RUN_NAME}}"

PYTHON="${PYTHON:-python3}"
BIN="${PIPELINE_DIR}/bin"

echo " Country    : ${COUNTRY} (${SPECIES}, age <= ${MAX_AGE_DAYS} d)"
echo " Steps      : ${STEPS}"
echo " Pipeline   : ${PIPELINE_DIR}"
echo " Terra table: ${TERRA_TSV}"
echo " Country dir: ${COUNTRY_DIR}"
echo "================================================="

for s in ${STEPS}; do
    case "$s" in download|filter|samplesheet|run) ;; *) fail "unknown step '$s' in STEPS (download filter samplesheet run)" ;; esac
done
[[ -f "${PIPELINE_DIR}/main.nf" ]] || fail "main.nf not found in ${PIPELINE_DIR}; run from the pipeline folder or set PIPELINE_DIR"
mkdir -p "${COUNTRY_DIR}"

#-------------------------------------------------------------------------------
# 1. Download the assemblies from Terra
#-------------------------------------------------------------------------------
if want download; then
    step "Download assemblies for ${COUNTRY}"
    [[ -s "${TERRA_TSV}" ]] || fail "Terra table not found: ${TERRA_TSV} (set TERRA_TSV)"
    if [[ -n "${SLURM_JOB_ID:-}" ]]; then
        echo "[NOTE] downloading inside a Slurm job: compute nodes often have no internet and no gcloud" >&2
        echo "       credentials. If this step fails, run it on the login node:" >&2
        echo "       COUNTRY=${COUNTRY} STEPS=download bash submit_assemblies.sh" >&2
    fi
    # --by-country writes to <outdir>/<Country>/..., so the parent of the country folder is the outdir
    "${PYTHON}" "${BIN}/download_terra_outputs.py" \
        --terra-output "${TERRA_TSV}" \
        --outdir "$(dirname "${COUNTRY_DIR}")" \
        --files assembly_fasta \
        --by-country \
        --country-col "${COUNTRY_COL}" \
        --countries "${COUNTRY}" \
        ${COUNTRY_MAP:+--country-map "${COUNTRY_MAP}"}
fi

#-------------------------------------------------------------------------------
# 2. Keep the neonatal blood isolates of this species
#-------------------------------------------------------------------------------
if want filter; then
    step "Filter ${COUNTRY}: blood, age <= ${MAX_AGE_DAYS} d, ${ORGANISM}"
    [[ -s "${TERRA_TSV}" ]] || fail "Terra table not found: ${TERRA_TSV} (set TERRA_TSV)"
    "${PYTHON}" "${BIN}/filter_terra_table.py" \
        --terra-output "${TERRA_TSV}" \
        --out "${FILTERED}" \
        --specimen-col "${SPECIMEN_COL}" \
        --specimen-contains \
        --organism-col "${ORGANISM_COL}" \
        --organism "${ORGANISM}" \
        --country-col "${COUNTRY_COL}" \
        --countries "${COUNTRY}" \
        --age-col "${AGE_COL}" \
        --max-age-days "${MAX_AGE_DAYS}"
fi

#-------------------------------------------------------------------------------
# 3. Build the samplesheet from the filtered table and the downloaded contigs
#-------------------------------------------------------------------------------
if want samplesheet; then
    step "Samplesheet for ${COUNTRY} (${SPECIES})"
    [[ -s "${FILTERED}" ]] || fail "filtered table not found: ${FILTERED} (run the 'filter' step first)"
    [[ -d "${ASSEMBLY_DIR}" ]] || fail "contigs folder not found: ${ASSEMBLY_DIR} (run the 'download' step first)"
    # The Terra table's first column is 'entity:<table>_sample_id' and differs per workspace, so read it
    ID_COL="${ID_COL:-$(head -n 1 "${FILTERED}" | cut -f1)}"
    echo "       sample ID column: ${ID_COL}"
    "${PYTHON}" "${BIN}/make_assembly_samplesheet.py" \
        --metadata "${FILTERED}" \
        --contigs-dir "${ASSEMBLY_DIR}" \
        --id-col "${ID_COL}" \
        --species-col "${ORGANISM_COL}" \
        --species "${SPECIES}" \
        --group-col "${GROUP_COL}" \
        --out "${INPUT}" \
        --metadata-out "${METADATA_OUT}"
fi

# if ! want run; then
#     echo; echo "================================================="
#     echo " Preparation done. Submit the pipeline with:"
#     echo "   sbatch --export=ALL,COUNTRY=${COUNTRY},SPECIES=${SPECIES},STEPS=run submit_assemblies.sh"
#     echo "================================================="
#     exit 0
# fi

# #-------------------------------------------------------------------------------
# # Checks before anything is submitted (before loading modules, so a mistake is reported at once)
# #-------------------------------------------------------------------------------
# [[ -f "${PIPELINE_DIR}/conf/slurm.config" ]] || fail "conf/slurm.config not found in ${PIPELINE_DIR}"
# [[ -f "${INPUT}" ]]                          || fail "samplesheet not found: ${INPUT} (run the 'samplesheet' step first)"
# [[ "${PHYLO_GROUPING}" =~ ^(species|kleborate|samplesheet)$ ]] \
#     || fail "PHYLO_GROUPING must be species, kleborate or samplesheet (got '${PHYLO_GROUPING}')"
# grep -q "'species'" "${PIPELINE_DIR}/subworkflows/local/utils_neonatal_amr_pipeline/main.nf" \
#     || [[ "${PHYLO_GROUPING}" != species ]] \
#     || fail "this copy of the pipeline predates --phylo_grouping species: update ${PIPELINE_DIR}"
# [[ -e "${AMRFINDER_DB}" ]] || fail "AMRFinderPlus database not found: ${AMRFINDER_DB} (course/TRAINER_GUIDE.md, section 5.1)"
# head -n 1 "${INPUT}" | grep -q 'fasta' \
#     || fail "${INPUT} has no fasta column: it is not an assembly samplesheet (for reads use submit_neonatal_amr_reads.sbatch)"

# # Every contigs file in the samplesheet must exist and be readable from the compute nodes
# INPUT_DIR=$(cd "$(dirname "${INPUT}")" && pwd)
# missing=0
# while IFS= read -r f; do
#     [[ -z "$f" ]] && continue
#     [[ "$f" = /* ]] || f="${INPUT_DIR}/${f}"
#     [[ -r "$f" ]] || { echo "[ERROR] contigs file missing or unreadable: $f" >&2; missing=$((missing + 1)); }
# done < <(awk -F',' 'NR==1 { for (i = 1; i <= NF; i++) col[$i] = i; next } { print $col["fasta"] }' "${INPUT}")
# (( missing == 0 )) || fail "${missing} contigs file(s) missing"

# INPUT="${INPUT_DIR}/$(basename "${INPUT}")"

# #-------------------------------------------------------------------------------
# # Environment for the pipeline
# #-------------------------------------------------------------------------------
# set +u                                        # module scripts can use unset variables
# module load nextflow/25.04
# set -u

# export NXF_OPTS="-Xms512M -Xmx4G"
# export NXF_ANSI_LOG=false                     # plain progress lines in the log file
# # Singularity's own temporary files on scratch, not in your home quota.
# # (The pipeline's images come from the shared read-only folder; nothing is downloaded.)
# export SINGULARITY_CACHEDIR="${BASE_DIR}/.singularity_cache"
# export APPTAINER_CACHEDIR="${SINGULARITY_CACHEDIR}"
# export SINGULARITY_TMPDIR="${BASE_DIR}/.singularity_tmp"
# export APPTAINER_TMPDIR="${SINGULARITY_TMPDIR}"

# mkdir -p "${LAUNCH_DIR}" "${OUTDIR}" "${WORKDIR}" "${SINGULARITY_CACHEDIR}" "${SINGULARITY_TMPDIR}"

# step "Pipeline"
# echo " Input      : assemblies (${SPECIES})"
# echo " Samplesheet: ${INPUT} ($(($(wc -l < "${INPUT}") - 1)) samples)"
# echo " Outdir     : ${OUTDIR}"
# echo " Work dir   : ${WORKDIR}"
# echo " Launch dir : ${LAUNCH_DIR}"
# echo " AMRFinder  : ${AMRFINDER_DB}"
# echo " Phylogeny  : ${PHYLO_METHOD}, groups: ${PHYLO_GROUPING} (min ${MIN_GROUP} isolates)"
# echo "================================================="

# #-------------------------------------------------------------------------------
# # Run
# #-------------------------------------------------------------------------------
# cd "${LAUNCH_DIR}"

# status=0
# nextflow run "${PIPELINE_DIR}/main.nf" \
#     -profile singularity \
#     -c "${PIPELINE_DIR}/conf/slurm.config" \
#     --input "${INPUT}" \
#     --input_type assembly \
#     --outdir "${OUTDIR}" \
#     -work-dir "${WORKDIR}" \
#     --amrfinderplus_db "${AMRFINDER_DB}" \
#     --phylo_min_group_size "${MIN_GROUP}" \
#     --phylo_grouping "${PHYLO_GROUPING}" \
#     --phylo_method "${PHYLO_METHOD}" \
#     -resume \
#     || status=$?

# #-------------------------------------------------------------------------------
# # Summary
# #-------------------------------------------------------------------------------
# echo "================================================="
# echo " Nextflow exit status: ${status}"
# echo " End Time            : $(date)"
# echo " Nextflow log        : ${LAUNCH_DIR}/.nextflow.log"

# trace=$(ls -t "${OUTDIR}"/pipeline_info/execution_trace_*.txt 2>/dev/null | head -n 1 || true)
# if [[ -n "${trace}" ]]; then
#     echo " Trace               : ${trace}"
#     echo "-------------------------------------------------"
#     echo " Tasks by status:"
#     awk -F'\t' 'NR==1 { for (i = 1; i <= NF; i++) if ($i == "status") s = i; next } { n[$s]++ }
#                 END { for (k in n) printf "   %-10s %d\n", k, n[k] }' "${trace}"
#     failed=$(awk -F'\t' -v wd="${WORKDIR}" 'NR==1 { for (i = 1; i <= NF; i++) { if ($i == "status") s = i; if ($i == "name") m = i; if ($i == "hash") h = i }; next }
#                          $s == "FAILED" { print "   " $m "  ->  " wd "/" $h "*" }' "${trace}")
#     if [[ -n "${failed}" ]]; then
#         echo " Failed tasks (see .command.err in each work folder):"
#         echo "${failed}"
#     fi
# fi
# echo "================================================="

# exit "${status}"

# # After the run: figures and transmission pairs for this country and species. In bin/plot_results.sbatch set
# #   RESULTS="/var/scratch/global/${USER}/neonatal-amr/outdir_<TAG>"      (TAG = <Country>_le<N>d_<species>)
# #   METADATA="<COUNTRY_DIR>/<TAG>_metadata.tsv"
# #   AST=""                                                               (no AST table)
# # then:  bash bin/plot_results.sbatch
