#!/usr/bin/env bash
#
# Step 08 - AMR genes and point mutations with AMRFinderPlus  (pipeline: TYPING_AMR -> AMRFINDERPLUS_UPDATE, AMRFINDERPLUS_RUN)
#
#   bash step_by_step/08_amrfinderplus.sh [sample ... | all]
#
# Needs: an assembly, and an AMRFinderPlus database: set AMRFINDER_DB to a .tar.gz or folder
# (e.g. the pipeline's <outdir>/databases/amrfinderdb.tar.gz), or leave it empty to download the latest
# one into ${OUTDIR}/databases (needs internet: run that once on the login node).
# What it does: searches the assembly for acquired resistance genes, and, with --organism, for
# species-specific resistance mutations (e.g. gyrA/parC for fluoroquinolones, ompK36 porin changes).
# --plus adds stress, heat, biocide and virulence genes.
# Output (04_amr/amrfinderplus/): <sample>.tsv (one row per gene/mutation found),
#   <sample>-mutations.tsv (every known mutation position examined)
# Look at: Class/Subclass (e.g. BETA-LACTAM / CARBAPENEM), Element subtype (AMR vs POINT), % identity
#   and % coverage; compare with Kleborate's AMR columns.

source "$(dirname "$0")/config.sh"

dbdir="${OUTDIR}/databases/amrfinderdb"
if [ -z "${AMRFINDER_DB}" ]; then
    if [ ! -d "${dbdir}/latest" ]; then
        step "AMRFinderPlus database: downloading the latest (needs internet)"
        mkdir -p "${OUTDIR}/databases"
        run AMRFINDERPLUS_UPDATE amrfinder_update -d "${dbdir}"
    fi
    db="${dbdir}/latest"
elif [[ "${AMRFINDER_DB}" == *.tar.gz ]]; then
    db="${OUTDIR}/databases/amrfinderdb_unpacked"
    if [ ! -d "$db" ]; then mkdir -p "$db"; tar xzf "${AMRFINDER_DB}" -C "$db"; fi
else
    db="${AMRFINDER_DB}"
fi
info "database: $db"

out="${OUTDIR}/04_amr/amrfinderplus"; mkdir -p "$out"; cd "$out"
for s in $(samples_from_args "$@"); do
    asm=$(best_assembly "$s"); need_file "$asm" "run step 04 (Shovill) first"
    case "$(field "$s" species)" in
        ecoli) organism=Escherichia ;;
        *)     organism=Klebsiella_pneumoniae ;;
    esac
    step "AMRFinderPlus (--organism $organism): $s"
    run_sh AMRFINDERPLUS_RUN "amrfinder \
        -n $asm \
        --organism $organism --mutation_all ${s}-mutations.tsv \
        --plus \
        --database $db \
        --threads ${THREADS} > ${s}.tsv"
    info "$(($(wc -l < "${s}.tsv") - 1)) genes/mutations found; AMR genes and mutations:"
    awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) { if ($i ~ /^(Element symbol|Gene symbol)$/) g = i; if ($i ~ /^(Type|Element type)$/) t = i; if ($i == "Subclass") c = i }; next }
                $t == "AMR" { printf "       %-20s %s\n", $g, $c }' "${s}.tsv"
done
