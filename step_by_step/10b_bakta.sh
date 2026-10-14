#!/usr/bin/env bash
#
# Step 10b (optional) - genome annotation with Bakta           (pipeline: BAKTA_BAKTA, with --run_bakta)
#
#   bash step_by_step/10b_bakta.sh [sample ... | all]
#
# Needs: an assembly and the Bakta database (BAKTA_DB in config.sh: the trainer's light database, v6).
# An alternative to Prokka (step 10): Bakta uses a much larger database of known proteins, so more genes get a
# name and a function, and it marks known AMR genes (via AMRFinderPlus) and plasmid replicons. Not needed for
# typing, AMR or phylogeny. Takes ~10-30 minutes per isolate with 1 CPU.
# Output (02_assembly/bakta/<sample>/): <sample>.gff3, .gbff, .embl, .faa (proteins), .ffn (genes), .tsv
#   (one row per feature), .txt (summary), .json, and plots (.png/.svg)
# Look at: the .txt summary (compare the number of CDSs and hypothetical proteins with Prokka's), and the
#   .tsv around an AMR gene from step 08.

source "$(dirname "$0")/config.sh"

[ -s "${BAKTA_DB}/version.json" ] \
    || fail "no Bakta database in ${BAKTA_DB} (set BAKTA_DB; the trainer guide shows how to download it)"
info "Bakta database: ${BAKTA_DB} ($(tr -d '\n ' < "${BAKTA_DB}/version.json"))"

for s in $(samples_from_args "$@"); do
    asm=$(best_assembly "$s"); need_file "$asm" "run step 04 (Shovill) first"
    case "$(field "$s" species)" in
        ecoli) genus=Escherichia; species=coli ;;
        *)     genus=Klebsiella;  species=pneumoniae ;;
    esac
    step "Bakta ($genus $species): $s"
    out="${OUTDIR}/02_assembly/bakta/${s}"; mkdir -p "$out"; cd "$out"
    mkdir -p bakta_home     # as the pipeline: a writable HOME for fontconfig
    run_sh BAKTA_BAKTA "export HOME=\$PWD/bakta_home; bakta \
        $asm \
        --genus $genus --species $species --strain $s --min-contig-length 200 \
        --threads ${THREADS} \
        --prefix $s \
        --output . --force \
        --db ${BAKTA_DB}"
    cat "${s}.txt"
done
