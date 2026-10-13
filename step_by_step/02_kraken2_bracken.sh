#!/usr/bin/env bash
#
# Step 02 - species composition / contamination check     (pipeline: READ_QC -> KRAKEN2_KRAKEN2, BRACKEN_BRACKEN)
#
#   bash step_by_step/02_kraken2_bracken.sh [sample ... | all]
#
# Needs: step 01 (trimmed reads). Kraken2 loads the whole database into memory: the course's Standard-8
# database (KRAKEN2_DB) needs ~8 GB, so a session with --mem=16G is enough (the full Standard needs ~100 GB).
# What it does: Kraken2 assigns each read pair to a taxon; Bracken re-estimates species abundance
# from Kraken2's report (-l S: species level; -t 10: ignore species with fewer than 10 reads).
# Output (01_qc/kraken2/, 01_qc/bracken/): <sample>.kraken2.report.txt, <sample>.tsv (Bracken species table)
# Look at: is the expected species > 90 % of the reads? A second species with many reads means a
#   mixed culture or contamination; a different top species means the isolate may be mislabelled.

source "$(dirname "$0")/config.sh"
[ -d "${KRAKEN2_DB}" ] || fail "Kraken2 database not found: ${KRAKEN2_DB} (set KRAKEN2_DB)"

kout="${OUTDIR}/01_qc/kraken2"; bout="${OUTDIR}/01_qc/bracken"; mkdir -p "$kout" "$bout"
for s in $(samples_from_args "$@"); do
    step "Kraken2 + Bracken: $s"
    need_file "$(trimmed_r1 "$s")" "run step 01 (fastp) first"
    cd "$kout"
    run KRAKEN2_KRAKEN2 kraken2 \
        --db "${KRAKEN2_DB}" \
        --threads "${THREADS}" \
        --report "${s}.kraken2.report.txt" \
        --gzip-compressed \
        --output /dev/null \
        --paired \
        "$(trimmed_r1 "$s")" "$(trimmed_r2 "$s")"

    cd "$bout"
    run BRACKEN_BRACKEN bracken \
        -l S -t 10 \
        -d "${KRAKEN2_DB}" \
        -i "${kout}/${s}.kraken2.report.txt" \
        -o "${s}.tsv" \
        -w "${s}.kraken2.report_bracken.txt"

    info "top species (Bracken):"
    awk 'NR > 1' "${s}.tsv" | sort -t$'\t' -k7,7gr | awk -F'\t' 'NR <= 4 { printf "       %-40s %6.2f %%\n", $1, $7 * 100 }'
done
