#!/usr/bin/env bash
#
# Step 06 - assembly statistics with QUAST                  (pipeline: ASSEMBLY -> QUAST)
#
#   bash step_by_step/06_quast.sh [sample ... | all]
#
# Needs: step 04 (and 05 for hybrid samples), or an assembly samplesheet.
# Runs on each sample's best assembly: hybrid if there is one, otherwise Shovill's contigs
# (or the samplesheet's contigs for assembly input). No reference genome is used.
# Output (02_assembly/quast/<sample>/): report.txt, report.html, report.tsv
# Look at: # contigs, total length, N50, GC % (K. pneumoniae ~57 %, E. coli ~50.5 %), largest contig.

source "$(dirname "$0")/config.sh"

out="${OUTDIR}/02_assembly/quast"; mkdir -p "$out"; cd "$out"
for s in $(samples_from_args "$@"); do
    step "QUAST: $s"
    asm=$(best_assembly "$s"); need_file "$asm" "run step 04 (Shovill) first"
    info "assembly: $asm"
    run QUAST quast.py \
        --output-dir "$s" \
        --threads "${THREADS}" \
        "$asm"
    cat "${s}/report.txt"
done
