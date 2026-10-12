#!/usr/bin/env bash
#
# Step 09 - plasmid reconstruction and typing with MOB-suite  (pipeline: TYPING_AMR -> MOBSUITE_RECON)
#
#   bash step_by_step/09_mobsuite.sh [sample ... | all]
#
# Needs: an assembly. Works best on hybrid assemblies, where plasmids are complete circular contigs.
# What it does: assigns each contig to the chromosome or to a plasmid, groups plasmid contigs into
# plasmids, and types each plasmid (replicon, relaxase, mobility: conjugative / mobilizable / non-mobilizable).
# Output (05_plasmids/mobsuite/<sample>/results/): contig_report.txt, mobtyper_results.txt, chromosome.fasta, plasmid_*.fasta
# Look at: mobtyper_results.txt (one row per plasmid: size, rep_type, predicted_mobility); then find which
#   plasmid carries the carbapenemase/ESBL genes from step 08 (match contig IDs in contig_report.txt).

source "$(dirname "$0")/config.sh"

for s in $(samples_from_args "$@"); do
    asm=$(best_assembly "$s"); need_file "$asm" "run step 04 (Shovill) first"
    step "MOB-recon: $s"
    out="${OUTDIR}/05_plasmids/mobsuite/${s}"; mkdir -p "$out"; cd "$out"
    run MOBSUITE_RECON mob_recon \
        --infile "$asm" \
        --num_threads "${THREADS}" \
        --outdir results \
        --sample_id "$s" \
        --force
    if [ -s results/mobtyper_results.txt ]; then
        awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) col[$i] = i }
                    { printf "       %-28s %-10s %-35s %s\n", $col["sample_id"], $col["size"], $col["rep_type(s)"], $col["predicted_mobility"] }' \
            results/mobtyper_results.txt
    else
        info "no plasmids found"
    fi
done
