#!/usr/bin/env bash
#
# Step 04 - short-read assembly with Shovill (SPAdes)       (pipeline: ASSEMBLY -> SHOVILL)
#
#   bash step_by_step/04_shovill.sh [sample ... | all]
#
# Needs: step 01. Every isolate gets a short-read assembly, also those with ONT reads (step 05).
# What it does: Shovill estimates the genome size, subsamples to ~100x, corrects reads, assembles with
# SPAdes, corrects the contigs with Pilon and removes contigs shorter than 200 bp (--minlen 200).
# Output (02_assembly/shovill/): <sample>.contigs.fa, <sample>.contigs.gfa (open in Bandage), <sample>.shovill.log
# Look at: number of contigs (K. pneumoniae typically 50-150), total length (~5.2-5.9 Mb), the log's warnings.

source "$(dirname "$0")/config.sh"

out="${OUTDIR}/02_assembly/shovill"; mkdir -p "$out"
for s in $(samples_from_args "$@"); do
    step "Shovill: $s"
    need_file "$(trimmed_r1 "$s")" "run step 01 (fastp) first"
    work="${out}/${s}_work"; mkdir -p "$work"; cd "$work"
    run SHOVILL shovill \
        --R1 "$(trimmed_r1 "$s")" \
        --R2 "$(trimmed_r2 "$s")" \
        --minlen 200 \
        --cpus "${THREADS}" \
        --ram "${MEM_GB}" \
        --outdir ./ \
        --force
    # the pipeline publishes these three files with the sample name in front
    for f in contigs.fa contigs.gfa shovill.log; do [ -e "$f" ] && cp "$f" "${out}/${s}.${f}"; done
    info "assembly: ${out}/${s}.contigs.fa ($(grep -c '>' "${out}/${s}.contigs.fa") contigs)"
done
