#!/usr/bin/env bash
#
# Step 01 - Illumina read trimming and QC with fastp            (pipeline: READ_QC -> FASTP)
#
#   bash step_by_step/01_fastp.sh              # first sample in the samplesheet
#   bash step_by_step/01_fastp.sh CHI11E       # one or more samples
#   bash step_by_step/01_fastp.sh all          # every sample
#
# What it does: removes adapters (detected automatically for paired reads), trims low-quality bases
# (--qualified_quality_phred 20) and drops reads shorter than 50 bp (--length_required 50).
# Output (01_qc/fastp/): <sample>_R1/_R2.fastp.fastq.gz (trimmed reads, used by every later step),
#   <sample>.fastp.html (open it in a browser), .json (read by MultiQC), .log
# Look at: reads and bases before/after filtering, Q30 rate, duplication, insert size, adapter content.

source "$(dirname "$0")/config.sh"
is_assembly_sheet && fail "the samplesheet has assemblies (fasta column): there are no reads to trim"

out="${OUTDIR}/01_qc/fastp"; mkdir -p "$out"
for s in $(samples_from_args "$@"); do
    step "fastp: $s"
    r1=$(path "$s" fastq_1); r2=$(path "$s" fastq_2)
    need_file "$r1" "check the samplesheet"; need_file "$r2" "check the samplesheet"
    cd "$out"
    run_sh FASTP "fastp \
        --in1 $r1 \
        --in2 $r2 \
        --out1 ${s}_R1.fastp.fastq.gz \
        --out2 ${s}_R2.fastp.fastq.gz \
        --json ${s}.fastp.json \
        --html ${s}.fastp.html \
        --thread ${THREADS} \
        --detect_adapter_for_pe \
        --qualified_quality_phred 20 --length_required 50 \
        2> ${s}.fastp.log"
    info "report: ${out}/${s}.fastp.html"
done
