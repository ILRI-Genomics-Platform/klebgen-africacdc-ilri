#!/usr/bin/env bash
#
# Step 03 - Oxford Nanopore read filtering and QC          (pipeline: READ_QC -> FILTLONG, NANOPLOT)
#
#   bash step_by_step/03_long_read_qc.sh [sample ... | all]
#
# Only for samples with a long_fastq in the samplesheet (others are skipped).
# What it does: Filtlong removes reads shorter than 1 kb and keeps the best 90 % of bases
# (--min_length 1000 --keep_percent 90); NanoPlot summarises the filtered reads.
# Output: 01_qc/filtlong/<sample>.filtlong.fastq.gz (used by steps 05), .log;
#         01_qc/nanoplot/<sample>/ (<sample>_NanoPlot-report.html, <sample>_NanoStats.txt)
# Look at: read N50 and total yield (aim for > 50x coverage of a ~5.5 Mb genome, i.e. > 275 Mb), mean quality.

source "$(dirname "$0")/config.sh"

fout="${OUTDIR}/01_qc/filtlong"; mkdir -p "$fout"
for s in $(samples_from_args "$@"); do
    long=$(path "$s" long_fastq)
    if [ -z "$long" ]; then info "$s: no long reads; skipped"; continue; fi
    step "Filtlong + NanoPlot: $s"
    need_file "$long" "check the samplesheet"
    cd "$fout"
    run_sh FILTLONG "filtlong \
        --min_length 1000 --keep_percent 90 \
        $long \
        2> ${s}.filtlong.log \
        | gzip -n > ${s}.filtlong.fastq.gz"

    nout="${OUTDIR}/01_qc/nanoplot/${s}"; mkdir -p "$nout"; cd "$nout"
    run NANOPLOT NanoPlot \
        --prefix "${s}_" \
        -t "${THREADS}" \
        --fastq "${fout}/${s}.filtlong.fastq.gz" \
        -o .
    info "report: ${nout}/${s}_NanoPlot-report.html"
done
