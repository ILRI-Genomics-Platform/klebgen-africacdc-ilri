#!/usr/bin/env bash
#
# Step 12 - one summary report with MultiQC                 (pipeline: MULTIQC)
#
#   bash step_by_step/12_multiqc.sh
#
# What it does: scans the results folder for reports it recognises (fastp, Kraken2, NanoPlot, QUAST, ...)
# and puts them in one interactive HTML page.
# Output (multiqc/): multiqc_report.html, multiqc_data/

source "$(dirname "$0")/config.sh"

out="${OUTDIR}/multiqc"; mkdir -p "$out"; cd "$out"
step "MultiQC"
inputs=(); for d in "${OUTDIR}/01_qc" "${OUTDIR}/02_assembly/quast"; do [ -d "$d" ] && inputs+=("$d"); done
[ ${#inputs[@]} -gt 0 ] || fail "no QC or QUAST results yet: run steps 01 and 06 first"
run MULTIQC multiqc --force --title "Neonatal AMR genomic surveillance" --outdir . "${inputs[@]}"
info "report: ${out}/multiqc_report.html"
