#!/usr/bin/env bash
#
# Step 00 - check the set-up before running anything
#
#   bash step_by_step/00_check_setup.sh
#
# Shows the settings from config.sh, the samples, and the version of every tool (which also proves that
# each container image, or module, works).

source "$(dirname "$0")/config.sh"

step "Settings"
info "samplesheet : ${SAMPLESHEET} ($(all_samples | wc -l) samples: $(all_samples | tr '\n' ' '))"
info "input type  : $(is_assembly_sheet && echo assemblies || echo reads)"
info "results     : ${OUTDIR}"
info "tools from  : ${TOOL_MODE}$([ "${TOOL_MODE}" = singularity ] && echo " images in ${CONTAINERS}")"
info "threads/RAM : ${THREADS} / ${MEM_GB} GB"
info "Kraken2 DB  : ${KRAKEN2_DB} $([ -d "${KRAKEN2_DB}" ] && echo '(found)' || echo '(NOT FOUND)')"
info "references  : ${REFERENCES} $([ -d "${REFERENCES}" ] && echo "($(ls "${REFERENCES}"/*.fasta 2>/dev/null | wc -l) FASTA files)" || echo '(NOT FOUND)')"
info "AMRFinder DB: ${AMRFINDER_DB:-(will be downloaded in step 08)}"

step "Input files"
missing=0
for s in $(all_samples); do
    for c in fastq_1 fastq_2 long_fastq fasta; do
        f=$(path "$s" "$c"); [ -z "$f" ] && continue
        if [ -r "$f" ]; then info "ok       $s $c"; else warn "missing  $s $c: $f"; missing=$((missing + 1)); fi
    done
done

step "Tool versions"
check() {   # check PROCESS 'version command(s)'
    local out status=0
    out=$(in_tool_env "$1" bash -c "$2" 2>&1) || status=$?
    out=$(printf '%s\n' "$out" | grep -v '^\s*$' | head -n 3 | tr '\n' ' ')
    if [ "$status" -eq 0 ] && [ -n "$out" ]; then
        printf '     %-22s %s\n' "$1" "$out"
    else
        printf '\033[1;31m     %-22s FAILED: %s\033[0m\n' "$1" "$out"; missing=$((missing + 1))
    fi
}
if ! is_assembly_sheet; then                       # read steps 01-05
    check FASTP                "fastp --version 2>&1"
    check KRAKEN2_KRAKEN2      "kraken2 --version"
    check BRACKEN_BRACKEN      "bracken -v 2>&1"
    check SHOVILL              "shovill --version"
    if [ -n "$(awk -F',' 'NR == 1 { for (i = 1; i <= NF; i++) if ($i == "long_fastq") c = i; next } c && $c != "" { print; exit }' "${SAMPLESHEET}")" ]; then
        check FILTLONG             "filtlong --version"
        check NANOPLOT             "NanoPlot --version"
        check AUTOCYCLER_ASSEMBLE  "autocycler --version; flye --version; raven --version"
        check MEDAKA               "medaka --version"
        check POLYPOLISH_SHORTREAD "polypolish --version"
        check PYPOLCA_RUN          "pypolca --version"
    fi
fi
check QUAST                "quast.py --version 2>&1"
check KLEBORATE            "kleborate --version"
check KAPTIVE              "kaptive --version"
has_ecoli=no
for s in $(all_samples); do [ "$(field "$s" species)" = ecoli ] && has_ecoli=yes; done
if [ "$has_ecoli" = yes ]; then
    check ECTYPER              "ectyper --version"
    check EZCLERMONT           "ezclermont --version 2>&1 || pip show ezclermont 2>/dev/null | grep Version"
fi
check AMRFINDERPLUS_RUN    "amrfinder --version"
check MOBSUITE_RECON       "mob_recon --version"
check SNIPPY_RUN           "snippy --version 2>&1"
check GUBBINS              "run_gubbins.py --version 2>&1"
check IQTREE               "iqtree --version"
check SNPDISTS             "snp-dists -v"
check MULTIQC              "multiqc --version"
required=$missing
check PROKKA               "prokka --version 2>&1"           # optional (step 10)
[ "$missing" -gt "$required" ] && warn "Prokka is not available: step 10 (optional) will not run"
missing=$required

echo
if [ "$missing" -eq 0 ]; then info "All set."; else warn "$missing problem(s) above"; exit 1; fi
