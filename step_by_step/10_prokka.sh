#!/usr/bin/env bash
#
# Step 10 (optional) - genome annotation with Prokka           (pipeline: PROKKA, with --run_prokka)
#
#   bash step_by_step/10_prokka.sh [sample ... | all]
#
# Needs: an assembly. Not needed for typing, AMR or phylogeny; useful to look at genes around an AMR gene
# or to submit a genome. Takes ~5-10 minutes per isolate.
# Output (02_assembly/prokka/<sample>/): <sample>.gff, .gbk, .faa (proteins), .ffn (genes), .tsv, .txt (summary)

source "$(dirname "$0")/config.sh"

for s in $(samples_from_args "$@"); do
    asm=$(best_assembly "$s"); need_file "$asm" "run step 04 (Shovill) first"
    case "$(field "$s" species)" in
        ecoli) genus=Escherichia; species=coli ;;
        *)     genus=Klebsiella;  species=pneumoniae ;;
    esac
    step "Prokka ($genus $species): $s"
    out="${OUTDIR}/02_assembly/prokka"; mkdir -p "$out"; cd "$out"
    run PROKKA prokka \
        --kingdom Bacteria --genus "$genus" --species "$species" --strain "$s" \
        --mincontiglen 200 --addgenes \
        --cpus "${THREADS}" \
        --prefix "$s" \
        --outdir "$s" --force \
        "$asm"
    cat "${s}/${s}.txt"
done
