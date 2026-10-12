#!/usr/bin/env bash
#
# Step 07 - species, MLST, virulence, AMR and antigen typing
#           (pipeline: TYPING_AMR -> KLEBORATE, KAPTIVE, ECTYPER, EZCLERMONT)
#
#   bash step_by_step/07_typing.sh [sample ... | all]
#
# Needs: an assembly (step 04/05, or an assembly samplesheet).
#   Kleborate  : every isolate. Preset from the samplesheet's species column:
#                kpneumoniae -> kpsc (K. pneumoniae species complex: species, MLST, virulence loci
#                (ybt, clb, iuc, iro, rmp), AMR genes and the K/O loci via Kaptive)
#                ecoli       -> escherichia (E. coli MLST, Achtman and Pasteur schemes)
#   Kaptive    : isolates that Kleborate places in the K. pneumoniae species complex: the full K (capsule)
#                and O (LPS) locus report, with the locus sequences
#   ECTyper, EzClermont : E. coli only: O:H serotype and Clermont phylogroup
# Output: 03_typing/kleborate/<sample>.kleborate.txt, 03_typing/kaptive/<sample>.kaptive_{K,O}.tsv,
#         03_typing/ectyper/<sample>.tsv, 03_typing/ezclermont/<sample>.ezclermont.tsv
# Look at: species and ST, virulence score, resistance score (0-3), K and O locus and their confidence.

source "$(dirname "$0")/config.sh"

# Species reported by Kleborate (column 'species', or '<preset>__species__species' in full headers)
kleborate_value() {   # kleborate_value TABLE species|ST
    awk -F'\t' -v want="$2" 'NR == 1 { for (i = 1; i <= NF; i++) {
                                  if (want == "species" && ($i == "species" || $i ~ /__species__species$/)) c = i
                                  if (want == "ST" && ($i == "ST" || $i ~ /__mlst(_achtman)?__ST$/)) c = i }
                              next }
                             NR == 2 && c { print $c }' "$1"
}

kout="${OUTDIR}/03_typing/kleborate"; mkdir -p "$kout"
for s in $(samples_from_args "$@"); do
    asm=$(best_assembly "$s"); need_file "$asm" "run step 04 (Shovill) first"
    species=$(field "$s" species)
    preset=$([ "$species" = ecoli ] && echo escherichia || echo kpsc)

    step "Kleborate (-p $preset): $s"
    work="${kout}/${s}_work"; rm -rf "$work"; mkdir -p "$work"; cd "$work"
    run KLEBORATE kleborate -a "$asm" -p "$preset" -o kleborate_out
    # Kleborate names its table after the species it found; take the main table (as the pipeline does)
    main_out=""
    for f in klebsiella_pneumo_complex_output.txt klebsiella_oxytoca_complex_output.txt escherichia_output.txt; do
        if [ -s "kleborate_out/$f" ]; then main_out="kleborate_out/$f"; break; fi
    done
    [ -n "$main_out" ] || fail "Kleborate wrote no table for $s: its species does not match the preset $preset"
    cp "$main_out" "${kout}/${s}.kleborate.txt"
    [ -s kleborate_out/klebsiella_pneumo_complex_hAMRonization_output.txt ] && \
        cp kleborate_out/klebsiella_pneumo_complex_hAMRonization_output.txt "${kout}/${s}.kleborate_hamronization.txt"
    kspecies=$(kleborate_value "${kout}/${s}.kleborate.txt" species)
    info "Kleborate: species '${kspecies}', ST '$(kleborate_value "${kout}/${s}.kleborate.txt" ST)'"

    if [[ "$kspecies" =~ ^Klebsiella\ (pneumoniae|quasipneumoniae|variicola|quasivariicola|africana) ]]; then
        step "Kaptive (K and O loci): $s"
        out="${OUTDIR}/03_typing/kaptive"; mkdir -p "$out"; cd "$out"
        run_sh KAPTIVE "
            # Kaptive 3.3 renamed the locus-FASTA option from --fasta to --loci
            if kaptive assembly --help 2>&1 | grep -- '--loci' > /dev/null; then fasta_opt=--loci; else fasta_opt=--fasta; fi
            for locus in K O; do
                db=kpsc_\$(echo \$locus | tr 'KO' 'ko')
                mkdir -p ${s}_\${locus}_locus
                kaptive assembly \$db $asm -o ${s}.kaptive_\${locus}.tsv \$fasta_opt ${s}_\${locus}_locus
            done"
    fi

    if [ "$species" = ecoli ]; then
        step "ECTyper (O:H serotype): $s"
        out="${OUTDIR}/03_typing/ectyper/${s}_work"; mkdir -p "$out"; cd "$out"
        run ECTYPER ectyper --cores "${THREADS}" --output ./ --input "$asm"
        cp output.tsv "../${s}.tsv"

        step "EzClermont (phylogroup): $s"
        out="${OUTDIR}/03_typing/ezclermont"; mkdir -p "$out"; cd "$out"
        # EzClermont can exit non-zero when no phylogroup is assigned; keep the result if it wrote one
        run_sh EZCLERMONT "ezclermont $asm > ${s}.ezclermont.tsv 2> ${s}.ezclermont.log || [ -s ${s}.ezclermont.tsv ]"
        cat "${s}.ezclermont.tsv"
    fi
done
