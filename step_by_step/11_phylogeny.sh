#!/usr/bin/env bash
#
# Step 11 - core-genome SNP phylogeny and SNP distances per species (or per species + ST) group
#           (pipeline: PHYLOGENY -> SNIPPY_RUN / SNIPPY_CONTIGS, SNIPPY_CORE and/or SKA2, then SNIPPY_CLEAN,
#            GUBBINS, IQTREE, SNPDISTS)
#
#   bash step_by_step/11_phylogeny.sh                          # every sample in the samplesheet
#   bash step_by_step/11_phylogeny.sh S1 S2 S3                 # only these samples
#   PHYLO_METHOD=both bash step_by_step/11_phylogeny.sh        # Snippy and SKA2, then compare them
#   PHYLO_GROUPING=st bash step_by_step/11_phylogeny.sh        # one group per species + ST instead of per species
#
# Needs: step 07 (Kleborate: species and ST) and step 01 (trimmed reads) or an assembly.
#
#   11a Group     : PHYLO_GROUPING (config.sh; the pipeline's --phylo_grouping)
#                   species (default): all isolates of a Kleborate species in one group (e.g. Kpneumoniae),
#                     whatever their ST, mapped to the species reference in conf/phylo_references.csv: one SNP
#                     matrix and tree for every isolate. Transmission is only plausible within an ST, so read
#                     the close pairs with their STs.
#                   st: one group per Kleborate species + ST (e.g. Kpneumoniae_ST39); an ST-specific row in
#                     conf/phylo_references.csv wins over the species row.
#   How SNPs are found: PHYLO_METHOD (config.sh; the pipeline's --phylo_method) = snippy (default), ska2 or both
#   snippy
#   11b Snippy    : map each isolate's trimmed reads to the reference and call SNPs (assemblies: snippy --ctgs)
#   11c snippy-core : one whole-genome alignment of the group (<group>.full.aln)
#   ska2
#   11b SKA2 build: split k-mers (31 bp with the middle base left out) of every isolate of the group at once,
#                   from the trimmed reads (or the assembly). No mapping, so no per-isolate step.
#   11c ska map   : place the split k-mers on the same reference -> an alignment with the same layout as
#                   snippy-core's (the reference is added as the first sequence, as snippy-core does)
#       ska distance : SNP distances between isolates without any reference (<group>.ska_distances.tsv)
#   then, for each alignment
#   11c snippy-clean_full_aln : ambiguous characters turned into N
#   11d Gubbins   : find and mask recombined regions (imported DNA would otherwise look like many mutations),
#                   keep only the remaining SNP sites
#   11e IQ-TREE   : maximum-likelihood tree on those sites (-m GTR+G, -B 1000 ultrafast bootstrap)
#   11f snp-dists : pairwise SNP distances between isolates
#   11g (both)    : bin/compare_snp_methods.py puts Snippy's and SKA2's SNP distances side by side
# Groups with fewer than MIN_GROUP isolates (default ${MIN_GROUP:-2}) are skipped. Gubbins and IQ-TREE need
# at least 4 sequences (reference + 3 isolates): a group of 2 gets SNP distances from the core alignment only.
#
# Output (06_phylogeny/), as the pipeline's:
#   snippy: <group>/snippy/<sample>/ (snps.tab, snps.vcf, ...), <group>/<group>.full.aln
#   ska2:   <group>/ (or <group>_ska2/ with both): .ska_inputs.tsv, .skf, .full.aln, .ska_distances.tsv
#   both:   <name>.treefile, <name>.snp_distances.tsv in each folder; gubbins/<name>.* (recombination);
#           snippy_vs_ska2_pairs.tsv and snippy_vs_ska2_summary.tsv (with both)
# Look at: the SNP distance matrix (isolates < ~20 SNPs apart are candidate transmission pairs), the tree
#   (open the .treefile in FigTree or iTOL), Gubbins' recombination predictions, and with both, whether
#   Snippy and SKA2 call the same pairs close.

source "$(dirname "$0")/config.sh"

phy="${OUTDIR}/06_phylogeny"; mkdir -p "$phy"
case "${PHYLO_METHOD}" in snippy|ska2|both) ;; *) fail "PHYLO_METHOD must be snippy, ska2 or both (got '${PHYLO_METHOD}')" ;; esac
case "${PHYLO_GROUPING}" in species|st) ;; *) fail "PHYLO_GROUPING must be species or st (got '${PHYLO_GROUPING}')" ;; esac
info "SNP method: ${PHYLO_METHOD}; groups: ${PHYLO_GROUPING}"

# ---------------------------------------------------------------- 11a group isolates and pick references
kleborate_value() {   # kleborate_value TABLE species|ST   (same columns as the pipeline reads)
    awk -F'\t' -v want="$2" 'NR == 1 { for (i = 1; i <= NF; i++) {
                                  if (want == "species" && ($i == "species" || $i ~ /__species__species$/)) c = i
                                  if (want == "ST" && ($i == "ST" || $i ~ /__mlst(_achtman)?__ST$/)) c = i }
                              next }
                             NR == 2 && c { print $c }' "$1"
}
reference_for() {     # reference_for "Klebsiella pneumoniae" ST39 -> accession (ST row first, then species row)
    awk -F',' -v sp="$(echo "$1" | tr 'A-Z' 'a-z')" -v st="$2" '
        /^[[:space:]]*#/ || NF < 4 { next }
        !h { h = 1; next }
        { s = tolower($1); gsub(/^ +| +$/, "", s)
          if (index(sp, s) == 1) { if ($2 == st && !hit_st) hit_st = $4; if ($2 == "" && !hit_sp) hit_sp = $4 } }
        END { print (hit_st ? hit_st : hit_sp) }' "${REFERENCE_SHEET}"
}

step "11a Groups (Kleborate species$([ "${PHYLO_GROUPING}" = st ] && echo ' + ST')) and references"
declare -A group_of ref_of members
if [ $# -gt 0 ]; then samples=$(printf '%s\n' "$@"); else samples=$(all_samples); fi
for s in $samples; do
    kt="${OUTDIR}/03_typing/kleborate/${s}.kleborate.txt"
    if [ ! -s "$kt" ]; then warn "$s: no Kleborate result (run step 07); left out"; continue; fi
    sp=$(kleborate_value "$kt" species); st=$(kleborate_value "$kt" ST)
    if [ -z "$sp" ] || [ -z "$st" ] || [ "$st" = NA ] || [ "$st" = "-" ]; then
        warn "$s: no species/ST from Kleborate ('$sp', '$st'); left out"; continue
    fi
    [ "${PHYLO_GROUPING}" = species ] && st=""       # species group and species reference only
    acc=$(reference_for "$sp" "$st")
    if [ -z "$acc" ]; then warn "$s: '$sp' has no reference in ${REFERENCE_SHEET}; left out"; continue; fi
    read -r genus epithet _ <<< "$sp"
    g=$(echo "${genus:0:1}${epithet}${st:+_${st}}" | tr -cd 'A-Za-z0-9_-')
    group_of[$s]=$g; ref_of[$g]="${REFERENCES}/${acc}.fasta"; members[$g]="${members[$g]:-} $s"
done
for g in "${!members[@]}"; do info "$(printf '%-22s %-14s %s' "$g" "$(basename "${ref_of[$g]}" .fasta)" "${members[$g]}")"; done

# ---------------------------------------------------------------- one alignment -> tree and SNP distances
tree_and_distances() {   # tree_and_distances NAME DIR N_ISOLATES   (DIR holds NAME.full.aln)
    local name=$1 dir=$2 n=$3 dist_input
    step "11c Clean the alignment: $name"
    cd "$dir"
    run_sh SNIPPY_CLEAN "snippy-clean_full_aln ${name}.full.aln > ${name}.aln"

    if [ "$n" -ge 3 ]; then
        step "11d Gubbins (recombination): $name"
        mkdir -p "${phy}/gubbins"; cd "${phy}/gubbins"
        cp "${dir}/${name}.aln" .
        mkdir -p numba_cache_dir
        run_sh GUBBINS "export NUMBA_CACHE_DIR=./numba_cache_dir; run_gubbins.py --threads ${THREADS} ${name}.aln"

        step "11e IQ-TREE: $name"
        cd "$dir"
        run IQTREE iqtree \
            -st DNA -m GTR+G -B 1000 \
            -s "${phy}/gubbins/${name}.filtered_polymorphic_sites.fasta" \
            -pre "$name" \
            -nt AUTO -ntmax "${THREADS}" -mem "${MEM_GB}G"
        dist_input="${phy}/gubbins/${name}.filtered_polymorphic_sites.fasta"
    else
        warn "$name has $n isolates: too few for Gubbins and IQ-TREE (3 or more); SNP distances only, without recombination filtering"
        dist_input="${dir}/${name}.aln"
    fi

    step "11f snp-dists: $name"
    cd "$dir"
    run_sh SNPDISTS "snp-dists $dist_input > ${name}.snp_distances.tsv"
    expand -t 10 "${name}.snp_distances.tsv" | sed 's/^/       /'
}

# ---------------------------------------------------------------- Snippy: map each isolate, then snippy-core
snippy_group() {   # snippy_group GROUP REF SAMPLE...
    local g=$1 ref=$2; shift 2
    local ids=("$@") gdir="${phy}/${g}" s
    mkdir -p "${gdir}/snippy"
    for s in "${ids[@]}"; do
        if [ -s "${gdir}/snippy/${s}/${s}.vcf" ]; then info "$s: Snippy already done"; continue; fi
        cd "${gdir}/snippy"
        if ! is_assembly_sheet && [ -s "$(trimmed_r1 "$s")" ]; then
            step "11b Snippy (reads -> ${g} reference): $s"
            run SNIPPY_RUN snippy \
                --cpus "${THREADS}" --ram "${MEM_GB}" \
                --outdir "$s" --reference "$ref" --prefix "$s" \
                --R1 "$(trimmed_r1 "$s")" --R2 "$(trimmed_r2 "$s")" --force
        else
            step "11b Snippy (contigs -> ${g} reference): $s"
            run SNIPPY_CONTIGS snippy \
                --cpus "${THREADS}" --ram "${MEM_GB}" \
                --outdir "$s" --reference "$ref" --prefix "$s" \
                --ctgs "$(best_assembly "$s")" --force
        fi
    done

    step "11c Core alignment (snippy-core): $g (${#ids[@]} isolates)"
    cd "$gdir"
    run SNIPPY_CORE snippy-core --ref "$ref" --prefix "$g" "${ids[@]/#/snippy/}"
    tree_and_distances "$g" "$gdir" "${#ids[@]}"
}

# ---------------------------------------------------------------- SKA2: split k-mers of the whole group
ska2_group() {   # ska2_group NAME REF SAMPLE...   (NAME = group, or group_ska2 with PHYLO_METHOD=both)
    local name=$1 ref=$2; shift 2
    local ids=("$@") sdir="${phy}/${name}" s asm
    mkdir -p "$sdir"; cd "$sdir"

    step "11b SKA2 build (split k-mers, k=${SKA_K}): $name (${#ids[@]} isolates)"
    : > "${name}.ska_inputs.tsv"          # one line per isolate: name <tab> R1 <tab> R2  (or: name <tab> contigs)
    for s in "${ids[@]}"; do
        if ! is_assembly_sheet && [ -s "$(trimmed_r1 "$s")" ]; then
            printf '%s\t%s\t%s\n' "$s" "$(trimmed_r1 "$s")" "$(trimmed_r2 "$s")" >> "${name}.ska_inputs.tsv"
        else
            asm=$(best_assembly "$s"); need_file "$asm" "run step 04 (Shovill) first"
            printf '%s\t%s\n' "$s" "$asm" >> "${name}.ska_inputs.tsv"
        fi
    done
    sed 's/^/       /' "${name}.ska_inputs.tsv"
    run SKA2 ska build -k "${SKA_K}" --threads "${THREADS}" -f "${name}.ska_inputs.tsv" -o "$name"

    step "11c SKA2 map onto the ${name%_ska2} reference: $name"
    run SKA2 ska map "$ref" "${name}.skf" --threads "${THREADS}" -o "${name}.ska_map.aln"
    # ska map writes the isolates only: add the reference first (its contigs joined in file order, as ska map
    # does), so the alignment looks like snippy-core's
    {
        printf '>Reference\n'
        grep -v '^>' "$ref" | tr -d '\r\n ' | tr 'acgtn' 'ACGTN'
        printf '\n'
        cat "${name}.ska_map.aln"
    } > "${name}.full.aln"
    if [ "$(awk '/^>/ { if (seq != "") print length(seq); seq = ""; next } { seq = seq $0 } END { print length(seq) }' \
            "${name}.full.aln" | sort -u | wc -l)" -ne 1 ]; then
        fail "sequences in ${sdir}/${name}.full.aln differ in length"
    fi

    step "11c SKA2 distance (reference-free SNP distances): $name"
    run SKA2 ska distance --threads "${THREADS}" -o "${name}.ska_distances.tsv" "${name}.skf"
    expand -t 14 "${name}.ska_distances.tsv" | sed 's/^/       /'

    tree_and_distances "$name" "$sdir" "${#ids[@]}"
}

# ---------------------------------------------------------------- per group
for g in $(printf '%s\n' "${!members[@]}" | sort); do
    read -ra ids <<< "${members[$g]}"; n=${#ids[@]}
    if [ "$n" -lt "${MIN_GROUP}" ]; then warn "group $g has $n isolate(s) (< MIN_GROUP=${MIN_GROUP}); skipped"; continue; fi
    ref="${ref_of[$g]}"; need_file "$ref" "run bin/fetch_references.py (trainer) first"

    case "${PHYLO_METHOD}" in
        snippy) snippy_group "$g" "$ref" "${ids[@]}" ;;
        ska2)   ska2_group   "$g" "$ref" "${ids[@]}" ;;
        both)   snippy_group "$g" "$ref" "${ids[@]}"
                ska2_group   "${g}_ska2" "$ref" "${ids[@]}" ;;
    esac
done

# ---------------------------------------------------------------- 11g compare the two methods
if [ "${PHYLO_METHOD}" = both ]; then
    step "11g Snippy vs SKA2 SNP distances"
    if command -v python3 > /dev/null; then
        python3 "${PIPELINE_DIR}/bin/compare_snp_methods.py" --phylo "$phy" \
            || warn "no comparison (see above); it needs groups with Snippy and SKA2 SNP distances"
    else
        warn "python3 not found: run bin/compare_snp_methods.py --phylo $phy yourself"
    fi
fi
