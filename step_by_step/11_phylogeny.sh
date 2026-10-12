#!/usr/bin/env bash
#
# Step 11 - core-genome SNP phylogeny and SNP distances per species + ST group
#           (pipeline: PHYLOGENY -> SNIPPY_RUN / SNIPPY_CONTIGS, SNIPPY_CORE, SNIPPY_CLEAN, GUBBINS, IQTREE, SNPDISTS)
#
#   bash step_by_step/11_phylogeny.sh            # every sample in the samplesheet (a phylogeny needs several)
#   bash step_by_step/11_phylogeny.sh S1 S2 S3   # only these samples
#
# Needs: step 07 (Kleborate: species and ST) and step 01 (trimmed reads) or an assembly.
#
#   11a Group     : isolates are grouped by Kleborate species + ST (e.g. Kpneumoniae_ST39): transmission is
#                   only plausible between isolates of the same ST. Each group gets the reference genome
#                   listed for its species in conf/phylo_references.csv (an ST-specific row wins if there is one).
#   11b Snippy    : map each isolate's trimmed reads to the reference and call SNPs (assemblies: snippy --ctgs)
#   11c snippy-core + snippy-clean_full_aln : one whole-genome alignment of the group, with ambiguous
#                   characters turned into N
#   11d Gubbins   : find and mask recombined regions (imported DNA would otherwise look like many mutations),
#                   keep only the remaining SNP sites
#   11e IQ-TREE   : maximum-likelihood tree on those sites (-m GTR+G, -B 1000 ultrafast bootstrap)
#   11f snp-dists : pairwise SNP distances between isolates
# Groups with fewer than MIN_GROUP isolates (default ${MIN_GROUP:-2}) are skipped. Gubbins and IQ-TREE need
# at least 4 sequences (reference + 3 isolates): a group of 2 gets SNP distances from the core alignment only.
#
# Output (06_phylogeny/): <group>/snippy/<sample>/ (snps.tab, snps.vcf, ...), <group>/<group>.full.aln,
#   gubbins/<group>.* (recombination), <group>/<group>.treefile, <group>/<group>.snp_distances.tsv
# Look at: the SNP distance matrix (isolates < ~20 SNPs apart are candidate transmission pairs), the tree
#   (open the .treefile in FigTree or iTOL), and Gubbins' recombination predictions.

source "$(dirname "$0")/config.sh"

phy="${OUTDIR}/06_phylogeny"; mkdir -p "$phy"

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

step "11a Groups (Kleborate species + ST) and references"
declare -A group_of ref_of members
if [ $# -gt 0 ]; then samples=$(printf '%s\n' "$@"); else samples=$(all_samples); fi
for s in $samples; do
    kt="${OUTDIR}/03_typing/kleborate/${s}.kleborate.txt"
    if [ ! -s "$kt" ]; then warn "$s: no Kleborate result (run step 07); left out"; continue; fi
    sp=$(kleborate_value "$kt" species); st=$(kleborate_value "$kt" ST)
    if [ -z "$sp" ] || [ -z "$st" ] || [ "$st" = NA ] || [ "$st" = "-" ]; then
        warn "$s: no species/ST from Kleborate ('$sp', '$st'); left out"; continue
    fi
    acc=$(reference_for "$sp" "$st")
    if [ -z "$acc" ]; then warn "$s: '$sp' has no reference in ${REFERENCE_SHEET}; left out"; continue; fi
    read -r genus epithet _ <<< "$sp"
    g=$(echo "${genus:0:1}${epithet}_${st}" | tr -cd 'A-Za-z0-9_-')
    group_of[$s]=$g; ref_of[$g]="${REFERENCES}/${acc}.fasta"; members[$g]="${members[$g]:-} $s"
done
for g in "${!members[@]}"; do info "$(printf '%-22s %-14s %s' "$g" "$(basename "${ref_of[$g]}" .fasta)" "${members[$g]}")"; done

# ---------------------------------------------------------------- per group
for g in $(printf '%s\n' "${!members[@]}" | sort); do
    read -ra ids <<< "${members[$g]}"; n=${#ids[@]}
    if [ "$n" -lt "${MIN_GROUP}" ]; then warn "group $g has $n isolate(s) (< MIN_GROUP=${MIN_GROUP}); skipped"; continue; fi
    ref="${ref_of[$g]}"; need_file "$ref" "run bin/fetch_references.py (trainer) first"
    gdir="${phy}/${g}"; mkdir -p "${gdir}/snippy"

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

    step "11c Core alignment: $g ($n isolates)"
    cd "$gdir"
    run SNIPPY_CORE snippy-core --ref "$ref" --prefix "$g" "${ids[@]/#/snippy/}"
    run_sh SNIPPY_CLEAN "snippy-clean_full_aln ${g}.full.aln > ${g}.aln"

    if [ "$n" -ge 3 ]; then
        step "11d Gubbins (recombination): $g"
        mkdir -p "${phy}/gubbins"; cd "${phy}/gubbins"
        cp "${gdir}/${g}.aln" .
        mkdir -p numba_cache_dir
        run_sh GUBBINS "export NUMBA_CACHE_DIR=./numba_cache_dir; run_gubbins.py --threads ${THREADS} ${g}.aln"

        step "11e IQ-TREE: $g"
        cd "$gdir"
        run IQTREE iqtree \
            -m GTR+G -B 1000 \
            -s "${phy}/gubbins/${g}.filtered_polymorphic_sites.fasta" \
            -pre "$g" \
            -nt AUTO -ntmax "${THREADS}" -mem "${MEM_GB}G"
        dist_input="${phy}/gubbins/${g}.filtered_polymorphic_sites.fasta"
    else
        warn "group $g has $n isolates: too few for Gubbins and IQ-TREE (3 or more); SNP distances only, without recombination filtering"
        dist_input="${gdir}/${g}.aln"
    fi

    step "11f snp-dists: $g"
    cd "$gdir"
    run_sh SNPDISTS "snp-dists $dist_input > ${g}.snp_distances.tsv"
    expand -t 10 "${g}.snp_distances.tsv" | sed 's/^/       /'
done
