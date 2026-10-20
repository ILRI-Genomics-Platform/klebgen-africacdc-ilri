#!/usr/bin/env bash
#
# make_precomputed.sh - the trainers' results for every sample, in the same folders as a participant's $RESULTS.
#
# Participants run each slow step on a few isolates ($PRACTICE) and copy the rest with `use_precomputed <folder>`
# (course/tools.sh). This script fills those folders under $SHARED/precomputed. Where the training pipeline run
# (TRAINER_GUIDE section 8.1) already produced the same result, with the same tool, options and database, it is
# reused and only renamed to the course's file names. Everything else is made by running the course page's own
# command, so a copied result is what the participant's command would have produced.
#
#   cd $SHARED/neonatal-amr-nf && mkdir -p logs
#   sbatch course/make_precomputed.sh                                         # every step
#   sbatch --export=ALL,STEPS="checkm2 snippy" course/make_precomputed.sh     # only these steps
#
# Needs: the training pipeline run, in $SHARED/precomputed/pipeline_results (set PIPE for another folder).
# Steps already done are skipped, so the script can simply be submitted again after a failure.
#
# Step       Folder under $SHARED/precomputed       From the pipeline run                Run here
#   reads    01_qc/read_stats.tsv                   -                                    seqkit stats (fast)
#   fastp    01_qc/fastp/                           reports (same options), renamed      -
#   nanopore 01_qc/nanoplot/                        -                                    Filtlong + NanoPlot, 3 hybrid isolates
#   kraken2  01_qc/kraken2/                         Kraken2 reports (same database)      Bracken with the course's -r 150 (seconds)
#   assembly 02_assembly/{shovill,autocycler,hybrid}/  assemblies                        -
#   checkm2  02_assembly/checkm2/                   -                                    CheckM2, every Shovill assembly
#   kleborate 03_typing/kleborate/                  per-sample tables, merged            -
#   amr      04_amr/amrfinderplus/                  per-sample tables, + Name column     -
#   plasmids 05_plasmids/{hybrid,shortread}/        MOB-suite on each best assembly      MOB-suite, Shovill assembly of 3 hybrid isolates
#   snippy   06_phylogeny/ST39/snippy/              -  (it keeps no Snippy alignments)  fastp + Snippy, 12 ST39 isolates, against
#                                                                                        the reference for their Kleborate species
#                                                                                        (conf/phylo_references.csv: NTUH-K2044)
#   phylo    06_phylogeny/ST39/                     -                                    snippy-core, Gubbins, IQ-TREE, snp-dists (Day 6, 6.4-6.7)
#   ska2     06_phylogeny/ST39/ST39_ska*            -  (it used reads, the course uses assemblies)  ska build, ska distance (Day 6, 6.8)
#   kaptive  07_antigens/                           per-sample tables, merged            -
# A step whose pipeline result is missing for a sample runs the course's command for that sample instead.
# Expect about 2-3 hours with 16 CPUs (Snippy, CheckM2 and Gubbins take most of it).
#
#SBATCH --job-name=amr-precomputed
#SBATCH --partition=batch
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/%u_precomputed_%j.txt

set -o pipefail

SHARED="${SHARED:?set SHARED to the course shared folder, e.g. export SHARED=/var/scratch/global/jjuma/neonatal-amr/shared}"
P="${PRECOMPUTED:-$SHARED/precomputed}"                                      # what participants copy from
WORK="${WORK:-/var/scratch/global/$USER/neonatal-amr/precompute_work}"      # trimmed reads, links: not shared
PIPE="${PIPE:-$P/pipeline_results}"                                          # the training pipeline run (section 8.1)
SHEET="${SHEET:-$SHARED/samplesheet.csv}"
STEPS="${STEPS:-reads fastp nanopore kraken2 assembly checkm2 kleborate amr plasmids snippy phylo ska2 kaptive}"
THREADS="${SLURM_CPUS_PER_TASK:-${THREADS:-4}}"
HYBRIDS="CHI11E CNS1F1 1013757"
ST39="CAA65Z CAA99X CAA47J CAA4X5 CAA574 CAA5XX CAA5YM CAA4GD D29438 BKR8JH BHB7VE CHI11E"   # as on Day 6
CLUSTER="CAA65Z CAA99X CAA47J CAA4X5 CAA574 CAA5XX CAA5YM CAA4GD"                            # as on Day 6, 6.7
REF_SHEET=$SHARED/neonatal-amr-nf/conf/phylo_references.csv
KDB=$SHARED/databases/kraken2
AMRDB=$SHARED/databases/amrfinderplus/latest

# The course pages' variables, pointing at the shared result folders
export RESULTS=$P COURSE=$WORK READS=$WORK/reads THREADS
mkdir -p "$P" "$READS" "$WORK/trimmed" "$WORK/fastp" "$WORK/assemblies" "$WORK/best"

source "$SHARED/neonatal-amr-nf/course/tools.sh" || { echo "ERROR: cannot load course/tools.sh" >&2; exit 1; }
use_tools all >/dev/null 2>&1 || true                                         # missing images show up per step

step() { printf '\n==== %s  (%s)\n' "$*" "$(date '+%F %T')"; }
want() { [[ " $STEPS " == *" $1 "* ]]; }
failed=""
fail() { echo "FAILED: $*" >&2; failed="$failed $1"; }
reused=0; ran=0
from_pipeline() { reused=$((reused + 1)); }
ran_here()      { ran=$((ran + 1)); }

# cpn SRC... DEST : copy each SRC unless DEST (a file, or DEST/<name> for a folder) exists already
cpn() {
    local dest=${!#} f t
    for f in "${@:1:$#-1}"; do
        [ -e "$f" ] || continue
        if [ -d "$dest" ]; then t=$dest/$(basename "$f"); else t=$dest; fi
        [ -e "$t" ] || cp "$f" "$t"
    done
}

# Trimmed reads, only where a step needs them: fastp with the course's options, reports kept out of $P
ensure_trimmed() {   # ensure_trimmed SAMPLE
    local s=$1
    [ -s "$WORK/trimmed/${s}_R1.fastq.gz" ] && [ -s "$WORK/trimmed/${s}_R2.fastq.gz" ] && return 0
    fastp -i "$READS/${s}_R1.fastq.gz" -I "$READS/${s}_R2.fastq.gz" \
          -o "$WORK/trimmed/${s}_R1.fastq.gz" -O "$WORK/trimmed/${s}_R2.fastq.gz" \
          --qualified_quality_phred 20 --length_required 50 \
          -j "$WORK/fastp/${s}.fastp.json" -h "$WORK/fastp/${s}.fastp.html" \
          -w "$THREADS" > "$WORK/fastp/${s}.log" 2>&1
}

# merge_tables OUT FILE... : one table from per-sample tables with the same columns. The first column (the
# assembly name, which the pipeline takes from its own file names) becomes the sample name, as on the course
# pages, where each assembly is called <sample>.fasta.
merge_tables() {
    local out=$1; shift
    awk -F'\t' 'BEGIN { OFS = "\t" }
        FNR == 1 { s = FILENAME; sub(/.*\//, "", s); sub(/\..*/, "", s); if (NR == 1) print; next }
        { $1 = s; print }' "$@" > "$out"
}

# Samples and read links, as on Day 1
samples=$(tail -n +2 "$SHEET" | cut -d, -f1)
tail -n +2 "$SHEET" | while IFS=, read -r sample species r1 r2 ont group; do
    ln -sf "$r1" "$READS/${sample}_R1.fastq.gz"
    ln -sf "$r2" "$READS/${sample}_R2.fastq.gz"
    [ -n "$ont" ] && ln -sf "$ont" "$READS/${sample}_ont.fastq.gz"
done
echo "Samples: $(echo $samples | wc -w); results in $P; pipeline run in $PIPE; scratch in $WORK; $THREADS threads"
[ -d "$PIPE" ] || echo "WARNING: no pipeline run in $PIPE: every step will be run here (much longer)" >&2

# ------------------------------------------------------------------------- Day 1
if want reads && [ ! -s "$P/01_qc/read_stats.tsv" ]; then
    step "seqkit stats"
    mkdir -p "$P/01_qc"
    seqkit stats -a -T -j "$THREADS" "$READS"/*.fastq.gz > "$P/01_qc/read_stats.tsv" && ran_here || fail reads
fi

if want fastp; then
    step "fastp reports"
    mkdir -p "$P/01_qc/fastp"
    for s in $samples; do
        [ -s "$P/01_qc/fastp/${s}.fastp.json" ] && continue
        if [ -s "$PIPE/01_qc/fastp/${s}.fastp.json" ]; then           # same options (conf/modules.config)
            cp "$PIPE/01_qc/fastp/${s}.fastp.json" "$P/01_qc/fastp/${s}.fastp.json"
            [ -f "$PIPE/01_qc/fastp/${s}.fastp.html" ] && cp "$PIPE/01_qc/fastp/${s}.fastp.html" "$P/01_qc/fastp/${s}.fastp.html"
            [ -f "$PIPE/01_qc/fastp/${s}.fastp.log" ]  && cp "$PIPE/01_qc/fastp/${s}.fastp.log"  "$P/01_qc/fastp/${s}.log"
            from_pipeline
        else
            ensure_trimmed "$s" || { fail "fastp:$s"; continue; }
            cp "$WORK/fastp/${s}.fastp.json" "$WORK/fastp/${s}.fastp.html" "$WORK/fastp/${s}.log" "$P/01_qc/fastp/"
            ran_here
        fi
    done
fi

if want nanopore; then
    step "NanoPlot and Filtlong (hybrid isolates)"
    mkdir -p "$P/01_qc/nanoplot"
    for s in $HYBRIDS; do
        [ -e "$READS/${s}_ont.fastq.gz" ] || continue
        [ -s "$P/01_qc/nanoplot/${s}_filtered/NanoStats.txt" ] && continue
        NanoPlot --fastq "$READS/${s}_ont.fastq.gz" -o "$P/01_qc/nanoplot/${s}_raw" -t "$THREADS" || fail "nanoplot:$s"
        filtlong --min_length 1000 --keep_percent 90 "$READS/${s}_ont.fastq.gz" | gzip > "$WORK/trimmed/${s}_ont.fastq.gz"
        NanoPlot --fastq "$WORK/trimmed/${s}_ont.fastq.gz" -o "$P/01_qc/nanoplot/${s}_filtered" -t "$THREADS" \
            && ran_here || fail "nanoplot:$s"
    done
fi

if want kraken2; then
    step "Kraken2 (from the pipeline) and Bracken"
    cat "$KDB/VERSION" 2>/dev/null
    mkdir -p "$P/01_qc/kraken2"
    for s in $samples; do
        k=$P/01_qc/kraken2/${s}.k2report
        if [ ! -s "$k" ]; then
            if [ -s "$PIPE/01_qc/kraken2/${s}.kraken2.report.txt" ]; then   # same shared database
                cp "$PIPE/01_qc/kraken2/${s}.kraken2.report.txt" "$k"; from_pipeline
            else
                ensure_trimmed "$s" || { fail "kraken2:$s (fastp)"; continue; }
                kraken2 --db "$KDB" --threads "$THREADS" --paired --report "$k" --output /dev/null \
                    "$WORK/trimmed/${s}_R1.fastq.gz" "$WORK/trimmed/${s}_R2.fastq.gz" && ran_here || { fail "kraken2:$s"; continue; }
            fi
        fi
        # Bracken with the course's options (the pipeline's differ: -t 10, default read length)
        [ -s "$P/01_qc/kraken2/${s}.bracken.tsv" ] && continue
        bracken -d "$KDB" -i "$k" -o "$P/01_qc/kraken2/${s}.bracken.tsv" -l S -r 150 >/dev/null || fail "bracken:$s"
    done
fi

# ------------------------------------------------------------------------- Day 2
if want assembly; then
    step "Assemblies from the pipeline run"
    R=$PIPE/02_assembly
    [ -d "$R/shovill" ] || fail "assembly (no $R/shovill: run the training pipeline first, section 8.1)"
    for f in "$R"/shovill/*.contigs.fa; do
        [ -e "$f" ] || continue
        s=$(basename "$f" .contigs.fa)
        mkdir -p "$P/02_assembly/shovill/$s"
        cpn "$f" "$P/02_assembly/shovill/$s/contigs.fa"
        [ -f "$R/shovill/$s.contigs.gfa" ] && cpn "$R/shovill/$s.contigs.gfa" "$P/02_assembly/shovill/$s/contigs.gfa"
    done
    for f in "$R"/autocycler/*.autocycler.fasta; do
        [ -e "$f" ] || continue
        s=$(basename "$f" .autocycler.fasta)
        mkdir -p "$P/02_assembly/autocycler/$s/autocycler_out"
        cpn "$f" "$P/02_assembly/autocycler/$s/autocycler_out/consensus_assembly.fasta"
        [ -f "$R/autocycler/$s.autocycler.gfa" ] && cpn "$R/autocycler/$s.autocycler.gfa" "$P/02_assembly/autocycler/$s/autocycler_out/consensus_assembly.gfa"
    done
    mkdir -p "$P/02_assembly/hybrid" && cpn "$R"/hybrid/*.hybrid.fasta "$P/02_assembly/hybrid/" 2>/dev/null
fi

# Assemblies named as on the course pages: $COURSE/assemblies/<s>.fasta (Shovill), <s>.hybrid.fasta, and best/
for d in "$P"/02_assembly/shovill/*/; do
    [ -s "$d/contigs.fa" ] || continue
    s=$(basename "$d"); cpn "$d/contigs.fa" "$WORK/assemblies/$s.fasta"
done
cpn "$P"/02_assembly/hybrid/*.hybrid.fasta "$WORK/assemblies/" 2>/dev/null
for f in "$WORK"/assemblies/*.fasta; do
    [ -e "$f" ] || continue
    s=$(basename "$f" .fasta); s=${s%.hybrid}
    if [ -f "$WORK/assemblies/$s.hybrid.fasta" ]; then cp "$WORK/assemblies/$s.hybrid.fasta" "$WORK/best/$s.fasta"
    else cp "$WORK/assemblies/$s.fasta" "$WORK/best/$s.fasta"; fi
done
best_samples=$(for f in "$WORK"/best/*.fasta; do [ -e "$f" ] && basename "$f" .fasta; done)
echo "Best assemblies: $(echo $best_samples | wc -w)"

if want checkm2 && [ ! -s "$P/02_assembly/checkm2/quality_report.tsv" ]; then
    step "CheckM2"
    shovill_asm=$(ls "$WORK"/assemblies/*.fasta | grep -v '\.hybrid\.fasta$')     # Shovill only, as on Day 2
    checkm2 predict --input $shovill_asm -x fasta \
        --output-directory "$P/02_assembly/checkm2" --threads "$THREADS" --force \
        --database_path "$SHARED/databases/checkm2/CheckM2_database/uniref100.KO.1.dmnd" && ran_here || fail checkm2
fi

# ------------------------------------------------------------------------- Day 3
if want kleborate && [ ! -s "$P/03_typing/kleborate/klebsiella_pneumo_complex_output.txt" ]; then
    step "Kleborate"
    mkdir -p "$P/03_typing/kleborate"
    tables=""; missing=""
    for s in $best_samples; do
        f=$PIPE/03_typing/kleborate/${s}.kleborate.txt
        if [ -s "$f" ]; then tables="$tables $f"; else missing="$missing $s"; fi
    done
    if [ -z "$missing" ] && [ -n "$tables" ]; then                  # the pipeline ran Kleborate on the same best assemblies
        merge_tables "$P/03_typing/kleborate/klebsiella_pneumo_complex_output.txt" $tables && from_pipeline || fail kleborate
    else
        [ -n "$missing" ] && echo "No pipeline Kleborate table for:$missing; running Kleborate on every best assembly"
        kleborate -a "$WORK"/best/*.fasta -o "$P/03_typing/kleborate" -p kpsc && ran_here || fail kleborate
    fi
fi

if want amr; then
    step "AMRFinderPlus"
    mkdir -p "$P/04_amr/amrfinderplus"
    for s in $best_samples; do
        out=$P/04_amr/amrfinderplus/${s}.tsv
        [ -s "$out" ] && continue
        f=$PIPE/04_amr/amrfinderplus/${s}.tsv
        if [ -s "$f" ]; then
            # same organism and --plus; the course adds --name, which puts the sample in a first 'Name' column
            awk -v s="$s" 'BEGIN { OFS = "\t" } NR == 1 { print "Name", $0; next } { print s, $0 }' "$f" > "$out" && from_pipeline
        else
            amrfinder -n "$WORK/best/$s.fasta" --organism Klebsiella_pneumoniae --plus --name "$s" \
                -d "$AMRDB" --threads "$THREADS" -o "$out" && ran_here || fail "amr:$s"
        fi
    done
fi

# ------------------------------------------------------------------------- Days 4-5
if want plasmids; then
    step "MOB-suite"
    pipeline_mob() {   # pipeline_mob SAMPLE -> folder with the pipeline's MOB-suite files, or nothing
        local d=$PIPE/05_plasmids/mobsuite/$1
        if [ -s "$d/results/contig_report.txt" ]; then echo "$d/results"; elif [ -s "$d/contig_report.txt" ]; then echo "$d"; fi
    }
    for s in $best_samples; do
        hyb=no; [[ " $HYBRIDS " == *" $s "* ]] && hyb=yes
        kind=shortread; [ $hyb = yes ] && kind=hybrid                  # the pipeline ran MOB-suite on the best assembly
        out=$P/05_plasmids/$kind/$s
        if [ ! -s "$out/contig_report.txt" ]; then
            src=$(pipeline_mob "$s")
            if [ -n "$src" ]; then mkdir -p "$out" && cp -r "$src"/. "$out"/ && from_pipeline
            else
                asm=$WORK/assemblies/$s.fasta; [ $hyb = yes ] && asm=$WORK/assemblies/$s.hybrid.fasta
                mob_recon --infile "$asm" --outdir "$out" -n "$THREADS" --force && ran_here || fail "mob_recon:$kind:$s"
            fi
        fi
        # Day 4 also compares the hybrid isolates' Shovill assemblies: run those here
        if [ $hyb = yes ] && [ ! -s "$P/05_plasmids/shortread/$s/contig_report.txt" ]; then
            mob_recon --infile "$WORK/assemblies/$s.fasta" --outdir "$P/05_plasmids/shortread/$s" -n "$THREADS" --force \
                && ran_here || fail "mob_recon:shortread:$s"
        fi
    done
fi

# ------------------------------------------------------------------------- Day 6
D=$P/06_phylogeny/ST39
# An earlier version of this script mapped to the CHI11E chromosome: never mix two references in one alignment
if [ -e "$D/ref_CHI11E_chromosome.fasta" ] && { want snippy || want phylo; }; then
    echo "Removing $D (made with the earlier CHI11E reference)"; rm -rf "$D"
fi
if want snippy; then
    step "Snippy, ST39 against the reference for their species"
    mkdir -p "$D/snippy"
    # Species of each ST39 isolate from Kleborate, as on Day 6 (6.2): only K. pneumoniae goes in this tree
    K=$P/03_typing/kleborate/klebsiella_pneumo_complex_output.txt
    [ -s "$K" ] || fail "snippy (no Kleborate table: run the kleborate step first)"
    kp=""
    for s in $ST39; do
        sp=$(awk -F'\t' -v s="$s" 'NR == 1 { for (i = 1; i <= NF; i++) if ($i ~ /(^|__)species$/) c = i; next }
                                   $1 == s && c { print $c }' "$K" 2>/dev/null)
        case "$sp" in
            "Klebsiella pneumoniae"*) kp="$kp $s" ;;
            "") echo "  $s: not in the Kleborate table: left out" ;;
            *)  echo "  $s: Kleborate species '$sp', not K. pneumoniae: left out" ;;
        esac
    done
    # The species-wide reference from the pipeline's sheet (an ST39-specific row would win, as in the pipeline)
    read -r ref_name ref_acc < <(awk -F',' '/^[[:space:]]*#/ || NF < 4 { next } $1 == "Klebsiella pneumoniae" && $2 == "ST39" { print $3, $4; f = 1; exit }
             $1 == "Klebsiella pneumoniae" && $2 == "" && !sp { sp = $3 " " $4 } END { if (!f && sp) print sp }' "$REF_SHEET")
    ref=$D/ref_${ref_name}.fasta
    echo "Reference: $ref_name ($ref_acc); isolates:$kp"
    if [ ! -s "$ref" ]; then
        cp "$SHARED/references/${ref_acc}.fasta" "$ref" || fail "snippy reference ($SHARED/references/${ref_acc}.fasta: run bin/fetch_references.py, section 7.4)"
    fi
    for s in $kp; do
        [ -s "$D/snippy/$s/snps.tab" ] && continue
        ensure_trimmed "$s" || { fail "snippy:$s (fastp)"; continue; }
        snippy --outdir "$D/snippy/$s" --ref "$ref" \
            --R1 "$WORK/trimmed/${s}_R1.fastq.gz" --R2 "$WORK/trimmed/${s}_R2.fastq.gz" \
            --cpus "$THREADS" --force && ran_here || fail "snippy:$s"
        rm -f "$D/snippy/$s"/*.bam "$D/snippy/$s"/*.bam.bai      # large; snippy-core does not need them
    done
fi

# The rest of Day 6, 6.4-6.7, with the page's commands and file names (in $D, like the participants)
if want phylo && [ ! -s "$D/ST39_pairs.tsv" ]; then
    step "Core alignment, Gubbins, IQ-TREE and SNP distances (ST39)"
    ref=$(ls "$D"/ref_*.fasta 2>/dev/null | head -1)
    if [ -z "$ref" ] || [ -z "$(ls -d "$D"/snippy/*/ 2>/dev/null)" ]; then
        fail "phylo (no Snippy results: run the snippy step first)"
    else
        ( cd "$D" || exit 1
          snippy-core --ref "$(basename "$ref")" --prefix core snippy/* &&
          snippy-clean_full_aln core.full.aln > clean.full.aln &&
          run_gubbins.py --prefix ST39 --threads "$THREADS" clean.full.aln &&
          snp-sites -c ST39.filtered_polymorphic_sites.fasta > ST39.core_snps.fasta &&
          iqtree -s ST39.core_snps.fasta -st DNA -m GTR+G+ASC -B 1000 -T "$THREADS" --prefix ST39 &&
          snp-dists ST39.core_snps.fasta > ST39_snp_dists_filtered.tsv &&
          snp-dists core.aln > ST39_snp_dists_raw.tsv &&
          awk -F'\t' -v c="$CLUSTER" 'BEGIN{n=split(c,a," "); for(i=1;i<=n;i++) inC[a[i]]=1}
              NR==1 {for(i=2;i<=NF;i++) h[i]=$i; next}
              {for(i=2;i<=NF;i++) if ($1 < h[i]) {
                  t = (inC[$1] && inC[h[i]]) ? "within_cluster" : ((inC[$1] || inC[h[i]]) ? "cluster_vs_other" : "other_vs_other")
                  print $1"\t"h[i]"\t"$i"\t"t }}' ST39_snp_dists_filtered.tsv > ST39_pairs.tsv
        ) && ran_here || fail phylo
    fi
fi

if want ska2 && [ ! -s "$D/ST39_ska_distances.tsv" ]; then
    step "SKA2 (ST39, from the best assemblies)"
    mkdir -p "$D"
    ( cd "$D" || exit 1
      for s in $ST39; do [ -s "$WORK/best/$s.fasta" ] && printf '%s\t%s\n' "$s" "$WORK/best/$s.fasta"; done > ska_inputs.tsv
      ska build -f ska_inputs.tsv -o ST39_ska -k 31 --threads "$THREADS" &&
      ska distance ST39_ska.skf -o ST39_ska_distances.tsv
    ) && ran_here || fail ska2
    # The input list names this script's scratch folder; participants write their own (6.8)
    rm -f "$D/ska_inputs.tsv"
fi

# ------------------------------------------------------------------------- Day 8
if want kaptive; then
    step "Kaptive"
    mkdir -p "$P/07_antigens"
    for locus in K O; do
        lc=$(echo $locus | tr 'KO' 'ko')
        out=$P/07_antigens/kaptive_${lc}.tsv
        [ -s "$out" ] && continue
        tables=""; missing=""
        for s in $best_samples; do
            f=$PIPE/03_typing/kaptive/${s}.kaptive_${locus}.tsv
            if [ -s "$f" ]; then tables="$tables $f"; else missing="$missing $s"; fi
        done
        if [ -z "$missing" ] && [ -n "$tables" ]; then
            merge_tables "$out" $tables && from_pipeline || fail "kaptive_$lc"
        else
            kaptive assembly "kpsc_$lc" "$WORK"/best/*.fasta -o "$out" -t "$THREADS" && ran_here || fail "kaptive_$lc"
        fi
    done
fi

# ------------------------------------------------------------------------- done
step "Summary"
echo "Results taken from the pipeline run: $reused; made here: $ran"
( cd "$P" && find . -mindepth 1 -maxdepth 2 -type d ! -path './pipeline_results*' | sort | while read -r d; do
      printf '  %-40s %6s files\n' "${d#./}" "$(find "$d" -type f | wc -l)"; done )
if [ -n "$failed" ]; then
    echo "Some steps failed:$failed  (fix them and submit the script again: finished steps are skipped)" >&2
    exit 1
fi
echo "All done. Make it read-only for participants:  chmod -R o+rX,go-w $P"
