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
#   snippy   06_phylogeny/ST39/snippy/              -  (the pipeline uses another reference)  fastp + Snippy, 12 ST39 isolates
#   kaptive  07_antigens/                           per-sample tables, merged            -
# A step whose pipeline result is missing for a sample runs the course's command for that sample instead.
# Expect about 2-3 hours with 16 CPUs (Snippy and CheckM2 take most of it).
#
#SBATCH --job-name=amr-precomputed
#SBATCH --partition=batch
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/%u_precomputed_%j.txt

set -o pipefail

SHARED="${SHARED:?set SHARED to the course shared folder, e.g. export SHARED=/var/scratch/global/jjuma/neonatal-amr/shared}"
P="${PRECOMPUTED:-$SHARED/precomputed}"                                      # what participants copy from
WORK="${WORK:-/var/scratch/global/$USER/neonatal-amr/precompute_work}"      # trimmed reads, links: not shared
PIPE="${PIPE:-$P/pipeline_results}"                                          # the training pipeline run (section 8.1)
SHEET="${SHEET:-$SHARED/samplesheet.csv}"
STEPS="${STEPS:-reads fastp nanopore kraken2 assembly checkm2 kleborate amr plasmids snippy kaptive}"
THREADS="${SLURM_CPUS_PER_TASK:-${THREADS:-4}}"
HYBRIDS="CHI11E CNS1F1 1013757"
ST39="CAA65Z CAA99X CAA47J CAA4X5 CAA574 CAA5XX CAA5YM CAA4GD D29438 BKR8JH BHB7VE CHI11E"   # as on Day 6
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
        cp -n "$f" "$P/02_assembly/shovill/$s/contigs.fa"
        [ -f "$R/shovill/$s.contigs.gfa" ] && cp -n "$R/shovill/$s.contigs.gfa" "$P/02_assembly/shovill/$s/contigs.gfa"
    done
    for f in "$R"/autocycler/*.autocycler.fasta; do
        [ -e "$f" ] || continue
        s=$(basename "$f" .autocycler.fasta)
        mkdir -p "$P/02_assembly/autocycler/$s/autocycler_out"
        cp -n "$f" "$P/02_assembly/autocycler/$s/autocycler_out/consensus_assembly.fasta"
        [ -f "$R/autocycler/$s.autocycler.gfa" ] && cp -n "$R/autocycler/$s.autocycler.gfa" "$P/02_assembly/autocycler/$s/autocycler_out/consensus_assembly.gfa"
    done
    mkdir -p "$P/02_assembly/hybrid" && cp -n "$R"/hybrid/*.hybrid.fasta "$P/02_assembly/hybrid/" 2>/dev/null
fi

# Assemblies named as on the course pages: $COURSE/assemblies/<s>.fasta (Shovill), <s>.hybrid.fasta, and best/
for d in "$P"/02_assembly/shovill/*/; do
    [ -s "$d/contigs.fa" ] || continue
    s=$(basename "$d"); cp -n "$d/contigs.fa" "$WORK/assemblies/$s.fasta"
done
cp -n "$P"/02_assembly/hybrid/*.hybrid.fasta "$WORK/assemblies/" 2>/dev/null
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
if want snippy; then
    step "Snippy, ST39 against the CHI11E chromosome"
    D=$P/06_phylogeny/ST39
    mkdir -p "$D/snippy"
    ref=$D/ref_CHI11E_chromosome.fasta
    if [ ! -s "$ref" ]; then
        seqkit sort -l -r "$P/02_assembly/hybrid/CHI11E.hybrid.fasta" | seqkit head -n 1 > "$ref" || fail "snippy reference"
    fi
    for s in $ST39; do
        [ -s "$D/snippy/$s/snps.tab" ] && continue
        ensure_trimmed "$s" || { fail "snippy:$s (fastp)"; continue; }
        snippy --outdir "$D/snippy/$s" --ref "$ref" \
            --R1 "$WORK/trimmed/${s}_R1.fastq.gz" --R2 "$WORK/trimmed/${s}_R2.fastq.gz" \
            --cpus "$THREADS" --force && ran_here || fail "snippy:$s"
        rm -f "$D/snippy/$s"/*.bam "$D/snippy/$s"/*.bam.bai      # large; snippy-core does not need them
    done
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
