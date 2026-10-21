[← Day 1](day01.md) · [Course home](README.md) · **Day 2** · [Day 3 →](day03.md)

# Day 2: Genome assembly

## At a glance

| Session | Time | Topic | Module |
| --- | --- | --- | --- |
| 1 | 09:00–10:30 | How assemblers work; short-read assembly with Shovill | 2 |
| 2 | 11:00–12:30 | Assembly quality: QUAST, CheckM2, Bandage | 2 |
| 3 | 13:30–15:00 | Long-read assembly with Autocycler | 2 |
| 4 | 15:30–17:00 | Polishing to a hybrid assembly; comparing all three | 2 |

**You will produce:** an assembly QC table with a fit-for-purpose decision for every isolate.

## Learning objectives

By the end of today you can:

- Explain why repeats break short-read assemblies into many contigs.
- Assemble Illumina reads with Shovill, and Nanopore reads with Autocycler.
- Polish a long-read assembly with Medaka, then with Illumina reads (Polypolish and Pypolca).
- Judge an assembly with QUAST, CheckM2 and Bandage, and decide what it is good enough for.

---

## Module 2: Genome assembly

### 2.1 Background in one table

| | Short-read assembly | Long-read assembly | Hybrid (long-read first, short-read polished) |
| --- | --- | --- | --- |
| Reads | Illumina, 150 bp | Nanopore, 1–100+ kb | Both |
| Typical result for *K. pneumoniae* | 30–100 contigs | Chromosome + plasmids, often circular | Chromosome + plasmids, circular, base-accurate |
| Good for | MLST, gene detection, SNPs | Structure: which genes sit on which plasmid | Everything, including plasmid tracking |
| Weakness | Plasmids broken at repeats (IS elements) | Small plasmids can be missed; base errors before polishing | Needs both data types |

Expected genome size: *K. pneumoniae* about **5.2–5.8 Mb** (GC about 57%); *E. coli* about **4.6–5.5 Mb**
(GC about 51%).

### 2.2 Short-read assembly with Shovill

```bash
use_tools assembly
mkdir -p $RESULTS/02_assembly/shovill

# One sample first (with 1 thread, allow roughly an hour; start it before a break). CHI11E is also one of the hybrid
# isolates, so you can compare this assembly with its hybrid assembly this afternoon.
s=CHI11E
shovill --R1 $COURSE/trimmed/${s}_R1.fastq.gz --R2 $COURSE/trimmed/${s}_R2.fastq.gz \
    --outdir $RESULTS/02_assembly/shovill/$s --cpus $THREADS --minlen 200
grep -c ">" $RESULTS/02_assembly/shovill/$s/contigs.fa
```

While it runs, copy the trainers' Shovill assemblies for the other isolates (yours are kept). If you have time
later, assemble `BHB7VE` and `CAA65Z` too: the rest of your practice set.

```bash
use_precomputed 02_assembly/shovill
ls $RESULTS/02_assembly/shovill | wc -l     # 29
```

### 2.3 Assembly quality: QUAST and CheckM2

```bash
mkdir -p $COURSE/assemblies
for d in $RESULTS/02_assembly/shovill/*/; do
    s=$(basename $d)
    cp $d/contigs.fa $COURSE/assemblies/${s}.fasta
done

quast.py $COURSE/assemblies/*.fasta -o $RESULTS/02_assembly/quast_shovill -t $THREADS
# One row per assembly: contigs, total length, N50 and GC are the columns to look at
column -t -s $'\t' $RESULTS/02_assembly/quast_shovill/transposed_report.tsv | less -S

# CheckM2 on your practice isolates (a few minutes each with 1 thread)
checkm2 predict --input $(for s in $PRACTICE; do echo $COURSE/assemblies/$s.fasta; done) -x fasta \
    --output-directory $RESULTS/02_assembly/checkm2_practice --threads $THREADS \
    --database_path $SHARED/databases/checkm2/CheckM2_database/uniref100.KO.1.dmnd
cut -f1-3 $RESULTS/02_assembly/checkm2_practice/quality_report.tsv | column -t -s $'\t'

# The trainers' table for all 29: check your three rows, then look for problems
use_precomputed 02_assembly/checkm2
cut -f1-3 $RESULTS/02_assembly/checkm2/quality_report.tsv | sort -t$'\t' -k3,3gr | head
```

> [!WARNING]
> Red flags in a short-read assembly: **hundreds of contigs**, a **total size far outside the expected range**,
> or **CheckM2 contamination above 5%**. Any of these makes SNP calling and gene detection unreliable.

### 2.4 Long-read assembly with Autocycler

Autocycler builds several assemblies from subsets of the reads, with different assemblers, then takes a consensus.
This removes errors that only one assembler makes.

```bash
use_tools assembly          # Autocycler runs Flye and Raven from inside its own image
s=CHI11E
mkdir -p $RESULTS/02_assembly/autocycler/$s && cd $RESULTS/02_assembly/autocycler/$s
reads=$COURSE/trimmed/${s}_ont.fastq.gz

genome_size=$(autocycler helper genome_size --reads $reads --threads $THREADS)
echo "Estimated genome size: $genome_size"

autocycler subsample --reads $reads --out_dir subsampled_reads --genome_size $genome_size

mkdir -p assemblies
for assembler in flye raven; do
    for i in 01 02 03 04; do
        autocycler helper $assembler --reads subsampled_reads/sample_$i.fastq \
            --out_prefix assemblies/${assembler}_$i --threads $THREADS --genome_size $genome_size
    done
done

autocycler compress -i assemblies -a autocycler_out
autocycler cluster -a autocycler_out
for c in autocycler_out/clustering/qc_pass/cluster_*; do
    autocycler trim -c $c
    autocycler resolve -c $c
done
autocycler combine -a autocycler_out -i autocycler_out/clustering/qc_pass/cluster_*/5_final.gfa

grep ">" autocycler_out/consensus_assembly.fasta
cd $COURSE
```

> [!NOTE]
> This builds 8 assemblies: with 1 thread it takes several hours, longer than the session. Start it, read 2.5
> while it runs. If it has not finished when you reach 2.5, stop it (`Ctrl+C`) and fill in the trainers'
> result: `use_precomputed 02_assembly/autocycler` adds the consensus assembly yours had not written yet, and
> keeps everything your run did write.

### 2.5 Polishing: Medaka, Polypolish, Pypolca

```bash
s=CHI11E
draft=$RESULTS/02_assembly/autocycler/$s/autocycler_out/consensus_assembly.fasta
out=$RESULTS/02_assembly/hybrid/$s
mkdir -p $out

# 1. Long-read polishing
use_tools medaka
medaka_consensus -i $COURSE/trimmed/${s}_ont.fastq.gz -d $draft -o $out/medaka -t $THREADS

# 2. Short-read polishing: align each read file separately, keep all alignments (-a)
use_tools assembly
cp $out/medaka/consensus.fasta $out/medaka.fasta
bwa index $out/medaka.fasta
bwa mem -t $THREADS -a $out/medaka.fasta $COURSE/trimmed/${s}_R1.fastq.gz > $out/aln_1.sam
bwa mem -t $THREADS -a $out/medaka.fasta $COURSE/trimmed/${s}_R2.fastq.gz > $out/aln_2.sam
polypolish filter --in1 $out/aln_1.sam --in2 $out/aln_2.sam --out1 $out/filt_1.sam --out2 $out/filt_2.sam
polypolish polish $out/medaka.fasta $out/filt_1.sam $out/filt_2.sam > $out/polypolish.fasta

pypolca run -a $out/polypolish.fasta \
    -1 $COURSE/trimmed/${s}_R1.fastq.gz -2 $COURSE/trimmed/${s}_R2.fastq.gz \
    -t $THREADS -o $out/pypolca --careful
cp $out/pypolca/pypolca_corrected.fasta $COURSE/assemblies/${s}.hybrid.fasta
rm $out/*.sam
```

> [!TIP]
> Medaka 2 can use a bacteria-specific model if your basecalling model supports it. Check
> `medaka_consensus -h` for the `--bacteria` option.

### 2.6 Compare short-read, long-read and hybrid

```bash
s=CHI11E
quast.py $COURSE/assemblies/${s}.fasta \
         $RESULTS/02_assembly/autocycler/$s/autocycler_out/consensus_assembly.fasta \
         $COURSE/assemblies/${s}.hybrid.fasta \
         -l "shovill,autocycler,hybrid" -o $RESULTS/02_assembly/quast_compare_$s -t $THREADS
cat $RESULTS/02_assembly/quast_compare_$s/report.txt
```

Now look at the assembly graphs in **Bandage**. The server has no screen, so Bandage draws each graph as an
image; copy the PNG files to your computer to look at them (or open the `.gfa` files in Bandage on your own
computer):

```bash
Bandage image $RESULTS/02_assembly/shovill/$s/contigs.gfa $RESULTS/02_assembly/${s}_shovill.png
Bandage image $RESULTS/02_assembly/autocycler/$s/autocycler_out/consensus_assembly.gfa $RESULTS/02_assembly/${s}_autocycler.png
```

Fill in this table for `CHI11E`:

| Assembly | Contigs | Largest contig (bp) | Total length (bp) | Circular contigs |
| --- | --- | --- | --- | --- |
| Shovill | | | | |
| Autocycler | | | | |
| Hybrid | | | | |

### Check your understanding

<details>
<summary><b>1. The study closed <code>CHI11E</code>'s genome. How does your hybrid assembly compare?</b></summary>

Published closed assembly (Table S6 of the study):

| Contig | Length (bp) | Plasmid replicons |
| --- | --- | --- |
| Chromosome | 5,273,539 | none |
| Plasmid | 143,302 | IncFII(K), IncR |
| Plasmid | 11,493 | none |
| Plasmid | 4,269 | none |

Your hybrid should have one chromosome of about 5.27 Mb and the larger plasmids. Small differences of a few base
pairs are expected. If you are missing the 4.3 kb plasmid, that is a known weakness of long-read assembly: small
plasmids are under-represented in long reads.

</details>

<details>
<summary><b>2. And for the other two hybrid isolates?</b></summary>

| Isolate | Chromosome (bp) | Plasmids (bp, replicons) |
| --- | --- | --- |
| `CNS1F1` | 5,301,730 | 187,438 IncFIB(K)(pCAV1099-114); 142,041 IncFII(K), IncR; 16,361 |
| `1013757` | 5,315,055 | 140,482 IncFIB(AP001918), IncFII; 81,013 IncY; 10,638 |

</details>

<details>
<summary><b>3. In the Shovill assembly, why might a resistance gene sit at the very end of a contig?</b></summary>

Resistance genes are often flanked by insertion sequences such as IS26. IS copies are repeats: the assembler
cannot tell which copy joins to which, so it breaks the contig there. That is why short reads alone rarely tell
you which plasmid carries a gene.

</details>

<details>
<summary><b>4. What do you expect for <code>D38119</code> and <code>CES13U</code>?</b></summary>

From the study's tables: `D38119` assembled to **9.46 Mb** (about twice a *K. pneumoniae* genome, likely a mixed
culture) and `CES13U` to **1,468 contigs**. Neither is fit for SNP-based transmission analysis. Compare with your
Day 1 Kraken2 and depth results: do they explain it?

</details>

---

## Assembly decision

Add a column to your QC table from Day 1:

| Sample | Contigs | Total size (Mb) | CheckM2 completeness / contamination (%) | Fit for |
| --- | --- | --- | --- | --- |
| CAA65Z | | | | typing / SNPs / plasmids / exclude |

## Before you leave: the best-assembly folder

From Day 3 onwards you use one assembly per isolate: the hybrid where there is one, otherwise Shovill.

```bash
# Hybrid assemblies for the other two isolates, if you did not make them
cp -n $PRECOMPUTED/02_assembly/hybrid/*.hybrid.fasta $COURSE/assemblies/

mkdir -p $COURSE/best
for f in $COURSE/assemblies/*.fasta; do
    s=$(basename $f .fasta); s=${s%.hybrid}
    if [ -f $COURSE/assemblies/${s}.hybrid.fasta ]; then
        cp $COURSE/assemblies/${s}.hybrid.fasta $COURSE/best/${s}.fasta
    else
        cp $COURSE/assemblies/${s}.fasta $COURSE/best/${s}.fasta
    fi
done
ls $COURSE/best | wc -l    # 29
```

## End-of-day checklist

- [ ] `$COURSE/assemblies/` has a Shovill assembly for all 29 isolates (yours or pre-computed).
- [ ] `$COURSE/assemblies/` has hybrid assemblies for `CHI11E`, `CNS1F1` and `1013757`.
- [ ] `$COURSE/best/` has 29 assemblies.
- [ ] You filled in the three-way comparison table for `CHI11E`.
- [ ] Your QC table has an assembly decision for every isolate.

## Further reading

- [Shovill](https://github.com/tseemann/shovill), [Autocycler](https://github.com/rrwick/Autocycler/wiki),
  [Medaka](https://github.com/nanoporetech/medaka), [Polypolish](https://github.com/rrwick/Polypolish/wiki),
  [QUAST](https://github.com/ablab/quast), [Bandage](https://rrwick.github.io/Bandage/)
- Wick RR, Howden BP, Stinear TP. Autocycler: long-read consensus assembly for bacterial genomes. *Bioinformatics* 2025.

---

[← Day 1](day01.md) · [Course home](README.md) · **Day 2** · [Day 3 →](day03.md)
