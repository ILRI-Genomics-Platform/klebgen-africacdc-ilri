[← Day 5](day05.md) · [Course home](README.md) · **Day 6** · [Day 7 →](day07.md)

# Day 6: Building trees for transmission analysis

## At a glance

| Session | Time | Topic | Module |
| --- | --- | --- | --- |
| 1 | 09:00–10:30 | From reads to a core-genome SNP alignment (Snippy) | 6 |
| 2 | 11:00–12:30 | Recombination and why it matters (Gubbins) | 6 |
| 3 | 13:30–15:00 | Maximum-likelihood trees (IQ-TREE) | 6 |
| 4 | 15:30–17:00 | SNP distances and what a threshold can and cannot tell you | 6 |

**You will produce:** a recombination-filtered tree and a SNP-distance matrix for the ST39 group.

## Learning objectives

By the end of today you can:

- Choose a reference genome and build a core-genome SNP alignment.
- Explain why recombination inflates SNP counts, and remove it with Gubbins.
- Build a maximum-likelihood tree with branch support.
- Calculate pairwise SNP distances and use a threshold as **evidence**, not proof, of transmission.

---

## Module 6: Phylogenetics and outbreak detection (part 1)

### 6.1 The ST39 group

Today you analyse all 12 ST39 isolates together:

| Role | Isolates |
| --- | --- |
| Candidate cluster: Chatinkha Nursery, 2017 | `CAA65Z` `CAA99X` `CAA47J` `CAA4X5` `CAA574` `CAA5XX` `CAA5YM` `CAA4GD` |
| Background ST39 | `D29438` (2005), `BKR8JH` (2016) |
| Other ST39 | `BHB7VE` (2016, NDM-1), `CHI11E` (2019, hybrid) |

> [!IMPORTANT]
> Build trees **within one ST** (or clonal group). Between distant STs, recombination dominates and SNP counts run
> into the tens of thousands, which says nothing about transmission.

### 6.2 Choose a reference

The closer the reference, the more of the genome you can compare. A closed genome from the same ST is ideal, and
you made one on Day 2: the `CHI11E` hybrid assembly. Use its **chromosome only**, because plasmids are gained and
lost too often to follow the same history.

```bash
use_tools phylo
mkdir -p $RESULTS/06_phylogeny/ST39 && cd $RESULTS/06_phylogeny/ST39

# The chromosome is the longest sequence
seqkit sort -l -r $COURSE/assemblies/CHI11E.hybrid.fasta | seqkit head -n 1 > ref_CHI11E_chromosome.fasta
seqkit stats ref_CHI11E_chromosome.fasta
```

> [!NOTE]
> **How the pipeline chooses references.** By default it groups isolates by Kleborate's species and ST
> (`Kpneumoniae_ST39`) and uses one reference per species from `conf/phylo_references.csv`: NTUH-K2044 (ST23)
> for *K. pneumoniae*. That works for every ST without preparation, but a reference from another ST leaves out
> genes your lineage has and NTUH-K2044 lacks, so the core genome is smaller than with CHI11E. Adding an
> `ST39` row to the sheet gives the pipeline the same choice you are making by hand here.

### 6.3 Map reads and call SNPs with Snippy

```bash
ST39="CAA65Z CAA99X CAA47J CAA4X5 CAA574 CAA5XX CAA5YM CAA4GD D29438 BKR8JH BHB7VE CHI11E"

for s in $ST39; do
    snippy --outdir snippy/$s --ref ref_CHI11E_chromosome.fasta \
        --R1 $COURSE/trimmed/${s}_R1.fastq.gz --R2 $COURSE/trimmed/${s}_R2.fastq.gz \
        --cpus $THREADS --force
done

# How many SNPs did each isolate have against the reference?
for s in $ST39; do echo -e "$s\t$(grep -c -v '^CHROM' snippy/$s/snps.tab)"; done
```

> [!NOTE]
> `CHI11E` mapped to its own chromosome should give close to **zero** SNPs. If it does not, something is wrong
> with the reads, the reference or the settings: a useful built-in check.

### 6.4 Core alignment

```bash
snippy-core --ref ref_CHI11E_chromosome.fasta --prefix core snippy/*
cat core.txt                    # per-isolate: aligned length and how much is missing
snippy-clean_full_aln core.full.aln > clean.full.aln
```

<details>
<summary><b>What is the difference between <code>core.aln</code> and <code>core.full.aln</code>?</b></summary>

`core.aln` holds only the variable sites in the core genome (positions covered in every isolate).
`core.full.aln` is the whole reference-length alignment for every isolate, including invariant sites and gaps.
Gubbins needs the full alignment, because it looks for regions where SNPs are unusually dense.

</details>

### 6.5 Remove recombination with Gubbins

```bash
run_gubbins.py --prefix ST39 --threads $THREADS clean.full.aln
ls ST39.*
```

Key outputs:

| File | Contents |
| --- | --- |
| `ST39.recombination_predictions.gff` | Where recombination was detected, and on which branch |
| `ST39.filtered_polymorphic_sites.fasta` | SNP alignment with recombination removed |
| `ST39.final_tree.tre` | Gubbins' own tree |
| `ST39.per_branch_statistics.csv` | SNPs inside and outside recombination, per branch |

```bash
# How many recombination events, and how many SNPs did they carry?
grep -vc "^#" ST39.recombination_predictions.gff
awk -F'\t' 'NR>1 {inr+=$3; outr+=$4} END {print "SNPs outside recombination:", outr, " inside:", inr}' \
    ST39.per_branch_statistics.csv
```

> [!TIP]
> Despite its name, this file is tab-separated. Check the column names with
> `head -1 ST39.per_branch_statistics.csv | tr '\t' '\n' | cat -n`: the awk above assumes column 3 is SNPs
> inside recombination and column 4 is SNPs outside. Adjust the numbers if your Gubbins version differs.

### 6.6 Build the tree with IQ-TREE

```bash
snp-sites -c ST39.filtered_polymorphic_sites.fasta > ST39.core_snps.fasta

# GTR+G with ascertainment-bias correction (+ASC), because only variable sites are used;
# 1000 ultrafast bootstraps for branch support
iqtree -s ST39.core_snps.fasta -m GTR+G+ASC -B 1000 -T $THREADS --prefix ST39
ls ST39.treefile ST39.contree
```

> [!NOTE]
> On older installs the command is `iqtree2` instead of `iqtree`.

View the tree: download `ST39.treefile` and open it in [iTOL](https://itol.embl.de) or
[FigTree](https://github.com/rambaut/figtree/releases), or draw it in R on Day 9.

### 6.7 SNP distances

```bash
snp-dists ST39.core_snps.fasta > ST39_snp_dists_filtered.tsv
snp-dists core.aln           > ST39_snp_dists_raw.tsv          # before removing recombination
less -S ST39_snp_dists_filtered.tsv
```

Put the distances in long format and compare within the candidate cluster with distances to the background:

```bash
CLUSTER="CAA65Z CAA99X CAA47J CAA4X5 CAA574 CAA5XX CAA5YM CAA4GD"
awk -F'\t' -v c="$CLUSTER" 'BEGIN{n=split(c,a," "); for(i=1;i<=n;i++) inC[a[i]]=1}
    NR==1 {for(i=2;i<=NF;i++) h[i]=$i; next}
    {for(i=2;i<=NF;i++) if ($1 < h[i]) {
        t = (inC[$1] && inC[h[i]]) ? "within_cluster" : ((inC[$1] || inC[h[i]]) ? "cluster_vs_other" : "other_vs_other")
        print $1"\t"h[i]"\t"$i"\t"t }}' ST39_snp_dists_filtered.tsv > ST39_pairs.tsv

# Smallest, median and largest distance for each kind of pair
for t in within_cluster cluster_vs_other other_vs_other; do
    awk -v t=$t '$4==t {print $3}' ST39_pairs.tsv | sort -n \
      | awk -v t=$t '{a[NR]=$1} END {if (NR) print t": n="NR", min="a[1]", median="a[int((NR+1)/2)]", max="a[NR]}'
done
```

### 6.8 What does a SNP threshold mean?

> [!WARNING]
> There is no universal SNP threshold for transmission. For *K. pneumoniae*, **fewer than about 20–25 core SNPs**
> is often used to flag possible recent transmission, but the right value depends on the lineage, the mutation
> rate, how much of the genome is compared, and the time between samples. Use a threshold to **flag pairs for
> epidemiological review**, never to declare an outbreak on its own.

Questions to discuss with your neighbour:

1. Do the eight cluster isolates fall within a few SNPs of each other?
2. Are the background isolates (2005, 2016) clearly further away?
3. Where do `BHB7VE` and `CHI11E` sit?
4. How did removing recombination change the distances (`_raw` versus `_filtered`)?

### Check your understanding

<details>
<summary><b>1. Why does recombination make two isolates look less related than they are?</b></summary>

A single recombination event can import a stretch of DNA from a distant strain, bringing dozens or hundreds of
SNPs at once. Those SNPs come from one event, not from gradual mutation over time. Counting them as mutations
makes the two isolates look as if they diverged long ago.

</details>

<details>
<summary><b>2. Suppose the cluster isolates are 0–5 SNPs apart. Does that prove they were transmitted from one baby to another?</b></summary>

No. It shows they share a very recent common ancestor. Transmission could have been direct (baby to baby, via
hands) or indirect, from a shared source such as a sink, an incubator or contaminated feeds. Distinguishing these
needs epidemiology: ward stays, dates, and environmental sampling. That is tomorrow.

</details>

<details>
<summary><b>3. Why use only the chromosome as reference?</b></summary>

Plasmids are gained, lost and exchanged between lineages. SNPs on plasmids do not follow the same history as the
chromosome, and plasmid regions missing from some isolates reduce the core genome.

</details>

---

## On your own: repeat for ST35

Build the ST35 tree yourself: `CAAGNM CAAGJF CAAGM3 CAAG9G CAAHBT D46705`. There is no closed ST35 genome, so use
the Shovill assembly with the fewest contigs as reference (check `$RESULTS/02_assembly/quast_shovill`).

<details>
<summary><b>What changes when the reference is a draft assembly?</b></summary>

Snippy can map to a multi-contig reference. Regions at contig ends map poorly, and anything missing from the
draft cannot be compared, so the core genome is smaller. Distances are still valid for comparing isolates
within the group.

</details>

---

## End-of-day checklist

- [ ] `ST39.treefile` and `ST39_snp_dists_filtered.tsv` exist.
- [ ] You compared raw and recombination-filtered distances.
- [ ] You have min, median and max distances for within-cluster and cluster-versus-other pairs.
- [ ] Optional: ST35 tree and distances.

## Further reading

- [Snippy](https://github.com/tseemann/snippy), [Gubbins](https://github.com/nickjcroucher/gubbins),
  [IQ-TREE](https://iqtree.github.io)
- Cornick J *et al.* Genomic investigation of a suspected *Klebsiella pneumoniae* outbreak in a neonatal care
  unit in sub-Saharan Africa. *Microbial Genomics* 2021. [Full text on PubMed Central](https://pmc.ncbi.nlm.nih.gov/articles/PMC8743538)

---

[← Day 5](day05.md) · [Course home](README.md) · **Day 6** · [Day 7 →](day07.md)
