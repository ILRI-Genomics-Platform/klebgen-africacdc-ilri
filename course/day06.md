[← Day 5](day05.md) · [Course home](README.md) · **Day 6** · [Day 7 →](day07.md)

# Day 6: Building trees for transmission analysis

## At a glance

| Session | Time | Topic | Module |
| --- | --- | --- | --- |
| 1 | 09:00–10:30 | From reads to a core-genome SNP alignment (Snippy) | 6 |
| 2 | 11:00–12:30 | Recombination and why it matters (Gubbins) | 6 |
| 3 | 13:30–15:00 | Maximum-likelihood trees (IQ-TREE) | 6 |
| 4 | 15:30–17:00 | SNP distances, a reference-free check with SKA2, and what a threshold can and cannot tell you | 6 |

**You will produce:** a recombination-filtered tree and a SNP-distance matrix for the ST39 group.

## Learning objectives

By the end of today you can:

- Choose a reference genome and build a core-genome SNP alignment.
- Explain why recombination inflates SNP counts, and remove it with Gubbins.
- Build a maximum-likelihood tree with branch support.
- Calculate pairwise SNP distances and use a threshold as **evidence**, not proof, of transmission.
- Compare reference-based (Snippy) and reference-free (SKA2) SNP distances.

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

Snippy needs one reference genome for the whole group, from the same species. Choose it the way the pipeline
does: take each isolate's **species from Kleborate** (Day 3), then the reference listed for that species in the
pipeline's reference sheet, `conf/phylo_references.csv`. For *Klebsiella pneumoniae* that is **NTUH-K2044**
(accession AP006725.1): a complete, well-annotated chromosome (no plasmids, which are gained and lost too often
to follow the same history). Using the same reference as the pipeline also means you can compare your tree with
the pipeline's `Kpneumoniae_ST39` results.

```bash
use_tools phylo
mkdir -p $RESULTS/06_phylogeny/ST39 && cd $RESULTS/06_phylogeny/ST39
ST39="CAA65Z CAA99X CAA47J CAA4X5 CAA574 CAA5XX CAA5YM CAA4GD D29438 BKR8JH BHB7VE CHI11E"

# 1. Species of the ST39 isolates, from Kleborate: all should be Klebsiella pneumoniae
K=$RESULTS/03_typing/kleborate/klebsiella_pneumo_complex_output.txt
pick '(^|__)species$' $K | grep -wE "$(echo $ST39 | tr ' ' '|')"

# 2. The reference for that species (an ST39 row, if there were one, would be used instead)
grep -v '^#' $SHARED/neonatal-amr-nf/conf/phylo_references.csv | grep '^Klebsiella pneumoniae,'
#    Klebsiella pneumoniae,,NTUH-K2044,AP006725.1

# 3. The trainers downloaded it once into $SHARED/references, as FASTA and as GenBank (the same sequence,
#    plus its annotation: genes and what they code for)
cp $SHARED/references/AP006725.1.fasta ref_NTUH-K2044.fasta
cp $SHARED/references/AP006725.1.gbk   ref_NTUH-K2044.gbk
seqkit stats ref_NTUH-K2044.fasta          # one sequence, about 5.25 Mb
grep -c "^     CDS " ref_NTUH-K2044.gbk   # annotated genes
```

Snippy gets the **GenBank** file: with the annotation it reports, for every variant, the gene it falls in and its
effect (synonymous, missense, ...). The alignment steps (`snippy-core`, SKA2) use the FASTA, which is the same
sequence, so their coordinates match.

<details>
<summary><b>Would a reference from ST39 itself be better?</b></summary>

Usually, yes. NTUH-K2044 belongs to ST23. Genes that ST39 has and NTUH-K2044 lacks cannot be compared, so the core
genome is smaller than with an ST39 reference, such as the chromosome of your `CHI11E` hybrid assembly. Within one
ST the SNPs that separate isolates mostly fall in genes both share, so cluster conclusions rarely change, though
individual distances can. The pipeline takes an ST-specific reference from a row with that ST in
`conf/phylo_references.csv` (the sheet ends with an example); the trainers then download it with
`bin/fetch_references.py`.

</details>

### 6.3 Map reads and call SNPs with Snippy

```bash
# Your practice isolates (all three are ST39): 20-45 minutes each with 1 thread
for s in $PRACTICE; do
    snippy --outdir snippy/$s --ref ref_NTUH-K2044.gbk \
        --R1 $COURSE/trimmed/${s}_R1.fastq.gz --R2 $COURSE/trimmed/${s}_R2.fastq.gz \
        --cpus $THREADS --force
done
# The other nine ST39 isolates, from the trainers (same reference, same command)
use_precomputed 06_phylogeny/ST39/snippy
ls snippy | wc -l      # 12

# Variants with their gene and effect (columns: position, type, effect, gene, product)
cut -f2,3,11,13,14 snippy/CHI11E/snps.tab | head

# How many SNPs did each isolate have against the reference?
for s in $ST39; do echo -e "$s\t$(grep -c -v '^CHROM' snippy/$s/snps.tab)"; done
```

> [!NOTE]
> Every isolate has many thousands of SNPs against NTUH-K2044, because the reference belongs to another lineage
> (ST23). That is expected: what matters is how the ST39 isolates differ **from each other**, after the core
> alignment. All 12 should have similar counts; an isolate with far more or far fewer than the rest needs a
> closer look (coverage, contamination, wrong ST).

### 6.4 Core alignment

```bash
snippy-core --ref ref_NTUH-K2044.fasta --prefix core snippy/*
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

> [!TIP]
> **Fell behind, or a step failed?** The trainers ran 6.3–6.8 on all 12 isolates with the same reference and
> commands. `use_precomputed 06_phylogeny/ST39` copies every file you do not have yet (core alignment, Gubbins,
> IQ-TREE, SNP distances, SKA2) into this folder and keeps the ones you made.

### 6.8 A reference-free check: SKA2

Snippy compares isolates only where they map to the reference. **SKA2** (split k-mer analysis) needs no
reference: it cuts every genome into *split k-mers*, 31 bp pieces whose middle base is left out, and finds pairs
of isolates where the same 30 flanking bases have a different middle base. That is a SNP. It is fast, works from
reads or assemblies, and is well suited to closely related isolates such as one ST.

```bash
use_tools phylo
cd $RESULTS/06_phylogeny/ST39
ST39="CAA65Z CAA99X CAA47J CAA4X5 CAA574 CAA5XX CAA5YM CAA4GD D29438 BKR8JH BHB7VE CHI11E"   # as in 6.2

# 1. Split k-mers of every isolate, from its best assembly: one line per isolate, name <tab> file
for s in $ST39; do printf '%s\t%s\n' $s $COURSE/best/$s.fasta; done > ska_inputs.tsv
ska build -f ska_inputs.tsv -o ST39_ska -k 31 --threads $THREADS
#    (Short of time? use_precomputed 06_phylogeny/ST39/ST39_ska.skf 06_phylogeny/ST39/ST39_ska_distances.tsv)

# 2. SNP distances between every pair of isolates (columns: Sample1, Sample2, Distance, Mismatches)
ska distance ST39_ska.skf -o ST39_ska_distances.tsv
head ST39_ska_distances.tsv

# 3. Put the SKA2 distance next to the Snippy + Gubbins distance for each pair (from 6.7)
awk -F'\t' 'NR == FNR { if (FNR > 1) { k = ($1 < $2) ? $1" "$2 : $2" "$1; ska[k] = $3 }; next }
            { k = ($1 < $2) ? $1" "$2 : $2" "$1; print $1"\t"$2"\t"$3"\t"ska[k]"\t"$4 }' \
    ST39_ska_distances.tsv ST39_pairs.tsv > ST39_snippy_vs_ska.tsv
echo -e "isolate1\tisolate2\tsnippy_gubbins\tska2\tpair_type" | cat - ST39_snippy_vs_ska.tsv | expand -t 16 | less
```

<details>
<summary><b>Are the SKA2 distances larger or smaller than the Snippy + Gubbins ones? Why?</b></summary>

Several things pull in different directions:

- **Recombination**: SKA2's distances are raw; Gubbins removed recombined regions from the Snippy distances. If
  a pair differs by a recombination event, SKA2's distance is larger.
- **More genome**: SKA2 compares everything two isolates share, including genes and plasmids the reference
  (NTUH-K2044's chromosome) lacks, so it can find SNPs Snippy cannot see.
- **Fewer SNPs in dense regions**: a split k-mer needs 15 identical bases on each side, so SNPs close to another
  variant, and indels, are missed.

What matters for transmission is whether the **ranking** of pairs and the within-cluster distances agree: if both
methods put the cluster isolates a few SNPs apart and the background far away, the conclusion does not depend
on the method.

</details>

**The same comparison, from the pipeline.** The pipeline can build every group's phylogeny with both methods
(`--phylo_method both`): Snippy's results go to `06_phylogeny/<group>/`, SKA2's to `06_phylogeny/<group>_ska2/`.
For SKA2, `ska map` places the split k-mers on the same reference, so Gubbins, IQ-TREE and snp-dists run on its
alignment exactly as on Snippy's. Your trainers ran the training set this way; compare the two for every group:

```bash
P=$SHARED/precomputed/pipeline_results/06_phylogeny
ls $P                                          # Kpneumoniae_ST39/ and Kpneumoniae_ST39_ska2/, ...
python3 $SHARED/neonatal-amr-nf/bin/compare_snp_methods.py --phylo $P \
    --threshold 20 --out $RESULTS/06_phylogeny/snippy_vs_ska2
expand -t 14 $RESULTS/06_phylogeny/snippy_vs_ska2_pairs.tsv | less -S
```

For each group, the summary gives the **Spearman correlation** of the two sets of distances (1 = both methods
rank the pairs the same way), the typical and largest **difference**, and how many pairs each method calls
**close** (20 SNPs or fewer). The pairs table has one row per pair of isolates, with both distances after
Gubbins and SKA2's reference-free distance (`ska2_ref_free`).

<details>
<summary><b>A pair is close by one method but not by the other. What do you do?</b></summary>

Treat the pair as uncertain, and find out why before deciding. Look at the pair in both trees, and at Gubbins'
recombination predictions for each (`06_phylogeny/gubbins/<group>.recombination_predictions.gff` and
`<group>_ska2.recombination_predictions.gff`): a recombination block detected by one alignment but not the other
moves many SNPs at once. Check the isolates' QC too: low coverage or a mixed sample affects the two methods
differently. A conclusion about transmission should hold with either method; if it does not, report that it
depends on the method.

</details>

### 6.9 What does a SNP threshold mean?

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
- [ ] You compared Snippy + Gubbins distances with SKA2's (`ST39_snippy_vs_ska.tsv`), and the pipeline's two
      methods for every group (`snippy_vs_ska2_summary.tsv`).
- [ ] Optional: ST35 tree and distances.

## Further reading

- [Snippy](https://github.com/tseemann/snippy), [Gubbins](https://github.com/nickjcroucher/gubbins),
  [IQ-TREE](https://iqtree.github.io), [SKA2](https://github.com/bacpop/ska.rust)
- Derelle R *et al.* Seamless, rapid, and accurate analyses of outbreak genomic data using split k-mer analysis.
  *Genome Research* 2024;34:1661–1673. [doi:10.1101/gr.279449.124](https://doi.org/10.1101/gr.279449.124)
- Cornick J *et al.* Genomic investigation of a suspected *Klebsiella pneumoniae* outbreak in a neonatal care
  unit in sub-Saharan Africa. *Microbial Genomics* 2021. [Full text on PubMed Central](https://pmc.ncbi.nlm.nih.gov/articles/PMC8743538)

---

[← Day 5](day05.md) · [Course home](README.md) · **Day 6** · [Day 7 →](day07.md)
