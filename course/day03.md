[← Day 2](day02.md) · [Course home](README.md) · **Day 3** · [Day 4 →](day04.md)

# Day 3: Who is causing infection, and what resistance do they carry?

## At a glance

| Session | Time | Topic | Module |
| --- | --- | --- | --- |
| 1 | 09:00–10:30 | Species identification within the *K. pneumoniae* species complex | 3 |
| 2 | 11:00–12:30 | MLST, high-risk clones, global context in Pathogenwatch | 3 |
| 3 | 13:30–15:00 | AMR gene detection with AMRFinderPlus and Kleborate | 4 |
| 4 | 15:30–17:00 | Virulence loci and convergent strains | 4 |

**You will produce:** a one-page lineage summary, and a gene-by-isolate resistance table.

## Learning objectives

By the end of today you can:

- Confirm species within the *K. pneumoniae* species complex (KpSC), and explain why lab methods often miss
  *K. quasipneumoniae* and *K. variicola*.
- Assign MLST sequence types (STs) and recognise high-risk clones.
- Detect acquired resistance genes and chromosomal mutations, and tell intrinsic genes from acquired ones.
- Read Kleborate's resistance and virulence scores, and spot convergent (resistant and hypervirulent) strains.

> [!IMPORTANT]
> Today uses the **best assembly** for each isolate, in `$COURSE/best/` (made at the end of Day 2).

---

## Module 3: Species identification and typing

### 3.1 Kleborate

Kleborate screens KpSC assemblies for species, MLST, resistance and virulence loci, and capsule (K) and O-antigen
types in a single run.

```bash
use_tools typing
kleborate --version
mkdir -p $RESULTS/03_typing

kleborate -a $COURSE/best/*.fasta -o $RESULTS/03_typing/kleborate -p kpsc
ls $RESULTS/03_typing/kleborate
```

Kleborate writes one wide table. The `pick` function from your Day 1 setup prints just the columns whose names
match a pattern:

```bash
K=$RESULTS/03_typing/kleborate/klebsiella_pneumo_complex_output.txt
head -1 $K | tr '\t' '\n' | head -30      # list the column names
pick 'species$|__ST$' $K
```

### 3.2 Exercise: species

Compare Kleborate's species call with the label each isolate arrived with. In this study every isolate was
labelled *K. pneumoniae* by the lab.

<details>
<summary><b>Which isolates are not <i>K. pneumoniae</i>?</b></summary>

From the study's published typing (`metadata_compact.tsv`, column `published_species`):

| Isolate | Species |
| --- | --- |
| `1013757` | *K. variicola* subsp. *variicola* (one of the three hybrid isolates) |
| `CAAYXM` | *K. variicola* subsp. *variicola* |
| `CAAQN6` | *K. quasipneumoniae* subsp. *similipneumoniae* |

All three belong to the KpSC, which is why Kleborate types them too. MALDI-TOF and biochemical tests often call
them *K. pneumoniae*. The difference matters for surveillance counts and, for *K. variicola*, for some intrinsic
resistance genes (LEN instead of SHV beta-lactamases).

</details>

### 3.3 Exercise: lineages

```bash
pick '__ST$' $K | tail -n +2 | cut -f2 | sort | uniq -c | sort -rn
```

Now compare with the published STs:

```bash
cut -f1,7 $COURSE/training_data/metadata_compact.tsv
```

<details>
<summary><b>Which ST is most common, and is that a fair estimate of what causes neonatal infection in Blantyre?</b></summary>

**ST39** (12 of 29), then **ST35** (6). But this is **not** a fair estimate: the training set was chosen to
contain two clusters, so those STs are over-represented on purpose. Prevalence must come from a representative
sample (all isolates, or a random subset), never from a hand-picked one. This is the most common error in
genomic surveillance reports.

</details>

**High-risk clones to recognise:** *K. pneumoniae* ST307, ST147, ST15, ST17, ST39, ST14; *E. coli* ST131, ST410,
ST167, ST648.

### 3.4 Global context in Pathogenwatch

1. Go to [pathogen.watch](https://pathogen.watch) and create a free account.
2. Upload 5 assemblies from `$COURSE/best/` (drag and drop the FASTA files).
3. Compare Pathogenwatch's ST and species with yours.
4. Open the **LIN code** and **cgMLST** results: which of your isolates fall closest to genomes from other countries?

> [!WARNING]
> Upload only data you are allowed to share. The training data here is public; your own surveillance data may not
> be until your data-sharing agreements allow it.

### 3.5 *E. coli* typing (optional)

Only if your trainers have added *E. coli* genomes:

```bash
kleborate -a $COURSE/best/EC_*.fasta -o $RESULTS/03_typing/kleborate_ecoli -p escherichia   # MLST
ectyper -i $COURSE/best/EC_01.fasta -o $RESULTS/03_typing/ectyper_EC_01                    # O:H serotype
ezclermont $COURSE/best/EC_01.fasta                                                         # phylogroup
```

### Check your understanding

<details>
<summary><b>1. Kleborate reports an ST with a suffix such as <code>ST14-1LV</code>. What does that mean?</b></summary>

The allele profile matches ST14 at all but one of the seven MLST genes: a **single-locus variant**. It is closely
related to ST14 but is not an exact match to a defined ST, perhaps because the allele is new or the gene sits at a
contig break.

</details>

<details>
<summary><b>2. Why might two isolates with the same ST still be unrelated in a hospital outbreak?</b></summary>

MLST looks at seven genes. Isolates of one ST can differ by hundreds of SNPs across the genome and be years
apart. Transmission needs whole-genome resolution (Day 6).

</details>

---

## Module 4: Resistance and virulence (part 1: detection)

### 4.1 AMRFinderPlus

```bash
use_tools amr
amrfinder --version
AMRDB=$SHARED/databases/amrfinderplus/latest
amrfinder -d $AMRDB --database_version
mkdir -p $RESULTS/04_amr/amrfinderplus

for f in $COURSE/best/*.fasta; do
    s=$(basename $f .fasta)
    amrfinder -n $f --organism Klebsiella_pneumoniae --plus --name $s \
        -d $AMRDB --threads $THREADS -o $RESULTS/04_amr/amrfinderplus/${s}.tsv
done

# All results in one table (one header line)
awk 'FNR==1 && NR!=1 {next} {print}' $RESULTS/04_amr/amrfinderplus/*.tsv > $RESULTS/04_amr/amrfinderplus_all.tsv
head -3 $RESULTS/04_amr/amrfinderplus_all.tsv
```

> [!NOTE]
> `--organism Klebsiella_pneumoniae` turns on checks for species-specific point mutations (for example in *gyrA*,
> *parC*, *ompK35/36*). `--plus` adds virulence, stress and metal-resistance genes.

### 4.2 Exercise: the key genes

For each group below, count how many of the 29 isolates carry at least one:

| Group | Examples | Isolates with ≥ 1 |
| --- | --- | --- |
| ESBL | *bla*CTX-M-15, SHV ESBL variants | |
| Carbapenemase | *bla*NDM, *bla*OXA-48-like, *bla*KPC, *bla*VIM | |
| Aminoglycoside-modifying enzymes | *aac(3)-II*, *aac(6')-Ib*, *aph(3')-VI* | |
| 16S rRNA methylases | *armA*, *rmtB* | |
| Colistin | *mcr*, *mgrB* disruption | |
| Fluoroquinolone | *qnr*, *gyrA*/*parC* mutations | |

```bash
A=$RESULTS/04_amr/amrfinderplus_all.tsv
# Sample, gene and drug class for every hit
pick 'Element symbol|^Class$|^Subclass$' $A | head
# Example: carbapenemases
pick 'Element symbol|^Subclass$' $A | grep -i carbapenem | sort -u
```

<details>
<summary><b>Answers from the published Kleborate results</b></summary>

- **CTX-M-15**: 21 of 29 isolates.
- **NDM-1**: 2 isolates, `BHB7VE` (ST39) and `CES1C7` (ST45).
- **Acquired aminoglycoside resistance genes**: 27 of 29; 23 carry an *aac(3)* gene, which typically confers
  gentamicin resistance.

Your counts may differ slightly: AMRFinderPlus and Kleborate use different databases and versions.

</details>

### 4.3 Exercise: AMRFinderPlus versus Kleborate

```bash
pick 'Bla_ESBL_acquired$|Bla_Carb_acquired$|AGly_acquired$|resistance_score$' $K
```

Choose 5 isolates and compare the two tools gene by gene. Write down every difference and suggest a reason
(different database, gene name, truncated gene, contig break).

### 4.4 Intrinsic versus acquired

*K. pneumoniae* always carries a chromosomal *bla*SHV (usually SHV-1 or SHV-11), and *K. variicola* a *bla*LEN.
These make the species naturally resistant to ampicillin. Kleborate lists them in a separate column (`Bla_chr`)
so they are not mistaken for acquired resistance.

<details>
<summary><b>Why does it matter that ampicillin resistance is intrinsic?</b></summary>

The WHO-recommended first-line empiric regimen for neonatal sepsis is ampicillin plus gentamicin. Against
*K. pneumoniae*, ampicillin adds nothing: the regimen effectively relies on **gentamicin alone**. So gentamicin
resistance (for example *aac(3)-II*) is what decides whether first-line treatment can work.

</details>

### 4.5 Virulence and convergence

```bash
pick 'virulence_score$|Yersiniabactin$|Aerobactin$|Salmochelin$|Colibactin$|RmpADC$|rmpA2$' $K
```

| Locus | What it does |
| --- | --- |
| Yersiniabactin (*ybt*) | Siderophore; common in hospital lineages |
| Aerobactin (*iuc*) | Siderophore strongly linked to hypervirulence |
| Salmochelin (*iro*) | Siderophore, often with *iuc* on the virulence plasmid |
| Colibactin (*clb*) | Genotoxin |
| *rmpADC* / *rmpA2* | Hypermucoidy |

> [!WARNING]
> A **convergent** strain carries both a carbapenemase and aerobactin or *rmp*. Report it immediately: it combines
> very limited treatment options with high virulence.

<details>
<summary><b>Which isolates have the highest virulence scores, and is either convergent?</b></summary>

From the published Kleborate (v2) results:

| Isolate | Score | Loci |
| --- | --- | --- |
| `CAAUM3` (ST1049, KL5) | 4 | yersiniabactin (ybt 2, ICEKp1); aerobactin iuc 3 **truncated**; salmochelin iro 3 **truncated**; *rmp* 3 (ICEKp1) |
| `1013757` (*K. variicola*) | 3 | aerobactin iuc 5; salmochelin iro 5 |

Neither carries a carbapenemase, so neither is convergent. Two points to discuss:

- In `CAAUM3` the aerobactin genes are **truncated**. Does a broken locus still count towards hypervirulence?
- `1013757` is a hybrid isolate with closed plasmids. On Day 4 you can find which plasmid carries *iuc* 5.

</details>

---

## End-of-day checklist

- [ ] Kleborate results for all 29 isolates in `$RESULTS/03_typing/kleborate/`.
- [ ] You found the three non-*K. pneumoniae* isolates.
- [ ] One-page lineage summary: an ST count table (or bar chart) and three sentences of interpretation, including
      why these counts are not prevalence.
- [ ] AMRFinderPlus results combined in `$RESULTS/04_amr/amrfinderplus_all.tsv`.
- [ ] The key-genes table (4.2) is filled in.

## Further reading

- Lam MMC *et al.* A genomic surveillance framework and genotyping tool for *Klebsiella pneumoniae* and its related
  species complex. *Nat Commun* 2021. [Kleborate docs](https://github.com/klebgenomics/Kleborate)
- [AMRFinderPlus wiki](https://github.com/ncbi/amr/wiki)

---

[← Day 2](day02.md) · [Course home](README.md) · **Day 3** · [Day 4 →](day04.md)
