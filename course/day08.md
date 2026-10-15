[← Day 7](day07.md) · [Course home](README.md) · **Day 8** · [Day 9 →](day09.md)

# Day 8: Vaccine targets, and bringing all the results together

## At a glance

| Session | Time | Topic | Module |
| --- | --- | --- | --- |
| 1 | 09:00–10:30 | Capsule (K) and O-antigen typing with Kaptive | 7 |
| 2 | 11:00–12:30 | What would a vaccine cover? | 7 |
| 3 | 13:30–15:00 | One table for everything: joining results in R | 8 |
| 4 | 15:30–17:00 | Checking the joined table; summary statistics | 8 |

**You will produce:** a vaccine-coverage figure, and a master table with one row per isolate.

## Learning objectives

By the end of today you can:

- Explain K-locus (KL) and O-locus (OL) nomenclature, and the difference between a locus and a predicted type.
- Run Kaptive 3 and judge its confidence calls.
- Calculate the theoretical coverage of a multivalent vaccine, and explain how sampling and time affect it.
- Join typing, AMR, plasmid and metadata tables into one tidy table in R.

---

## Module 7: Antigen typing and vaccine relevance

### 7.1 Background

*K. pneumoniae* has two main surface polysaccharides, both candidates for maternal vaccines against neonatal sepsis:

| | Capsule (K antigen) | O antigen (lipopolysaccharide) |
| --- | --- | --- |
| Diversity | High: well over 100 known K loci | Low: about a dozen types |
| Locus name | KL1, KL2, ... | O1/O2v1, O3b, ... |
| Predicted type | K1, K2, ... or "unknown" | O1, O2, O3b, O4, O5, ... |
| Vaccine implication | Needs many valencies | A few valencies can cover most strains |

### 7.2 Kaptive 3

```bash
use_tools typing
mkdir -p $RESULTS/07_antigens

# Your practice isolates
practice=$(for s in $PRACTICE; do echo $COURSE/best/$s.fasta; done)
kaptive assembly kpsc_k $practice -o $RESULTS/07_antigens/kaptive_k_practice.tsv -t $THREADS
kaptive assembly kpsc_o $practice -o $RESULTS/07_antigens/kaptive_o_practice.tsv -t $THREADS
pick 'Best match locus|Best match type|Match confidence|Problems' $RESULTS/07_antigens/kaptive_k_practice.tsv

# The trainers' tables for all 29 isolates (kaptive_k.tsv, kaptive_o.tsv), used below
use_precomputed 07_antigens
pick 'Best match locus|Best match type|Match confidence|Problems' $RESULTS/07_antigens/kaptive_k.tsv
```

> [!NOTE]
> Kaptive 3 searches gene by gene, which makes it much better than earlier versions at typing loci split across
> several contigs. The K locus is hard to assemble from short reads because its GC content is lower than the rest
> of the genome, so always look at the **confidence** column.

### 7.3 Exercise: check your calls

1. Compare your K and O calls with the published ones (`published_K_locus`, `published_O_type` in
   `metadata_compact.tsv`).
2. For the three hybrid isolates, run Kaptive on both the Shovill and the hybrid assembly. Does the call or the
   confidence change?

```bash
for s in CHI11E CNS1F1 1013757; do
    kaptive assembly kpsc_k $COURSE/assemblies/${s}.fasta $COURSE/assemblies/${s}.hybrid.fasta -t $THREADS \
        | cut -f1-4
done
```

<details>
<summary><b>What K locus do the eight ST39 cluster isolates carry, and why might that matter?</b></summary>

**KL2** in the published typing (O type O1b), the same K locus as many ST14 isolates in this hospital. A clone can
also **switch capsule** by recombination at the K locus: the ST39 isolates in this set carry KL2, KL23, KL51,
KL62 and KL149. A capsule-based vaccine targets the capsule, not the ST, so capsule switching could let a lineage
escape it.

</details>

### 7.4 What would a vaccine cover?

Coverage from 29 hand-picked isolates would be meaningless. Use the published K and O types for **all 1,373
isolates** in the study instead (`training_data/published_K_O_all.tsv`).

```bash
use_tools r
cd $COURSE && R
```

```r
library(tidyverse)
ko <- read_tsv("training_data/published_K_O_all.tsv")

# Group O types. Careful: "O12" starts with "O1" but is a different O type.
ko <- ko |> mutate(O_group = case_when(
  str_starts(O_type, "unknown")     ~ "unknown",
  O_type %in% c("O1a", "O1b", "O1") ~ "O1",
  str_starts(O_type, "O2")          ~ "O2",
  TRUE                              ~ O_type
))

neo <- ko |> filter(neonate == "yes")

# O-antigen distribution in neonatal isolates
neo |> count(O_group, sort = TRUE) |> mutate(percent = round(100 * n / sum(n), 1))

# Coverage of a four-valent O-antigen vaccine (O1, O2, O3b, O5)
valencies <- c("O1", "O2", "O3b", "O5")
neo |> summarise(n = n(), covered = sum(O_group %in% valencies), percent = round(100 * covered / n, 1))

# How many K loci are needed to cover 50% and 70% of neonatal isolates?
k_curve <- neo |>
  filter(!str_detect(K_locus, "unknown")) |>
  count(K_locus, sort = TRUE) |>
  mutate(rank = row_number(), cumulative_percent = 100 * cumsum(n) / nrow(neo))
k_curve |> filter(cumulative_percent >= 50) |> slice(1)
k_curve |> filter(cumulative_percent >= 70) |> slice(1)

ggplot(k_curve, aes(rank, cumulative_percent)) +
  geom_step() + geom_hline(yintercept = c(50, 70), linetype = 2) +
  labs(x = "Number of K loci in vaccine (most common first)",
       y = "Neonatal isolates covered (%)",
       title = "Capsule-based vaccine: valencies needed") +
  theme_minimal()
ggsave("results/07_antigens/k_coverage_curve.png", width = 7, height = 4)

# Does coverage change over time?
neo |>
  filter(!is.na(collection_year)) |>
  mutate(period = cut(collection_year, c(1995, 2009, 2015, 2020),
                      labels = c("1996-2009", "2010-2015", "2016-2020"))) |>
  group_by(period) |>
  summarise(n = n(), O4valent_percent = round(100 * mean(O_group %in% valencies), 1),
            KL2_percent = round(100 * mean(K_locus == "KL2"), 1))
```

<details>
<summary><b>Answers (from the published data)</b></summary>

Neonatal isolates: **388**.

| Measure | Result |
| --- | --- |
| O1 | 75.3% |
| O2 | 10.8% |
| Four-valent O1/O2/O3b/O5 coverage | **89.9%** (83.6% across all 1,373 isolates) |
| Distinct K loci (excluding unknown) | 43 |
| K loci needed for ≥ 50% coverage | **3** (KL2 alone is 39.7%) |
| K loci needed for ≥ 70% coverage | **8** |

By period (neonatal isolates):

| Period | n | Four-valent O coverage | KL2 |
| --- | --- | --- | --- |
| 1996–2009 | 35 | 94.3% | 22.9% |
| 2010–2015 | 30 | 73.3% | 26.7% |
| 2016–2020 | 321 | 91.0% | 43.0% |

Isolates with an unknown O type count as **not covered**. If your O1 percentage came out slightly higher, check
whether you grouped O12 into O1.

</details>

<details>
<summary><b>Discussion: how much should you trust these coverage figures?</b></summary>

- **One hospital.** Antigen distributions differ between sites and countries.
- **Time.** Most neonatal isolates are from 2016–2020, and the earlier periods have only about 30 isolates each.
  The study found repeated replacement of dominant clones over 20 years, so coverage can shift as clones change.
- **KL2 is driven partly by clusters.** If one KL2 clone caused an outbreak, KL2 looks more common than it is
  across independent infections. One option is to count each cluster once (deduplicate by cluster).

</details>

---

## Module 8: Integration, visualisation and reporting (part 1)

### 8.1 One row per isolate

```r
library(tidyverse)
setwd(Sys.getenv("COURSE"))

meta <- read_tsv("training_data/metadata_compact.tsv")

# Kleborate: look at the column names first
kleb_raw <- read_tsv("results/03_typing/kleborate/klebsiella_pneumo_complex_output.txt")
names(kleb_raw)
kleb <- kleb_raw |>
  rename(sample = 1) |>
  select(sample,
         species          = matches("__species$"),
         ST               = matches("mlst__ST$"),
         K_locus          = matches("K_locus$"),
         O_type           = matches("O_type$"),
         resistance_score = matches("resistance_score$"),
         virulence_score  = matches("virulence_score$"))

# AMRFinderPlus: key resistance flags per isolate
amr <- read_tsv("results/04_amr/amrfinderplus_all.tsv") |>
  rename(sample = Name, gene = any_of(c("Element symbol", "Gene symbol")))
amr_flags <- amr |>
  group_by(sample) |>
  summarise(
    ESBL_or_AmpC   = any(str_detect(Subclass, "CEPHALOSPORIN"), na.rm = TRUE),
    carbapenemase  = any(str_detect(Subclass, "CARBAPENEM"), na.rm = TRUE),
    gentamicin_R   = any(str_detect(Subclass, "GENTAMICIN"), na.rm = TRUE),
    amr_genes      = paste(sort(unique(gene[Type == "AMR"])), collapse = ";")
  )

# Plasmids (MOB-suite): number of plasmids per isolate
plasmids <- read_tsv("results/05_plasmids/all_plasmids.tsv", col_names = FALSE) |>
  count(X1, name = "n_plasmids") |> rename(sample = X1)

master <- meta |>
  left_join(kleb, by = "sample") |>
  left_join(amr_flags, by = "sample") |>
  left_join(plasmids, by = "sample") |>
  mutate(n_plasmids = replace_na(n_plasmids, 0))

dir.create("results/08_integration", showWarnings = FALSE, recursive = TRUE)
write_tsv(master, "results/08_integration/master_table.tsv")
glimpse(master)
```

> [!WARNING]
> After every join, check the number of rows: `nrow(master)` should still be 29. More rows means a sample ID
> appears twice in one table; fewer rows (after an `inner_join`) means IDs do not match. A common cause is a
> different spelling of the sample name, for example a `.fasta` suffix left on in one table.

### 8.2 Check your table against the published results

```r
master |>
  mutate(ST_agrees = ST == published_ST,
         K_agrees  = K_locus == str_remove(published_K_locus, " .*")) |>
  count(ST_agrees, K_agrees)

master |> filter(ST != published_ST) |> select(sample, ST, published_ST)
```

Explain each disagreement: a new allele, a single-locus variant (`-1LV`), a contig break, or a tool version
difference?

### 8.3 Summary statistics

```r
master |>
  summarise(isolates = n(),
            neonatal = sum(neonate == "yes"),
            ESBL = sum(ESBL_or_AmpC, na.rm = TRUE),
            carbapenemase = sum(carbapenemase, na.rm = TRUE),
            gentamicin_R = sum(gentamicin_R, na.rm = TRUE))

master |> count(ST, sort = TRUE)
```

---

## End-of-day checklist

- [ ] Kaptive K and O results for all 29 isolates, compared with the published calls.
- [ ] Vaccine-coverage results and the K coverage curve (`results/07_antigens/k_coverage_curve.png`).
- [ ] `results/08_integration/master_table.tsv` with 29 rows.
- [ ] Every disagreement with the published typing explained.

## Further reading

- Stanton TD *et al.* Fast and accurate in silico antigen typing with Kaptive 3. *bioRxiv* 2025.
  [Kaptive docs](https://kaptive.readthedocs.io)
- Heinz E *et al.* Longitudinal analysis within one hospital in sub-Saharan Africa over 20 years reveals repeated
  replacements of dominant clones of *Klebsiella pneumoniae*. *Genome Med* 2024.
  [doi:10.1186/s13073-024-01342-3](https://doi.org/10.1186/s13073-024-01342-3)

---

[← Day 7](day07.md) · [Course home](README.md) · **Day 8** · [Day 9 →](day09.md)
