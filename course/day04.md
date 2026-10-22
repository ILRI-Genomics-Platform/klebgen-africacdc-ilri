[← Day 3](day03.md) · [Course home](README.md) · **Day 4** · [Day 5 →](day05.md)

# Day 4: Will treatment work? Where do the resistance genes live?

## At a glance

| Session | Time | Topic | Module |
| --- | --- | --- | --- |
| 1 | 09:00–10:30 | Genotype versus phenotype: how far can we trust predictions? | 4 |
| 2 | 11:00–12:30 | From genes to empiric therapy for neonatal sepsis | 4 |
| 3 | 13:30–15:00 | Plasmids, replicons and mobility; MOB-suite | 5 |
| 4 | 15:30–17:00 | Which plasmid carries which gene? Short-read versus hybrid | 5 |

**You will produce:** a genotype–phenotype agreement table for AMRFinderPlus and Kleborate, a regimen summary, and a plasmid table for the three
hybrid isolates.

## Learning objectives

By the end of today you can:

- Predict susceptibility from AMR genes (AMRFinderPlus and Kleborate) and measure agreement with lab AST (very major
  and major errors).
- Explain disagreements between genotype and phenotype.
- Estimate how often WHO-recommended empiric neonatal regimens would fail, and say why that estimate is limited.
- Reconstruct plasmids with MOB-suite and say which plasmid carries each resistance gene.
- Explain why short-read assemblies rarely give complete plasmids.

---

## Module 4: Resistance and virulence (part 2: interpretation)

### 4.6 Genotype versus phenotype in R

The study tested isolates by disc diffusion. Your AST table is `training_data/ast_compact.tsv` (S, I, R or NA).

```bash
use_tools r
cd $COURSE
R
```

Two tools predicted resistance genes on Day 3: AMRFinderPlus (`results/04_amr/`) and Kleborate
(`results/03_typing/kleborate/`). Turn each tool's genes into an S/R prediction per drug, then compare both with
the lab's AST in the same way.

```r
library(tidyverse)

meta <- read_tsv("training_data/metadata_compact.tsv")
ast  <- read_tsv("training_data/ast_compact.tsv", na = "NA")

ast_long <- ast |>
  pivot_longer(-sample, names_to = "ast_column", values_to = "phenotype") |>
  filter(phenotype %in% c("S", "R"))          # intermediate (I) results left out

# Genotype (S/R per sample and drug) against the AST, for either tool
compare_ast <- function(pred) {
  pred |>
    inner_join(ast_long, by = c("sample", "ast_column")) |>
    mutate(result = case_when(
      genotype == phenotype                  ~ "agree",
      genotype == "S" & phenotype == "R"     ~ "very major error (missed resistance)",
      genotype == "R" & phenotype == "S"     ~ "major error (false resistance)"
    ))
}

# ---- A. AMRFinderPlus: an isolate is predicted R if any hit's Subclass matches the pattern
amr <- read_tsv("results/04_amr/amrfinderplus_all.tsv") |>
  rename(sample = Name, gene = any_of(c("Element symbol", "Gene symbol"))) |>
  select(sample, gene, Class, Subclass)

rules <- tribble(
  ~drug,           ~ast_column,                    ~pattern,
  "Gentamicin",    "Gentamicin 10",                "GENTAMICIN",
  "Ceftriaxone",   "Cefpodoxime 10 (Ceftriaxone)", "CEPHALOSPORIN|CARBAPENEM",
  "Meropenem",     "Meropenem 10",                 "CARBAPENEM",
  "Amikacin",      "Amikacin 30",                  "AMIKACIN",
  "Ciprofloxacin", "Ciprofloxacin 5",              "QUINOLONE"
)

pred <- rules |>
  crossing(sample = meta$sample) |>
  left_join(amr, by = "sample", relationship = "many-to-many") |>
  group_by(drug, ast_column, sample) |>
  summarise(
    hit      = any(str_detect(Subclass, pattern), na.rm = TRUE),
    genes    = paste(unique(gene[which(str_detect(Subclass, pattern))]), collapse = ";"),
    .groups  = "drop"
  ) |>
  mutate(genotype = if_else(hit, "R", "S"))

# ---- B. Kleborate: one column per drug class (AGly_acquired, Bla_ESBL_acquired, Flq_mutations, ...), each
# holding that isolate's genes or mutations separated by ';' ('-' = none). An isolate is predicted R if a listed
# column holds a gene matching the pattern ('.' = any entry).
kleb_table <- read_tsv("results/03_typing/kleborate/klebsiella_pneumo_complex_output.txt",
                       col_types = cols(.default = "c")) |>
  rename(sample = 1)                                    # first column: the assembly (sample) name

kleb <- kleb_table |>
  select(sample, matches("_acquired$|_mutations$")) |>  # acquired genes and resistance mutations
  rename_with(\(x) sub(".*__", "", x), -sample) |>      # klebsiella_pneumo_complex__amr__AGly_acquired -> AGly_acquired
  pivot_longer(-sample, names_to = "column", values_to = "gene") |>
  separate_longer_delim(gene, ";") |>
  filter(!is.na(gene), gene != "-")

kleb_rules <- tribble(
  ~drug,           ~ast_column,                    ~column,                  ~pattern,
  "Gentamicin",    "Gentamicin 10",                "AGly_acquired",          "^(aac\\(3\\)|ant\\(2''\\)|armA|rmt|npmA)",
  "Ceftriaxone",   "Cefpodoxime 10 (Ceftriaxone)", "Bla_ESBL_acquired",      ".",
  "Ceftriaxone",   "Cefpodoxime 10 (Ceftriaxone)", "Bla_ESBL_inhR_acquired", ".",
  "Ceftriaxone",   "Cefpodoxime 10 (Ceftriaxone)", "Bla_Carb_acquired",      ".",
  "Ceftriaxone",   "Cefpodoxime 10 (Ceftriaxone)", "Bla_acquired",           "^(CMY|DHA|ACC|FOX|MOX)",   # plasmid AmpC
  "Meropenem",     "Meropenem 10",                 "Bla_Carb_acquired",      ".",
  "Amikacin",      "Amikacin 30",                  "AGly_acquired",          "^(aac\\(6'\\)-Ib|aph\\(3'\\)-VI|armA|rmt|npmA)",
  "Ciprofloxacin", "Ciprofloxacin 5",              "Flq_acquired",           ".",
  "Ciprofloxacin", "Ciprofloxacin 5",              "Flq_mutations",          "."
)

pred_kleb <- kleb_rules |>
  crossing(sample = kleb_table$sample) |>              # only isolates Kleborate typed (KpSC)
  left_join(kleb, by = c("sample", "column"), relationship = "many-to-many") |>
  group_by(drug, ast_column, sample) |>
  summarise(
    hit      = any(str_detect(gene, pattern), na.rm = TRUE),
    genes    = paste(unique(gene[which(str_detect(gene, pattern))]), collapse = ";"),
    .groups  = "drop"
  ) |>
  mutate(genotype = if_else(hit, "R", "S"))

# ---- Both tools against the AST
cmp <- bind_rows(AMRFinderPlus = compare_ast(pred), Kleborate = compare_ast(pred_kleb), .id = "tool")

cmp |>
  count(drug, tool, result) |>
  pivot_wider(names_from = result, values_from = n, values_fill = 0) |>
  arrange(drug, tool)
write_tsv(cmp, "results/04_amr/genotype_vs_phenotype.tsv")
```

The table has one row per drug and tool. Read the two rows of each drug together: where they differ, the tools
found different genes, or one reports a gene the other files under a different drug class.

> [!NOTE]
> Kleborate types only the *Klebsiella pneumoniae* species complex, so E. coli isolates (none in the training
> set, but common in your own data) get AMRFinderPlus predictions only. Kleborate's `Bla_chr` column (the intrinsic
> SHV, LEN or OKP of each species) is left out on purpose, as are the `truncated_resistance_hits` and
> `spurious_resistance_hits` columns: Kleborate puts partial or low-identity hits there instead of counting them.

> [!NOTE]
> **Very major error**: the genome says susceptible but the isolate is resistant. This is the dangerous one,
> because a patient could get an antibiotic that does not work. **Major error**: the genome says resistant but the
> isolate is susceptible, so a usable drug might be withheld.

### 4.7 Exercise: explain the disagreements

Pick two disagreements and investigate them:

```r
cmp |> filter(result != "agree") |> select(tool, sample, drug, genes, genotype, phenotype) |> arrange(drug, sample)

# Isolates where the two tools disagree with each other
inner_join(pred      |> select(sample, drug, AMRFinderPlus = genotype, amrfinder_genes = genes),
           pred_kleb |> select(sample, drug, Kleborate = genotype, kleborate_genes = genes),
           by = c("sample", "drug")) |>
  filter(AMRFinderPlus != Kleborate)
```

Possible explanations to consider:

- A gene that confers only **low-level** resistance, below the clinical breakpoint.
- Resistance from **porin loss** or **efflux**, which gene screening misses.
- A **truncated** or partial gene.
- An AST error: wrong breakpoints, wrong disc, mixed culture.
- A **database** or tool version difference: compare the two tools' rows for the same isolate and drug.

<details>
<summary><b>Hint: amikacin</b></summary>

Look at which genes drive the amikacin predictions. Genes such as *aac(6')-Ib* (and its *-cr* variant, common in
CTX-M-15 lineages) modify amikacin, but often not enough to cross the breakpoint. Expect some **major errors**
here. True high-level amikacin resistance usually comes from 16S rRNA methylases (*armA*, *rmtB*).

`BHB7VE` (ST39, NDM-1) is **amikacin-resistant** by AST. Which genes does it carry that could explain this? In the
published results it has *aph(3')-VI* and *aac(6')-Ib'*; check what your AMRFinderPlus run found.

</details>

<details>
<summary><b>What is missing from the AST for the two NDM-1 isolates?</b></summary>

`BHB7VE` has no meropenem result, and `CES1C7` has no AST at all. The genome is the only evidence of carbapenem
resistance for these isolates. In practice, a carbapenemase gene in a neonatal isolate should trigger phenotypic
confirmation and an infection-control alert, whatever the AST shows.

</details>

### 4.8 From genes to empiric therapy

```r
regimen <- pred |>
  select(sample, drug, genotype) |>
  pivot_wider(names_from = drug, values_from = genotype) |>
  transmute(
    sample,
    `Ampicillin + gentamicin` = Gentamicin,       # ampicillin is intrinsic R in KpSC
    `Third-generation cephalosporin` = Ceftriaxone,
    Meropenem,
    Amikacin
  )

regimen |>
  pivot_longer(-sample, names_to = "regimen", values_to = "predicted") |>
  group_by(regimen) |>
  summarise(n = n(), predicted_failure = sum(predicted == "R"),
            percent = round(100 * predicted_failure / n))
```

<details>
<summary><b>Why can't you report these percentages as "the regimen fails in X% of neonatal infections in Malawi"?</b></summary>

Three reasons:

1. The 29 isolates were **chosen**, not sampled at random: two clusters and two carbapenemase isolates are
   deliberately over-represented.
2. Only 15 of the 29 are from neonates.
3. Isolates come from one hospital.

A defensible estimate needs all (or a random sample of) neonatal bloodstream isolates from defined facilities
and periods, with the denominator stated.

</details>

---

## Module 5: Plasmids and mobile genetic elements (part 1)

### 5.1 Background

| Term | Meaning |
| --- | --- |
| **Replicon (Inc) type** | The plasmid's replication machinery, used to classify it (IncFII, IncFIB(K), IncX3, IncL...) |
| **Conjugative** | Carries everything needed to transfer itself to another bacterium |
| **Mobilisable** | Can transfer only with help from a conjugative plasmid |
| **Insertion sequence (IS)**, **transposon** | Smaller mobile elements that move genes within and between plasmids (IS26, IS*Ecp1*, Tn*4401*) |

Common resistance plasmids in Enterobacterales: IncFII/IncFIB(K) with *bla*CTX-M-15; IncL with *bla*OXA-48;
IncX3 with *bla*NDM-5.

### 5.2 MOB-suite: hybrid versus short-read

```bash
use_tools plasmids
mkdir -p $RESULTS/05_plasmids

# CHI11E yourself (the hybrid isolate in your practice set), both assemblies
s=CHI11E
mob_recon --infile $COURSE/assemblies/${s}.hybrid.fasta --outdir $RESULTS/05_plasmids/hybrid/$s -n $THREADS
mob_recon --infile $COURSE/assemblies/${s}.fasta --outdir $RESULTS/05_plasmids/shortread/$s -n $THREADS

# CNS1F1 and 1013757 from the trainers (this also brings the short-read results for all isolates, used tomorrow)
use_precomputed 05_plasmids/hybrid 05_plasmids/shortread

for s in CHI11E CNS1F1 1013757; do
    for t in hybrid shortread; do
        echo "== $s $t"
        pick '^size$|^rep_type\(s\)$|^predicted_mobility$|^primary_cluster_id$' $RESULTS/05_plasmids/$t/$s/mobtyper_results.txt
    done
done
```

Fill in this table:

| Isolate | Plasmids found (hybrid) | Plasmids found (short-read) | Same replicons? |
| --- | --- | --- | --- |
| CHI11E | | | |
| CNS1F1 | | | |
| 1013757 | | | |

<details>
<summary><b>What should the hybrid results contain?</b></summary>

The study's closed plasmids (Table S6):

| Isolate | Plasmid sizes (bp) and replicons |
| --- | --- |
| `CHI11E` | 143,302 IncFII(K), IncR · 11,493 · 4,269 |
| `CNS1F1` | 187,438 IncFIB(K)(pCAV1099-114) · 142,041 IncFII(K), IncR · 16,361 |
| `1013757` | 140,482 IncFIB(AP001918), IncFII · 81,013 IncY · 10,638 |

From short reads, MOB-suite usually finds the replicons but splits each plasmid into many contigs, and may merge
contigs from two plasmids, or miss contigs, when they share repeats.

</details>

### 5.3 Which plasmid carries which gene?

```bash
use_tools amr
AMRDB=$SHARED/databases/amrfinderplus/latest
for s in CHI11E CNS1F1 1013757; do
    for p in $RESULTS/05_plasmids/hybrid/$s/plasmid_*.fasta $RESULTS/05_plasmids/hybrid/$s/chromosome.fasta; do
        [ -f "$p" ] || continue
        echo "== $s $(basename $p .fasta)"
        amrfinder -n $p --plus -d $AMRDB --threads $THREADS 2>/dev/null \
            | awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i ~ /symbol/) c=i} NR>1{print "   " $c}'
    done
done
```

<details>
<summary><b>Questions to answer from the output</b></summary>

1. In `CHI11E` and `CNS1F1`, is *bla*CTX-M-15 on a plasmid or the chromosome? Which replicon type?
2. Both `CHI11E` (ST39) and `CNS1F1` (ST45) have a ~142–143 kb IncFII(K)/IncR plasmid. Could it be the same
   plasmid in two different lineages? (You will test this with Mash tomorrow.)
3. In `1013757` (*K. variicola*), which plasmid carries the aerobactin (*iuc*) and salmochelin (*iro*) genes you
   saw on Day 3?

There is no answer key here: these are the questions a surveillance team would ask. Discuss your findings at the
debrief.

</details>

### Check your understanding

<details>
<summary><b>1. A gene sits on a 3 kb contig with an IS26 copy at each end in a short-read assembly. Plasmid or chromosome?</b></summary>

You cannot tell. The contig ends at the IS repeats, so its neighbours are unknown. This is exactly the situation
long reads resolve.

</details>

<details>
<summary><b>2. Why do small plasmids (under 10 kb) sometimes appear in the short-read assembly but not the long-read one?</b></summary>

Nanopore library preparation and read-length filtering favour long DNA molecules, so small plasmids are under-
represented in long reads. Assemblers may then drop them. Always check a hybrid assembly against the short-read
data for small plasmids.

</details>

---

## End-of-day checklist

- [ ] `results/04_amr/genotype_vs_phenotype.tsv` written, and you can explain at least two disagreements.
- [ ] Regimen table with predicted failure percentages, and a sentence on why they are not prevalence.
- [ ] MOB-suite results for the three hybrid isolates, hybrid and short-read.
- [ ] The plasmid table (5.2) is filled in.

## Further reading

- [MOB-suite](https://github.com/phac-nml/mob-suite)
- EUCAST [expert rules and intrinsic resistance tables](https://www.eucast.org/expert_rules_and_expected_phenotypes)

---

[← Day 3](day03.md) · [Course home](README.md) · **Day 4** · [Day 5 →](day05.md)
