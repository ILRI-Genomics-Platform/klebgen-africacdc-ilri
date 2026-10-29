[← Day 6](day06.md) · [Course home](README.md) · **Day 7** · [Day 8 →](day08.md)

# Day 7: Is it an outbreak? Clusters, timing and action

## At a glance

| Session | Time | Topic | Module |
| --- | --- | --- | --- |
| 1 | 09:00–10:30 | Tree plus metadata: classify clusters (Microreact) | 6 |
| 2 | 11:00–12:30 | When did the cluster start? Time-scaled trees (BactDating) | 6 |
| 3 | 13:30–15:00 | Writing the outbreak summary | 6 |
| 4 | 15:30–17:00 | Role play: presenting to an infection prevention and control (IPC) team | 6 |

**You will produce:** a one-page outbreak investigation summary with a tree figure.

## Learning objectives

By the end of today you can:

- Combine the tree, SNP distances, ward and dates to classify each cluster.
- Test for temporal signal and estimate when a cluster began, and explain the limits of that estimate.
- Communicate a genomic finding to an IPC team in plain language, with clear actions.

---

## Module 6: Phylogenetics and outbreak detection (part 2)

### 6.10 Tree and metadata together in Microreact

Microreact shows a tree, a timeline and a table side by side. It needs a metadata file in CSV with an `id`
column that matches the tree's tip names.

```bash
cd $RESULTS/06_phylogeny/ST39
# Metadata as CSV: first column renamed to 'id', year column to 'year', every field quoted
# (some K-locus values contain commas)
awk -F'\t' 'BEGIN{OFS=","} NR==1 {$1="id"; $2="year"} {for(i=1;i<=NF;i++) $i="\""$i"\""; print}' \
    $COURSE/training_data/metadata_compact.tsv > microreact_metadata.csv
head -3 microreact_metadata.csv
```

1. Go to [microreact.org](https://microreact.org) and choose **Upload**.
2. Drop in `ST39.treefile` and `microreact_metadata.csv`.
3. Colour the tree by `ward`, then by `year`. Add a timeline using `year`.

> [!NOTE]
> The tree has a tip called `Reference` (from Snippy). It has no metadata row; that is fine.

### 6.11 Classify the cluster

Use this framework for each group of closely related isolates:

| Evidence | Probable transmission | Possible transmission | Unlinked |
| --- | --- | --- | --- |
| Core SNP distance | Within threshold | Near threshold | Well above threshold |
| Place | Same ward | Same hospital, different ward | Different facilities |
| Time | Overlapping or close dates | Months apart | Years apart |

Fill this in for the ST39 group:

| Isolate(s) | SNPs to nearest cluster member | Ward | Year | Classification | Reason |
| --- | --- | --- | --- | --- | --- |
| Cluster (8) | | Chatinkha Nursery | 2017 | | |
| `BKR8JH` | | Chatinkha Nursery | 2016 | | |
| `D29438` | | unknown | 2005 | | |
| `BHB7VE` | | 4A Female Medical Ward | 2016 | | |
| `CHI11E` | | unknown | 2019 | | |

Before you settle a classification, check that it does not depend on the method: look up the same pairs in the
Snippy and SKA2 comparison you made on Day 6 (`$RESULTS/06_phylogeny/snippy_vs_ska2_pairs.tsv`). A pair that is
close by one method and not the other is at most **possible** transmission until you know why.

<details>
<summary><b>Why is <code>BKR8JH</code> the interesting one?</b></summary>

It is from the same neonatal ward, one year earlier. If it falls within a few SNPs of the 2017 cluster, the lineage
may have persisted in the unit for a year or more, which points towards a **reservoir** (sinks, equipment,
colonised staff or babies) rather than one short outbreak. If it is far away, the 2017 cluster was probably a new
introduction. Note also its different capsule locus (KL23, versus KL2 for the cluster) in the published typing.
What does that suggest before you even look at the SNPs?

</details>

### 6.12 When did the cluster start? BactDating

First check whether there is enough **temporal signal**: do isolates collected later have accumulated more
mutations from the root?

```bash
use_tools r
cd $RESULTS/06_phylogeny/ST39
R
```

```r
source(file.path(Sys.getenv("COURSE"), "bin/r_setup.R"))   # your own library, and install-on-demand
use_packages(cran = c("ape", "readr", "coda"), github = c(BactDating = "xavierdidelot/BactDating"))

library(ape)
library(BactDating)
library(readr)

meta <- read_tsv(file.path(Sys.getenv("COURSE"), "training_data/metadata_compact.tsv"))

# Gubbins tree: branch lengths in substitutions, plus the share of each branch NOT affected by recombination
# (reads ST39.final_tree.tre, ST39.node_labelled.final_tree.tre and ST39.per_branch_statistics.csv)
tree <- loadGubbins("ST39")
tree <- drop.tip.useRec(tree, "Reference")      # the reference has no date; drop.tip.useRec keeps the
                                                # recombination information in step (ape's drop.tip would not)

dates <- meta$collection_year[match(tree$tip.label, meta$sample)]
names(dates) <- tree$tip.label
dates

# 1. Temporal signal: root the tree where dates fit best, then root-to-tip regression
#    (the p-value comes from randomly permuting the dates)
rooted <- initRoot(tree, dates, useRec = TRUE)
roottotip(rooted, dates)

# 2. Time-scaled tree, accounting for recombination (this takes a few minutes)
res <- bactdate(tree, dates, nbIts = 1e5, useRec = TRUE, showProgress = TRUE)
res                                   # root date and substitution rate, with credible intervals
plot(res, "treeCI")
plot(res, "trace")                    # the chains should look like flat, hairy caterpillars
```

> [!TIP]
> For your own data, `bin/run_bactdating.R` in the pipeline does all of this for one group from the pipeline's
> Gubbins output, treats year-only dates as a whole year rather than a point, checks convergence, and dates the
> common ancestor of chosen isolates (`--mrca`). See the pipeline README.

> [!WARNING]
> Collection dates here are **years only**, and most isolates are from 2016–2019. Expect wide credible intervals.
> If root-to-tip regression shows no clear positive slope, do not report a dated tree: say that the data cannot
> date the cluster. Honest uncertainty is a valid result.

> [!NOTE]
> `loadGubbins()` and `useRec = TRUE` both need Gubbins output. If the pipeline ran with `--skip_gubbins` (or
> Gubbins was skipped for this group because too few isolates had enough data), there is none: read the tree
> with `ape::read.tree()` instead and run `bactdate(tree, dates, useRec = FALSE)`. The rate and root date are
> then **recombination-naive** — recombination imports many substitutions at once, so an uncorrected rate is
> usually too high and the root date too recent. Report that limitation with the estimate.

<details>
<summary><b>What does the date of the cluster's most recent common ancestor tell an IPC team?</b></summary>

Roughly when the lineage started diversifying in the unit. If that date is well before the first detected case,
the lineage may have been circulating undetected (colonised babies, an environmental reservoir) before infections
were noticed. With year-only dates, treat it as "around year X", not a month.

</details>

### 6.13 Write the outbreak summary

One page, in this order:

1. **Key finding** (one or two sentences): what you found and how confident you are.
2. **Isolates**: how many, which ward, which period.
3. **Genomic evidence**: ST, SNP distances within and outside the cluster, tree figure.
4. **Interpretation**: probable or possible transmission, reservoir or new introduction.
5. **Limitations**: sampling (only cultured, sequenced isolates), year-only dates, no environmental samples.
6. **Recommended actions**: two or three, specific.

> [!TIP]
> Write "probable transmission within the neonatal unit" rather than "outbreak confirmed". Genomics supports or
> weakens an epidemiological hypothesis; it rarely settles it alone.

---

## Role play: the IPC meeting

**Scenario.** You are the genomics lead. The neonatal unit's IPC team has asked whether the 2017 *Klebsiella*
bloodstream infections are linked. You have 5 minutes to present and 10 minutes for questions.

**Roles** (groups of 4): genomics lead, IPC nurse, neonatologist, microbiology lab manager.

**The panel will ask:**

- "So is it an outbreak or not?"
- "Did it come from one baby, or from the unit itself?"
- "What should we do on Monday morning?"
- "Will treatment change?" (Use your Day 4 regimen results: these isolates carry CTX-M-15 and gentamicin
  resistance genes.)

**Agree on two actions**, for example: environmental sampling of sinks and incubators; a hand-hygiene and
cleaning audit; screening of newly admitted babies; review of empiric antibiotic guidance.

> [!NOTE]
> Trainers: model answers for each of the four questions, the weak answers to watch for, and how to mark the
> session are in `course/TRAINER_GUIDE.md`, section 9.

---

## End-of-day checklist

- [ ] Microreact project with the ST39 tree, coloured by ward and year.
- [ ] Classification table completed for the ST39 group.
- [ ] Root-to-tip result recorded; dated tree reported only if the temporal signal supports it.
- [ ] One-page outbreak summary written.
- [ ] Two actions agreed in the role play.

## Further reading

- Didelot X *et al.* Bayesian inference of ancestral dates on bacterial phylogenetic trees. *Nucleic Acids Res*
  2018. [BactDating](https://github.com/xavierdidelot/BactDating)
- [Microreact documentation](https://docs.microreact.org)

---

[← Day 6](day06.md) · [Course home](README.md) · **Day 7** · [Day 8 →](day08.md)
