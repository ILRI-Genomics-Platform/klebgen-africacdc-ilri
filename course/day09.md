[← Day 8](day08.md) · [Course home](README.md) · **Day 9** · [Day 10 →](day10.md)

# Day 9: Figures, global context and reporting

## At a glance

| Session | Time | Topic | Module |
| --- | --- | --- | --- |
| 1 | 09:00–10:30 | An annotated tree figure in R (ggtree) | 8 |
| 2 | 11:00–12:30 | Your isolates in global context | 8 |
| 3 | 13:30–15:00 | Writing for decision-makers: the one-page brief | 8 |
| 4 | 15:30–17:00 | Sharing data responsibly: ENA submission; peer review of briefs | 8 |

**You will produce:** an annotated tree figure and a peer-reviewed one-page brief.

## Learning objectives

By the end of today you can:

- Draw a tree annotated with metadata and resistance genes.
- Place local isolates in global context.
- Write a one-page brief that leads with the finding and states its limits.
- Describe what is needed to submit reads and metadata to a public archive.

---

## Module 8: Integration, visualisation and reporting (part 2)

### 8.4 An annotated tree with ggtree

```r
library(tidyverse)
library(ape)
library(ggtree)
setwd(Sys.getenv("COURSE"))

master <- read_tsv("results/08_integration/master_table.tsv")
tree   <- read.tree("results/06_phylogeny/ST39/ST39.treefile")
tree   <- drop.tip(tree, "Reference")
tree   <- root(tree, outgroup = "D29438", resolve.root = TRUE)   # the 2005 isolate as outgroup

p <- ggtree(tree) %<+% master +
  geom_tiplab(align = TRUE, size = 3, offset = 0.5) +
  geom_tippoint(aes(colour = ward), size = 2.5) +
  geom_treescale() +
  theme(legend.position = "bottom")

# Heatmap strips beside the tips: one column per feature, values as text
strips <- master |>
  filter(sample %in% tree$tip.label) |>
  transmute(sample,
            Year          = as.character(collection_year),
            Neonate       = neonate,
            ESBL          = if_else(ESBL_or_AmpC, "yes", "no"),
            Carbapenemase = if_else(carbapenemase, "yes", "no"),
            `K locus`     = K_locus) |>
  column_to_rownames("sample")

gheatmap(p, strips, offset = 6, width = 1.2, colnames_angle = 45, colnames_offset_y = -0.5,
         font.size = 3, hjust = 1) +
  scale_fill_viridis_d(option = "D", na.value = "grey90", name = "")

ggsave("results/08_integration/ST39_tree.png", width = 10, height = 6, dpi = 300)
```

> [!TIP]
> If the strips overlap the tip labels, increase `offset`. If the legend is unreadable, plot fewer strips: a
> figure for a report should carry one message.

**Questions for your figure:**

1. Does the cluster form one clade with short branches?
2. Do ward colours line up with the clade?
3. Which isolate breaks the pattern, and what does that suggest?

### 8.5 Global context

Choose one of these approaches:

**A. Pathogenwatch collections (no command line).**

1. Upload your ST39 assemblies to [pathogen.watch](https://pathogen.watch).
2. Open the *K. pneumoniae* public genomes view and filter to ST39.
3. Compare your isolates' LIN codes with ST39 genomes from other countries. Are your cluster isolates closer to
   each other than to any genome from elsewhere?

**B. Public genomes in your own tree (command line).**

1. Your trainers have downloaded a set of public ST39 assemblies to `$SHARED/public_genomes/ST39/`.
2. Add them to the analysis by mapping simulated reads, or use a whole-genome alignment tool such as
   [SKA](https://github.com/bacpop/ska.rust), which works directly from assemblies.
3. Rebuild the tree. Where does the Chatinkha cluster sit?

> [!NOTE]
> Global context answers "is this lineage new here or long-established?" and "is it the same strain causing
> problems elsewhere?". It does not change the local transmission conclusion, which depends on SNP distances and
> epidemiology within the unit.

### 8.6 The one-page brief

A brief for an IPC team, hospital management or a national AMR committee. Use this template:

```markdown
# Klebsiella pneumoniae in the neonatal unit, [facility], [period]

**Key finding.** [One or two sentences: what you found, how confident you are, and why it matters.]

**What we did.** Whole-genome sequencing of [n] bloodstream isolates from [units], [dates].

**Results**
- Lineages: [main STs and counts]
- Resistance: [ESBL, carbapenemase, gentamicin resistance: counts and what they mean for empiric treatment]
- Transmission: [clusters, SNP distances, wards, dates; probable / possible / unlinked]
- [One figure: the annotated tree]

**Limitations.** [Sampling, dates, missing metadata, no environmental samples.]

**Recommended actions.**
1. [Specific action, owner, timeline]
2. [Specific action, owner, timeline]

*Contact: [name]. Data and methods: [link or reference].*
```

> [!IMPORTANT]
> Lead with the finding, not the methods. Keep it to one page. Use plain words: "very closely related (fewer than
> 10 genetic differences)" rather than "< 10 core SNPs after recombination filtering".

**Peer review.** Swap briefs with another group. Check:

- [ ] Could a nurse or hospital manager understand the key finding in 30 seconds?
- [ ] Is every number in the brief traceable to a result file?
- [ ] Are the limitations stated honestly, including the sampling?
- [ ] Are the actions specific and realistic?

### 8.7 Sharing data responsibly

Public sequence data is what made this course possible. Sharing your own surveillance data lets others place it in
context, but it must be done under your institution's agreements.

| Step | What is involved |
| --- | --- |
| 1. Agreements | Data-sharing and ethics approvals cover public release |
| 2. Register | A Webin account at the [European Nucleotide Archive (ENA)](https://www.ebi.ac.uk/ena/browser/submit) |
| 3. Study | Register a study (BioProject) describing the surveillance programme |
| 4. Samples | Register samples with a pathogen checklist: species, collection date, country, host, isolation source |
| 5. Reads | Upload reads and link them to samples (Webin-CLI or the web interface) |
| 6. Record accessions | Put the run and sample accessions in your metadata and publications |

> [!WARNING]
> Public sample metadata must not identify patients. Use country and year (or month) rather than exact date and
> facility, unless your agreements allow more.

---

## End-of-day checklist

- [ ] `results/08_integration/ST39_tree.png` saved.
- [ ] You have compared your isolates with genomes from elsewhere (Pathogenwatch or your own tree).
- [ ] One-page brief written and peer-reviewed.
- [ ] You can list the steps and metadata needed for an ENA submission.

## Further reading

- [ggtree book](https://yulab-smu.top/treedata-book/), [Pathogenwatch](https://pathogen.watch),
  [ENA submission guide](https://ena-docs.readthedocs.io)

---

[← Day 8](day08.md) · [Course home](README.md) · **Day 9** · [Day 10 →](day10.md)
