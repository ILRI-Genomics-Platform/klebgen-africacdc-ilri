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

A published tree figure is a tree plus *rings* of metadata outside the tips, each with its own colour scale, and
the clades of interest shaded and labelled. Three packages do the work: `ggtree` draws the tree, `ggtreeExtra`
adds the rings (`geom_fruit()`), and `ggnewscale` lets one plot carry more than one `fill` scale
(`new_scale_fill()`) — without it, the second ring would overwrite the first ring's colours.

> [!IMPORTANT]
> `ggtreeExtra` is new to the course stack. If you installed the R packages before today you will get
> `there is no package called 'ggtreeExtra'`. Install it once, **on a login node** (compute nodes have no
> internet):
>
> ```bash
> cd $COURSE && bin/run_R.sh install_r_packages.R
> ```
>
> That re-runs the whole list but only fetches what is missing, so it is quick the second time. The preamble
> below also installs it on demand, which is the belt to that braces.

```r
setwd(Sys.getenv("COURSE"))
source("bin/r_setup.R")      # your own library, and install-on-demand
use_packages(cran = c("tidyverse", "ape", "ggnewscale"),
             bioc = c("ggtree", "ggtreeExtra"))

library(tidyverse)
library(ape)
library(ggtree)
library(ggtreeExtra)   # geom_fruit(): rings of data outside the tips
library(ggnewscale)    # new_scale_fill(): a second, third, fourth fill scale on one plot
source("bin/dviz_theme.R")   # dviz_font_family, already checked against the fonts installed here

master <- read_tsv("results/08_integration/master_table.tsv")
tree   <- read.tree("results/06_phylogeny/ST39/ST39.treefile")
tree   <- drop.tip(tree, "Reference")
tree   <- root(tree, outgroup = "D29438", resolve.root = TRUE)   # the 2005 isolate as outgroup

span <- max(node.depth.edgelength(tree))   # tree depth, for offsets below
```

**Step 1 — one annotation row per tip, in the tree's own order.** Start from the tip labels, not from the
master table: that guarantees one row per tip even when an isolate is missing from the metadata, and it puts
the rows in the order `ggtree` expects.

```r
ann <- tibble(sample = tree$tip.label) |>
  left_join(master, by = "sample") |>
  mutate(
    # 15 K loci means 15 colours and no message. Keep those seen >= 3 times, pool the rest.
    K_group = fct_lump_min(factor(replace_na(K_locus, "unknown")), min = 3, other_level = "Other"),
    # One ring instead of two yes/no rings: the strongest beta-lactam mechanism found
    Betalactam = case_when(
      replace_na(carbapenemase, FALSE) ~ "Carbapenemase",
      replace_na(ESBL_or_AmpC,  FALSE) ~ "ESBL / AmpC",
      TRUE                             ~ "Neither"
    ) |> factor(levels = c("Carbapenemase", "ESBL / AmpC", "Neither")),
    Year = factor(collection_year)
  )

missing <- ann$sample[is.na(ann$ST)]
if (length(missing)) message("no metadata for: ", paste(missing, collapse = ", "))
```

**Step 2 — the clades to shade.** Find the most recent common ancestor of each ST with at least two tips.
`getMRCA()` takes tip *labels* and returns a node *number*, which is what `geom_highlight()` and
`geom_cladelab()` want.

```r
clade_sts <- ann |> filter(!is.na(ST)) |> count(ST) |> filter(n >= 2) |> pull(ST)

clades <- tibble(ST = clade_sts, clade = paste0("ST", clade_sts)) |>
  mutate(node = map_int(ST, \(st) as.integer(getMRCA(tree, ann$sample[ann$ST == st & !is.na(ann$ST)]))))

clade_cols <- setNames(
  colorRampPalette(c("lightblue3", "lightpink3", "mediumpurple2"))(nrow(clades)),
  clades$clade)
clades
```

**Step 3 — the tree, then the rings.** Each ring is `new_scale_fill()` followed by `geom_fruit()` and the
scale for that ring. `offset` and `pwidth` in `geom_fruit()` are fractions of the tree's width; the offsets in
`geom_tiplab()` and `geom_cladelab()` are in branch-length units, which is why they are written as multiples of
`span`.

```r
p <- ggtree(tree, layout = "rectangular", linewidth = 0.5) %<+% ann +
  geom_highlight(data = clades, aes(node = node, fill = clade),
                 alpha = 0.15, align = "right", extend = 0.6 * span, show.legend = FALSE) +
  geom_cladelab(data = clades, mapping = aes(node = node, label = clade),
                barcolour = clade_cols[clades$clade], barsize = 1.5,
                align = TRUE, offset = 0.45 * span, fontsize = 3) +
  geom_tippoint(aes(colour = ward), size = 2.5) +
  geom_tiplab(align = TRUE, linesize = 0.2, size = 2.6, offset = 0.03 * span) +
  geom_treescale(fontsize = 2.6) +
  scale_colour_brewer(palette = "Dark2", na.value = "grey70", name = "Ward") +
  scale_fill_manual(values = clade_cols, breaks = clades$clade, guide = "none")

fig <- p +
  new_scale_fill() +
  geom_fruit(geom = geom_tile, mapping = aes(y = sample, fill = Betalactam),
             offset = 0.10, pwidth = 0.08, alpha = 0.9, linewidth = 0,
             position = position_identityx()) +
  scale_fill_manual(name = "Beta-lactam",
                    values = c("Carbapenemase" = "#b2182b",
                               "ESBL / AmpC"   = "#ef8a62",
                               "Neither"       = "#f0f0f0"),
                    na.value = "grey90",
                    guide = guide_legend(order = 2, override.aes = list(size = 4))) +
  new_scale_fill() +
  geom_fruit(geom = geom_tile, mapping = aes(y = sample, fill = K_group),
             offset = 0.05, pwidth = 0.08, alpha = 0.9, linewidth = 0,
             position = position_identityx()) +
  scale_fill_brewer(name = "K locus", palette = "Set3", na.value = "grey90",
                    guide = guide_legend(order = 3, ncol = 2)) +
  new_scale_fill() +
  geom_fruit(geom = geom_tile, mapping = aes(y = sample, fill = Year),
             offset = 0.05, pwidth = 0.08, alpha = 0.9, linewidth = 0,
             position = position_identityx()) +
  scale_fill_viridis_d(name = "Year", option = "D", na.value = "grey90",
                       guide = guide_legend(order = 4)) +
  theme_tree() +
  theme(text              = element_text(family = dviz_font_family),
        legend.position   = "right",
        legend.background = element_blank(),
        legend.key        = element_blank(),
        legend.title      = element_text(size = 10),
        legend.text       = element_text(size = 9)) +
  guides(colour = guide_legend(order = 1, ncol = 1, override.aes = list(size = 3)))
fig

dir.create("results/08_integration", showWarnings = FALSE, recursive = TRUE)
ggsave("results/08_integration/ST39_tree.png", fig, width = 11, height = 7, dpi = 300)
ggsave("results/08_integration/ST39_tree.pdf", fig, width = 11, height = 7,
       device = if (capabilities("cairo")) cairo_pdf else pdf)
```

> [!TIP]
> Tuning the figure, in the order worth trying:
> - Rings overlapping the tip labels, or each other: raise the first ring's `offset`; the later ones are
>   measured from the ring before them.
> - Rings too thin or too fat: `pwidth`.
> - Clade labels off the edge of the panel: lower `offset` in `geom_cladelab()`, or add
>   `+ xlim(0, 2 * span)`.
> - Legend crowding the tree: `ncol` inside `guide_legend()`, or drop a ring. A figure for a report should
>   carry one message; three rings is usually the limit before it carries none.

> [!WARNING]
> `getMRCA()` returns the ancestor of the tips you give it, which is the clade you want **only if that ST is
> monophyletic**. If an ST is scattered across the tree, the highlight will swallow unrelated isolates. Check
> each shaded block against the tip labels before the figure leaves your screen, and drop any ST that does not
> hold together — then say so in the caption.

> [!NOTE]
> Column names here are the ones built by hand on Day 8 (`ESBL_or_AmpC`, `carbapenemase`, `K_locus`). The
> master table written by `bin/aggregate_results.py` / `.R` uses the aggregator's own names — `esbl_or_ampc`,
> `carbapenemase`, `kleborate_K_locus`, `kaptive_K_locus`. Run `names(master)` first and adjust.

**Questions for your figure:**

1. Does the cluster form one clade with short branches?
2. Do ward colours line up with the clade?
3. Do the resistance and K locus rings change at the same point in the tree, or independently? A ring that
   changes *within* a clade is a gene moving between related isolates; a ring that changes *between* clades is
   two separate lineages.
4. Which isolate breaks the pattern, and what does that suggest?

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
