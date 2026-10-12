[← Day 9](day09.md) · [Course home](README.md) · **Day 10**

# Day 10: Capstone

## At a glance

| Session | Time | Activity |
| --- | --- | --- |
| 1 | 09:00–10:30 | Briefing; start the analysis in teams |
| 2 | 11:00–12:30 | Analysis |
| 3 | 13:30–15:00 | Analysis; prepare report and slides |
| 4 | 15:30–17:00 | Team presentations to a mock national AMR committee; feedback; close |

**You will produce:** a capstone report and a 10-minute presentation, as a team of 3–4.

---

## The scenario

> A neonatal surveillance network has sent you **36 bloodstream isolates** from the same hospital, sequenced with
> Illumina. The paediatric high-dependency unit (**PICU/HDU**) reported a rise in *Klebsiella* infections in
> **2019**, and the neonatal nursery wants to know whether its lineages from earlier years have come back.
>
> The national AMR committee meets this afternoon. It wants answers to four questions.

1. **Is there transmission?** Is there one cluster, several, or only background? In which ward and year?
2. **Is resistance spreading by clone or by plasmid?** Are any resistance plasmids shared between STs?
3. **Will empiric treatment work?** Which regimens are likely to fail, and how confident are you given the AST?
4. **What should the unit do next?** Two or three specific actions.

## Your data

The capstone isolates are new to you. They include more isolates of lineages you have met (ST39, ST35), lineages
you have not, and some deliberately difficult samples.

| File | Contents |
| --- | --- |
| `$SHARED/capstone/samplesheet.csv` | Sample sheet for the pipeline (`phylo_group` filled for ST14, ST35, ST39; the pipeline groups by Kleborate species + ST unless you add `--phylo_grouping samplesheet`) |
| `$SHARED/capstone/metadata_capstone.tsv` | Year, ward, neonate, published typing |
| `$SHARED/capstone/ast_capstone.tsv` | Disc-diffusion AST |
| `$SHARED/capstone/precomputed/` | Full pipeline results, if your own run does not finish in time |

> [!TIP]
> Start the pipeline in the first 15 minutes, then work from the pre-computed results while it runs. Spend your
> time on interpretation, not on waiting.

```bash
conda activate nextflow
mkdir -p $COURSE/capstone && cd $COURSE/capstone
nextflow run $SHARED/neonatal-amr-nf/main.nf -profile singularity -c $SITE_CONFIG \
    --input $SHARED/capstone/samplesheet.csv --outdir results \
    --kraken2_db $SHARED/databases/kraken2 \
    --amrfinderplus_db $SHARED/databases/amrfinderplus/amrfinderdb.tar.gz
```

## Suggested team roles

| Role | Responsible for |
| --- | --- |
| QC and assembly lead | Which isolates are fit for which analysis (Modules 1–2) |
| Typing and AMR lead | Lineages, resistance, regimens, genotype–phenotype (Modules 3–4) |
| Transmission lead | Trees, SNP distances, cluster classification (Modules 5–6) |
| Reporting lead | Master table, figures, report and slides (Modules 7–8) |

Everyone checks everyone else's numbers before the presentation.

## Checklist for the analysis

- [ ] QC: which isolates did you exclude, and why?
- [ ] Species: are all isolates *K. pneumoniae*?
- [ ] Lineages: ST counts, by ward and year.
- [ ] Resistance: ESBL, carbapenemase and aminoglycoside genes; agreement with AST.
- [ ] Plasmids: any clusters shared between STs?
- [ ] Transmission: a tree and SNP distances for each group (e.g. `Kpneumoniae_ST39`); each cluster classified with ward and year.
- [ ] Antigens: K and O types of any cluster.
- [ ] Master table: one row per isolate, row count checked.

## What to hand in

1. **A report** of at most two pages, using the Day 9 brief template, plus an appendix with your master table.
2. **A 10-minute presentation**: key finding first, one tree figure, one resistance table, actions last.

## Assessment

| Component | Weight | What we look for |
| --- | --- | --- |
| Module deliverables (Days 1–9) | 50% | Completed checklists and outputs |
| Capstone report | 30% | Correct analysis, honest limitations, traceable numbers |
| Capstone presentation | 20% | Clear finding, sound reasoning, specific actions, good answers to questions |

A participant passes with **60% or more** and a completed pipeline run.

### Presentation scoring guide

| | Excellent | Adequate | Needs work |
| --- | --- | --- | --- |
| **Finding** | Stated first, with confidence level | Present but buried | Unclear |
| **Evidence** | SNP distances, ward and dates combined | Genomics only | Over-claims ("outbreak confirmed") |
| **Resistance** | Regimen implications, with AST caveats | Gene list only | Missing |
| **Limitations** | Sampling, dates, data gaps all named | Some named | None |
| **Actions** | Specific, owned, realistic | Generic | Missing |

---

## After the course

- Keep your `$COURSE` folder: every command you ran is reusable on your own data.
- The pipeline (`neonatal-amr-nf`) runs on your own sample sheets; see its README.
- Before using these methods for routine surveillance, validate them on reference isolates with known results.

Thank you for taking part.

---

[← Day 9](day09.md) · [Course home](README.md) · **Day 10**
