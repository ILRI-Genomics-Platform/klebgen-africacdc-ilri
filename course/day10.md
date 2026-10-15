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

> A neonatal surveillance network has sent you **five *K. pneumoniae* bloodstream isolates** from the same
> hospital, collected between **2005 and 2019**, sequenced with Illumina; one also with Nanopore. The neonatal
> nursery (**Chatinkha**) had a run of *Klebsiella* infections in **2017**, and the hospital wants to know
> whether isolates from other years and wards are linked to it, and what that means for treatment.
>
> The national AMR committee meets this afternoon. It wants answers to four questions.

1. **Is there transmission?** Is there one cluster, several, or only background? In which ward and year?
2. **Is resistance spreading by clone or by plasmid?** Do the isolates carry the same resistance plasmids, or different ones?
3. **Will empiric treatment work?** Which regimens are likely to fail, and how confident are you given the AST?
4. **What should the unit do next?** Two or three specific actions.

## Your data

The five capstone isolates are `CHI11E`, `BHB7VE`, `CAA65Z`, `CAA99X` and `D29438`: the pipeline's test set.
You met them during the week among the 29 training isolates; today you analyse them **end to end**, from raw
reads to a recommendation, with the pipeline, and bring all the evidence together yourselves. `CHI11E` also has
Nanopore reads, so it gets a hybrid assembly.

| File | Contents |
| --- | --- |
| `$SHARED/capstone/samplesheet.csv` | Sample sheet for the pipeline (5 isolates; the pipeline groups them by Kleborate species + ST) |
| `$SHARED/capstone/metadata_capstone.tsv` | Year, ward, neonate, published typing |
| `$SHARED/capstone/ast_capstone.tsv` | Disc-diffusion AST |
| `$SHARED/capstone/precomputed/` | Full pipeline results, if your own run does not finish in time |

> [!TIP]
> Start the pipeline in the first 15 minutes, then work from the pre-computed results while it runs. Spend your
> time on interpretation, not on waiting. The hybrid assembly of `CHI11E` takes longest (1–3 hours).
>
> With five isolates, every number counts: say how confident you are, and what more sampling would change.

```bash
use_tools nextflow
mkdir -p $COURSE/capstone && cd $COURSE/capstone
nextflow run $SHARED/neonatal-amr-nf/main.nf -profile singularity -c $SITE_CONFIG \
    --input $SHARED/capstone/samplesheet.csv --outdir results \
    --kraken2_db $SHARED/databases/kraken2 \
    --amrfinderplus_db $SHARED/databases/amrfinderplus/amrfinderdb.tar.gz
```

> [!TIP]
> **Disconnected during the run?** The run stops with your session, but nothing finished is lost. Log in,
> `start_session`, `cd $COURSE/capstone`, and run the same command again with `-resume` at the end:
> Nextflow reuses every finished step and only runs the rest.

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
- [ ] Plasmids: which resistance genes are on plasmids, and are the same plasmids found in several isolates?
- [ ] Transmission: a tree and SNP distances for each group (e.g. `Kpneumoniae_ST39`); each cluster classified with ward and year.
- [ ] Robustness: does each cluster call hold with both Snippy and SKA2? (`bin/compare_snp_methods.py --phylo $SHARED/capstone/precomputed/06_phylogeny --out $COURSE/capstone/snippy_vs_ska2`)
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
