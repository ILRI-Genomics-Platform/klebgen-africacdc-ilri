[← Day 4](day04.md) · [Course home](README.md) · **Day 5** · [Day 6 →](day06.md)

# Day 5: Plasmids across lineages, and the whole analysis as a pipeline

## At a glance

| Session | Time | Topic | Module |
| --- | --- | --- | --- |
| 1 | 09:00–10:30 | Plasmid clusters: is the same plasmid in different lineages? | 5 |
| 2 | 11:00–12:30 | Genetic context of a resistance gene (Bakta annotation) | 5 |
| 3 | 13:30–15:00 | The course pipeline: structure, test run | Pipeline |
| 4 | 15:30–17:00 | Run the pipeline on the training set; read its reports | Pipeline |

**You will produce:** a plasmid-sharing table and figure, and a complete pipeline run.

## Learning objectives

By the end of today you can:

- Group plasmids into clusters and find plasmids shared between lineages or species.
- Annotate a plasmid and describe the genetic context of a resistance gene.
- Explain how the course pipeline is organised (nf-core modules, local modules, subworkflows, configuration).
- Run the pipeline, resume it after a change, and find each result.

---

## Module 5: Plasmids and mobile genetic elements (part 2)

### 5.4 Plasmid clusters across all isolates

MOB-suite assigns every reconstructed plasmid to a **cluster** (`primary_cluster_id`) of closely related plasmids.
Two isolates with the same cluster ID probably carry the same plasmid, even if they belong to different STs.

```bash
use_tools plasmids
# Short-read MOB-suite results for all 29 isolates (pre-computed; or run mob_recon on $COURSE/best/*.fasta)
cp -rn $SHARED/precomputed/05_plasmids/shortread/* $RESULTS/05_plasmids/shortread/

# One table: sample, cluster, replicons, size, mobility
for f in $RESULTS/05_plasmids/shortread/*/mobtyper_results.txt; do
    s=$(basename $(dirname $f))
    pick '^primary_cluster_id$|^rep_type\(s\)$|^size$|^predicted_mobility$' $f | tail -n +2 \
        | awk -v s=$s 'BEGIN{OFS="\t"} {$1=s; print}'
done > $RESULTS/05_plasmids/all_plasmids.tsv

# Which plasmid clusters occur in more than one ST?
join -t$'\t' -1 1 -2 1 \
    <(sort -k1,1 $RESULTS/05_plasmids/all_plasmids.tsv) \
    <(cut -f1,7 $COURSE/training_data/metadata_compact.tsv | tail -n +2 | sort -k1,1) \
  | awk -F'\t' 'BEGIN{OFS="\t"} {st=$NF; NF--; print st, $0}' | sort -t$'\t' -k6,6 \
  > $RESULTS/05_plasmids/plasmids_by_cluster.tsv
less -S $RESULTS/05_plasmids/plasmids_by_cluster.tsv
```

> [!TIP]
> Columns: ST, sample, size, replicons, mobility, `primary_cluster_id`. Rows are sorted by cluster, so plasmids
> in the same cluster sit together. The question is which clusters appear with **more than one ST**.

### 5.5 Are the two IncFII(K)/IncR plasmids the same?

Yesterday you saw that `CHI11E` (ST39) and `CNS1F1` (ST45) each carry a ~142–143 kb IncFII(K)/IncR plasmid.
Compare them with Mash:

```bash
mkdir -p $RESULTS/05_plasmids/mash && cd $RESULTS/05_plasmids/mash
for d in $RESULTS/05_plasmids/hybrid/*/; do      # prefix each plasmid file with its isolate
    s=$(basename $d)
    for p in $d/plasmid_*.fasta; do cp $p ${s}_$(basename $p); done
done

mash sketch -o plasmids *.fasta
mash dist plasmids.msh plasmids.msh | awk '$1 < $2' | sort -k3,3g > plasmid_distances.tsv
head plasmid_distances.tsv
cd $COURSE
```

<details>
<summary><b>How to read a Mash distance</b></summary>

Mash distance approximates the proportion of differing bases: **0** means identical k-mer content, **0.01** is
roughly 99% identity. MOB-suite uses a Mash distance of **0.06** as its cut-off for grouping plasmids into the same
cluster. A very small distance between plasmids from two STs suggests the plasmid has moved between lineages; a
full answer needs an alignment of the two plasmids.

</details>

### 5.6 Genetic context: annotate a plasmid with Bakta

```bash
use_tools plasmids
s=CHI11E
# Find the plasmid file carrying blaCTX-M-15 (you found it on Day 4, section 5.3)
P=$RESULTS/05_plasmids/hybrid/$s/plasmid_XXXX.fasta          # <- replace XXXX

cat $SHARED/databases/bakta/db-light/version.json        # Bakta light database, version 6
bakta --db $SHARED/databases/bakta/db-light --output $RESULTS/05_plasmids/bakta_$s \
    --prefix ${s}_ctxm_plasmid --threads $THREADS $P

# Genes within 10 kb either side of blaCTX-M-15
T=$RESULTS/05_plasmids/bakta_$s/${s}_ctxm_plasmid.tsv
pos=$(grep -i "CTX-M" $T | head -1 | cut -f3)
awk -F'\t' -v p=$pos '!/^#/ && $3 > p-10000 && $4 < p+10000 {print $3"\t"$4"\t"$5"\t"$7"\t"$8}' $T
```

Look for insertion sequences (IS*Ecp1*, IS26) and other resistance genes next to *bla*CTX-M-15.

> [!TIP]
> For a picture of the region, upload the Bakta GenBank file (`.gbff`) to [Proksee](https://proksee.ca) or draw it
> with [clinker](https://github.com/gamcil/clinker) for several plasmids side by side.

### Check your understanding

<details>
<summary><b>1. Two isolates of different STs share a plasmid cluster. What does that mean for an outbreak investigation?</b></summary>

There may be a **plasmid outbreak**: the resistance is spreading by plasmid transfer between bacteria, not by one
clone. Typing by ST alone would miss it, and infection-control actions would need to target the shared
environment and practices rather than one strain.

</details>

<details>
<summary><b>2. Why is the genetic context (the genes around <i>bla</i>CTX-M-15) useful?</b></summary>

The same gene can sit in different contexts. An identical context (same IS element, same neighbouring genes, same
plasmid backbone) in different isolates is evidence of a shared mobile element. Different contexts suggest
independent acquisitions.

</details>

---

## The course pipeline: neonatal-amr-nf

### Structure

Everything you have run by hand so far is in one Nextflow pipeline:

```text
neonatal-amr-nf/
├── main.nf                       entry point
├── nextflow.config               parameters and profiles
├── conf/modules.config           tool arguments and output folders  <- you edit this
├── workflows/neonatal_amr.nf     the analysis
├── subworkflows/local/           read_qc, assembly, typing_amr, phylogeny
└── modules/
    ├── nf-core/                  19 community modules (fastp, shovill, amrfinderplus, gubbins...)
    └── local/                    5 course modules (kleborate v3, autocycler, polypolish, ezclermont, snippy/clean)
```

Each module folder has three files: `main.nf` (the process), `meta.yml` (inputs and outputs, documented) and
`environment.yml` (the exact software version).

```bash
cd $SHARED/neonatal-amr-nf
cat modules/local/kleborate/main.nf
cat modules/local/kleborate/environment.yml
```

### The software: shared Singularity images

On this server the pipeline runs every tool inside a **Singularity image**: a single file holding the tool and
everything it needs, at a fixed version. Your trainers prepared all the images in advance in
`$SHARED/containers/`. The folder is **read-only**: you use the images, you never download, build or change them.

```bash
ls $SHARED/containers/nf-core | head -5        # images for the nf-core modules
ls $SHARED/containers/local                   # images built for the course's own modules
cut -f1,2,5 $SHARED/containers/MANIFEST.tsv | less -S       # where each image came from
```

These are the images you have used all week: `use_tools` runs each command inside its image. Run one
directly, the way `use_tools` does it for you:

```bash
singularity exec $SHARED/containers/local/kleborate-3.2.4.sif kleborate --version
```

Two options choose how the pipeline runs on this server:

| Option | What it does |
| --- | --- |
| `-profile singularity` | Run every tool inside its Singularity image |
| `-c $SITE_CONFIG` | This server's settings: where the shared images are, and (on the HPC) which Slurm queue and nodes each step uses |

```bash
less $SITE_CONFIG
```

> [!NOTE]
> With these options the pipeline uses only the shared images. If one is missing, the run stops before it
> starts, with a message saying so: tell a trainer. You cannot fix it yourself, by design.

Every pipeline has a help page. It needs no other options:

```bash
nextflow run $SHARED/neonatal-amr-nf/main.nf --help
```

### Run 1: a test run (no tools executed)

```bash
use_tools nextflow
mkdir -p $COURSE/pipeline && cd $COURSE/pipeline
nextflow run $SHARED/neonatal-amr-nf/main.nf -profile test,singularity -c $SITE_CONFIG -stub \
    --outdir test_results
```

Stub mode checks that every step connects, using placeholder files. Look at what it created:

```bash
find test_results -maxdepth 2 | sort
```

### Run 2: the training set

```bash
cd $COURSE/pipeline
nextflow run $SHARED/neonatal-amr-nf/main.nf -profile singularity -c $SITE_CONFIG \
    --input $SHARED/samplesheet.csv \
    --outdir results \
    --kraken2_db $SHARED/databases/kraken2 \
    --amrfinderplus_db $SHARED/databases/amrfinderplus/amrfinderdb.tar.gz
```

> [!NOTE]
> The full run takes several hours on one machine. Your trainers ran it in advance; after a few minutes, stop
> yours with `Ctrl+C` and look at the pre-computed results in `$SHARED/precomputed/pipeline_results/`.

### Reading the outputs

| Folder | Contents | Day you did this by hand |
| --- | --- | --- |
| `01_qc/` | fastp, Filtlong, NanoPlot, Kraken2, Bracken | 1 |
| `02_assembly/` | Shovill, Autocycler, hybrid, QUAST | 2 |
| `03_typing/` | Kleborate (and *E. coli* typing) | 3 |
| `04_amr/` | AMRFinderPlus | 3–4 |
| `05_plasmids/` | MOB-suite | 4–5 |
| `06_phylogeny/` | One tree and SNP-distance matrix per ST group | 6–7 |
| `multiqc/` | Summary report | |
| `pipeline_info/` | Software versions, run time per step | |

Open `pipeline_info/execution_report_*.html`: which step took the longest? Which used the most memory?

### Exercise: change a setting and resume

Tool settings live in `conf/modules.config`. To change one, do not edit the pipeline: write a small config
file with only the setting you want to change, and pass it with `-c`.

1. Make a config that raises the fastp minimum read length from 50 to 75:

    ```bash
    cat > my_settings.config <<'CFG'
    process {
        withName: 'FASTP' {
            ext.args = '--qualified_quality_phred 20 --length_required 75'
        }
    }
    CFG
    ```

2. Re-run the test with `-resume` and your config:

    ```bash
    nextflow run $SHARED/neonatal-amr-nf/main.nf -profile test,singularity -c $SITE_CONFIG -stub \
        --outdir test_results -c my_settings.config -resume
    ```

3. Which steps re-ran, and which were reused from the cache? Why?

<details>
<summary><b>Answer</b></summary>

fastp re-runs, because its command changed. Every step that uses the trimmed Illumina reads, directly or
indirectly, also re-runs: Shovill, short-read polishing, typing, AMR, plasmids, phylogeny and MultiQC. Steps that
use only the Nanopore reads (Filtlong, NanoPlot, Autocycler, Medaka) and the AMRFinderPlus database download are
reused from the cache.

</details>

---

## End-of-day checklist

- [ ] `all_plasmids.tsv` built, and you have looked for plasmid clusters shared between STs.
- [ ] Mash distances between the hybrid plasmids computed.
- [ ] The genes around *bla*CTX-M-15 in `CHI11E` listed.
- [ ] Test (stub) pipeline run completed; you can find each result folder in the pre-computed run.

## Further reading

- [Mash](https://mash.readthedocs.io), [Bakta](https://github.com/oschwengers/bakta),
  [Nextflow training](https://training.nextflow.io), [nf-core](https://nf-co.re)

---

[← Day 4](day04.md) · [Course home](README.md) · **Day 5** · [Day 6 →](day06.md)
