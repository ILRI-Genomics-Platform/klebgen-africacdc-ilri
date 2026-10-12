[Course home](README.md) · **Trainer guide**

# Trainer guide

Everything to prepare before Day 1. Allow **two to three weeks**: downloads, environment builds and the
pre-computed pipeline runs each take hours, and every command in the course should be run once on the course
server before participants arrive.

## Contents

1. [Timeline](#1-timeline)
2. [Hardware](#2-hardware)
3. [The shared folder](#3-the-shared-folder)
4. [Training data](#4-training-data)
5. [Databases](#5-databases)
6. [Software for the hands-on sessions](#6-software-for-the-hands-on-sessions)
7. [Pipeline containers: shared, read-only Singularity images](#7-pipeline-containers-shared-read-only-singularity-images)
   (7.4: phylogeny reference genomes)
8. [Pre-computed results](#8-pre-computed-results)
9. [Known gaps and decisions for you](#9-known-gaps-and-decisions-for-you)

---

## 1. Timeline

| When | Task |
| --- | --- |
| 3 weeks before | Server accounts and a course Unix group; create `$SHARED`; set up the software (section 6) |
| 3 weeks before | Prepare and lock the pipeline's Singularity images (section 7) |
| 3 weeks before | Download data (section 4) and databases (section 5) |
| 2 weeks before | Run the pipeline on the training and capstone sets (section 8) |
| 2 weeks before | Pre-compute the manual steps (section 8) |
| 1 week before | Dry run: one trainer works through every page, Day 1 to Day 10, on a participant account |
| 1 day before | Check every participant can log in and run `use_tools qc && fastp --version`, and has installed their R packages |
| 1 day before | As a participant: `bin/prepare_shared_containers.sh --check` shows all images present |

## 2. Hardware

| Item | Minimum | Comfortable for 25 participants |
| --- | --- | --- |
| CPU cores | 32 | 64–100 (about 4 per participant) |
| RAM | 128 GB | 256 GB or more |
| Disk | 500 GB | 1–2 TB (reads, databases, everyone's results) |

The Kraken2 standard-8 database needs about 8 GB of RAM per running job. If 25 people run Kraken2 at once, stagger
the session or have half the room use pre-computed reports.

## 3. The shared folder

Create a read-only folder and set `SHARED` to it in the Day 1 setup block:

```text
$SHARED/                                   e.g. /var/scratch/global/jjuma/neonatal-amr/shared
├── neonatal-amr-nf/                       the pipeline (this repository)
├── containers/                            read-only Singularity images: pipeline (section 7) and course tools (section 6)
│   ├── nf-core/                           images for the nf-core modules
│   ├── local/                             images built for the course's local modules
│   └── MANIFEST.tsv                       source, size and checksum of every image
├── references/                            read-only phylogeny reference genomes (section 7.4)
├── samplesheet.csv                        training set (29 isolates), absolute paths to reads
├── data/                                  downloaded reads (training set)
├── training_data/                         metadata_compact.tsv, ast_compact.tsv, published_K_O_all.tsv, ...
├── databases/
│   ├── kraken2/                           Kraken2 standard-8 (with Bracken files)
│   ├── amrfinderplus/latest/              AMRFinderPlus database, plus amrfinderdb.tar.gz
│   ├── checkm2/uniref100.KO.1.dmnd
│   └── bakta/db-light/
├── precomputed/                           section 7
├── capstone/                              samplesheet.csv, metadata_capstone.tsv, ast_capstone.tsv, precomputed/
└── public_genomes/ST39/                   optional, Day 9
```

```bash
export SHARED=/var/scratch/global/jjuma/neonatal-amr/shared
mkdir -p $SHARED && cp -r neonatal-amr-nf $SHARED/
cp -r $SHARED/neonatal-amr-nf/training_data $SHARED/
```

`$SHARED` sits under your own scratch folder, so participants must be allowed to pass through the folders above it,
and to read (not write) everything inside it. Do this after copying the data and databases:

```bash
chmod o+x /var/scratch/global/$USER /var/scratch/global/$USER/neonatal-amr   # pass through only, no listing
chmod -R o+rX,go-w $SHARED                                                   # read-only for participants
ls -ld /var/scratch/global/$USER $SHARED                                     # check: ...x and r-x for others
```

The container folder is locked separately by `bin/prepare_shared_containers.sh` (section 7), which also checks the
folders above it. If your scratch folder is reset by HPC policy (`chmod 700`) or purged, run these again.

> [!NOTE]
> `$SHARED` uses your username written out (`jjuma`) in `conf/*.config` and `course/day01.md`, not `$USER`: for
> participants, `$USER` would be their own name. A different trainer changes those three files and this section.

## 4. Training data

The training sets are already in `training_data/`. To regenerate them from the study's supplementary spreadsheet
(`13073_2024_1342_MOESM2_ESM.xlsx` from [the paper](https://doi.org/10.1186/s13073-024-01342-3)):

```bash
cd $SHARED/neonatal-amr-nf
for size in compact capstone full; do
    python3 bin/select_training_set.py 13073_2024_1342_MOESM2_ESM.xlsx training_data --size $size
done
```

| Set | Isolates | Used for |
| --- | --- | --- |
| `compact` | 29 | Days 1–9 |
| 5 test isolates | `CHI11E`, `BHB7VE`, `CAA65Z`, `CAA99X`, `D29438` (from `compact`; the pipeline's test set) | Day 10 capstone |
| `capstone` | 36 (full minus compact) | Optional: a larger, unseen Day 10 data set |
| `full` | 65 | Optional larger course, or trainers' reference run |

Download the reads from ENA and write the sample sheets (check with `--dry-run` first):

```bash
cd $SHARED
python3 neonatal-amr-nf/bin/fetch_ena.py fetch neonatal-amr-nf/training_data/accessions_training_compact.tsv \
    -d data -o samplesheet.csv --dry-run
python3 neonatal-amr-nf/bin/fetch_ena.py fetch neonatal-amr-nf/training_data/accessions_training_compact.tsv \
    -d data -o samplesheet.csv --threads 4

# Day 10 capstone: the 5 test isolates, taken from the training set (no extra download)
mkdir -p capstone
TEST="CHI11E BHB7VE CAA65Z CAA99X D29438"
awk -F',' -v ids="$TEST" 'BEGIN { n = split(ids, a, " "); for (i = 1; i <= n; i++) keep[a[i]] = 1 }
                          NR == 1 || ($1 in keep)' samplesheet.csv > capstone/samplesheet.csv
for t in metadata ast; do
    awk -F'\t' -v ids="$TEST" 'BEGIN { n = split(ids, a, " "); for (i = 1; i <= n; i++) keep[a[i]] = 1 }
                               NR == 1 || ($1 in keep)' neonatal-amr-nf/training_data/${t}_compact.tsv > capstone/${t}_capstone.tsv
done
wc -l capstone/*.csv capstone/*.tsv          # 6 lines each: header + 5 isolates
```

For a larger, unseen capstone instead (36 isolates, including an ST14 PICU/HDU cluster in 2019), download the
`capstone` set and use its tables:

```bash
python3 neonatal-amr-nf/bin/fetch_ena.py fetch neonatal-amr-nf/training_data/accessions_training_capstone.tsv \
    -d capstone/data -o capstone/samplesheet.csv --threads 4
cp neonatal-amr-nf/training_data/metadata_capstone.tsv neonatal-amr-nf/training_data/ast_capstone.tsv capstone/
```
Day 10's scenario text describes the 5 test isolates; rewrite it if you use the 36.

Expect about 20–30 GB for the training set. The script checks every file against ENA's MD5 and skips files
already downloaded, so you can re-run it after an interruption.

> [!IMPORTANT]
> The Day 1 setup links reads by reading the sample sheet's columns in order
> (`sample,species,fastq_1,fastq_2,long_fastq,phylo_group`). Keep that column order, and keep absolute paths.

## 5. Databases

| Database | Used on | How to get it |
| --- | --- | --- |
| Kraken2 standard-8 (with Bracken files) | Day 1 | Download the latest `k2_standard_08gb` from the [Kraken2 index page](https://benlangmead.github.io/aws-indexes/k2) and unpack it into `databases/kraken2/` |
| AMRFinderPlus | Days 3–5 | `amrfinder_update -d $SHARED/databases/amrfinderplus`, then `tar czf $SHARED/databases/amrfinderplus/amrfinderdb.tar.gz -C $SHARED/databases/amrfinderplus/latest .` |
| CheckM2 | Day 2 | `checkm2 database --download --path $SHARED/databases/checkm2` |
| Bakta (light) | Day 5 | `bakta_db download --output $SHARED/databases/bakta --type light` |

> [!NOTE]
> Record each database's version or download date and give it to participants. AMRFinderPlus results in
> particular change between database releases.

## 6. Software for the hands-on sessions

Participants install nothing. In the hands-on sessions every tool runs from a **Singularity image** in the
shared, read-only image folder (`$SHARED/containers`): the same images the pipeline uses (section 7), so
manual results and pipeline results come from the same versions. `course/tools.sh` does this for them. Day 1
adds it to every participant's `~/.bashrc`, and each session switches a group of tools on with
`use_tools <group>` (in place of `conda activate`):

| Group | Tools | Image |
| --- | --- | --- |
| `qc` | seqkit, FastQC, fastp, NanoPlot, Filtlong, Kraken2, Bracken, MultiQC | pipeline images; seqkit and FastQC: course images |
| `assembly` | Shovill, Autocycler (with Flye, Raven), BWA, Polypolish, Pypolca, QUAST, CheckM2, Bandage | pipeline images; CheckM2, Bandage: course images |
| `medaka` | Medaka | pipeline image |
| `typing` | Kleborate (with Kaptive), ECTyper, EzClermont | pipeline images |
| `amr` | AMRFinderPlus | pipeline image |
| `plasmids` | MOB-suite (with Mash), Bakta | pipeline image; Bakta: course image |
| `phylo` | Snippy (with snp-sites), Gubbins, snp-dists, IQ-TREE, seqkit | pipeline images; seqkit: course image |
| `r` | R | module `R/4.5`; packages in each participant's own library |
| `nextflow` | Nextflow | module `nextflow/25.04` |

**Course-only images.** Five tools are used in the sessions but not by the pipeline: seqkit 2.8.2, FastQC
0.12.1, CheckM2 1.0.2, Bandage 0.8.1 and Bakta 1.9.4 (`course/envs/<tool>/environment.yml`).
`bin/prepare_shared_containers.sh` (section 7) builds them into `$SHARED/containers/course/` with the pipeline's
images, and lists them in `MANIFEST.tsv`, where `use_tools` finds them. Check that everything is in place:

```bash
bin/prepare_shared_containers.sh --check          # every image, pipeline and course: present / MISSING
source course/tools.sh && list_tools              # as a participant: every tool and its image
```

Keep the CheckM2 and Bakta databases in step with these versions: CheckM2 1.0.2 reads
`uniref100.KO.1.dmnd`; Bakta 1.9.4 needs a database of version 5.x (`db-light` is enough).

**Modules instead of images.** If your HPC provides the tools as modules, participants set
`export TOOL_MODE=module` in `~/.bashrc` (Day 1, section 0.1) and `use_tools <group>` runs `module load` for
the group. Set the module names for your HPC in `GROUP_MODULES` at the top of `course/tools.sh`, keeping the
versions in the table above.

**R.** `use_tools r` loads the R module (`R_MODULE`, default `R/4.5`). Each participant installs the R
packages once into their own library (`~/R/<platform>-library/4.5`) on Day 1, on the login node (it needs
internet; 20–40 minutes):

```bash
$SHARED/neonatal-amr-nf/bin/run_R.sh install_r_packages.R     # tidyverse, ape, ggtree, BactDating, ...
```

Other settings, all at the top of `course/tools.sh` (or set in `~/.bashrc` before it is loaded):
`CONTAINERS` (image folder, default `$SHARED/containers`), `TOOL_BINDS` (folders the images must see,
default `/var/scratch,/export/data,$HOME`), `NEXTFLOW_MODULE`, and `SINGULARITY_MODULE` if `singularity` is
not on the default `PATH`.

## 7. Pipeline containers: shared, read-only Singularity images

On the course server the pipeline runs with `-profile singularity -c <site config>`. Every tool runs inside a
Singularity (or Apptainer) image taken from one shared folder, `$SHARED/containers`, which participants can read
but not change. Nothing is downloaded or built during the course.

| Site config | Use it for | Set in it |
| --- | --- | --- |
| `conf/slurm.config` | ILRI HPC: Slurm executor, `batch`/`highmem` queues, node per step (Kraken2 on compute07, assembly on compute05, alignment and trees on compute06), `/var/scratch` and Kraken2 bind mounts, resource labels | `shared_container_dir`; `resourceLimits` (your largest node); the `module load` in `beforeScript` |
| `conf/local.config` | One shared server, no scheduler | `shared_container_dir`; `resourceLimits` per run |

Both include `conf/singularity_images.config`, which points Nextflow at the read-only image folder and names the
four course-built images. Participants set `SITE_CONFIG` to the one you choose in their Day 1 setup (it defaults
to `conf/slurm.config`).

| Folder | Images | Where they come from |
| --- | --- | --- |
| `containers/nf-core/` | 17, for the 20 nf-core module processes | Downloaded from the nf-core / BioContainers and Seqera image stores, named as Nextflow expects |
| `containers/local/` | 4: Kleborate 3.2.4; Autocycler + Flye + Raven; Polypolish + BWA; EzClermont | Built from each local module's `environment.yml` |
| `containers/course/` | 5: seqkit, FastQC, CheckM2, Bandage, Bakta | Hands-on sessions only (section 6); built from `course/envs/<tool>/environment.yml` |

### 7.1 Prepare the images (once, as the trainer who will own them)

You need Nextflow, Python 3, `curl`, network access, and Singularity or Apptainer able to build images (root,
`sudo`, or `--fakeroot`).

The image folder is `/var/scratch/global/jjuma/neonatal-amr/shared/containers` by default, in the script and in both site configs
(`shared_container_dir`). It must be a folder you, the trainer, can create and own.

```bash
cd $SHARED/neonatal-amr-nf

# What is needed, and what is already there (changes nothing)
bin/prepare_shared_containers.sh --check

# Download, build, write MANIFEST.tsv, then lock. With --group, only members of that
# Unix group can read the images; without it, every account on the server can.
BUILD_CMD="singularity build --fakeroot" \
    bin/prepare_shared_containers.sh --group amr-course

# A different folder: --dir /path/to/containers  (then set the same path as shared_container_dir in the site config)
```

If the script stops with "does not exist, and you cannot create it", you have no write permission in the parent
folder: check you are using your own scratch folder (`/var/scratch/global/$USER/...`), or pass a folder you can
create with `--dir` and set the same path as `shared_container_dir` in the site configs.

If `--fakeroot` is not enabled on your server, use `BUILD_CMD="sudo singularity build"`. On Apptainer-only
servers, `BUILD_CMD="apptainer build --fakeroot"`; Apptainer normally also installs a `singularity` command,
which is what Nextflow calls.

> [!IMPORTANT]
> The image list comes from the pipeline itself (`nextflow inspect`), so after any change to a module or to
> `conf/singularity_images.config`, run the script again. It only downloads or builds what is missing.
> If you use `--dir`, use the same folder as `shared_container_dir` in your site config.

### 7.2 What "read-only" means here

After locking:

| | Trainer (owner) | Participant (course group) | Anyone else |
| --- | --- | --- | --- |
| With `--group amr-course` | read; re-run the script to update | read only | no access |
| Without `--group` | read; re-run the script to update | read only | read only |

- Folders are `0550` (with a group) or `0555`; image files `0440` or `0444`. **Nobody has write permission**,
  including the owner, so an image cannot be deleted or overwritten by accident.
- The pipeline is configured never to download: `singularity.libraryDir` and `singularity.cacheDir` both point
  at the read-only folder. Before a run starts, the pipeline checks that the four course-built images exist; a
  missing image of either kind stops the run with a message telling the participant to ask a trainer.
- Participants can check the folder themselves with `--check`; any other use of the script stops with
  "Only its owner (a trainer) can add or update images".

> [!TIP]
> Put the pipeline folder itself under the same rule: `chmod -R a-w $SHARED/neonatal-amr-nf`. Participants run it,
> and change settings with their own `-c my_settings.config` (Day 5), never by editing the shared copy.

### 7.3 Where the image folder is set

Set `shared_container_dir` once, in the site config you use (`conf/slurm.config` or `conf/local.config`). It can
also be given on the command line with `--shared_container_dir`, which overrides the site config.

### 7.4 Phylogeny reference genomes

By default the pipeline groups isolates by the species and ST that Kleborate reports (`Kpneumoniae_ST39`, ...)
and maps each group to the reference genome listed for its species in `conf/phylo_references.csv`
(*K. pneumoniae*: NTUH-K2044; *K. quasipneumoniae* and *K. variicola*: ATCC 700603; *E. coli*: K-12 MG1655).
Download them once, on a node with internet access, into the folder the site configs name
(`phylo_reference_dir`, `$SHARED/references`):

```bash
cd $SHARED/neonatal-amr-nf
bin/fetch_references.py --check                # what is present / missing (changes nothing)
bin/fetch_references.py --group amr-course     # download from NCBI (ENA as fallback), then lock read-only
```

To use a closer reference for one ST, add a row with that ST (e.g. a closed ST39 genome) to
`conf/phylo_references.csv` and run the script again. Participants who want the course's own behaviour from
Day 6 (groups from the samplesheet, the group's best own assembly as reference) can run the pipeline with
`--phylo_grouping samplesheet`.

Prokka annotation is optional (`--run_prokka`); its image is prepared with the others in section 7.1.

## 8. Pre-computed results

Participants run some slow steps on one or two samples, then copy pre-computed results for the rest. Pages expect
these folders (all of `$SHARED/precomputed` can be locked read-only in the same way as the images):

| Folder | Contents | Used on |
| --- | --- | --- |
| `precomputed/02_assembly/shovill/<sample>/contigs.fa` (and `contigs.gfa`) | Shovill assembly and graph for all 29 | Day 2 |
| `precomputed/02_assembly/autocycler/<sample>/autocycler_out/consensus_assembly.fasta` | 3 hybrid isolates | Day 2 |
| `precomputed/02_assembly/hybrid/<sample>.hybrid.fasta` | 3 polished hybrid assemblies | Day 2 |
| `precomputed/05_plasmids/shortread/<sample>/` | MOB-suite on each Shovill assembly | Day 5 |
| `precomputed/pipeline_results/` | Full pipeline run, training set | Day 5 |
| `capstone/precomputed/` | Full pipeline run, capstone set | Day 10 |

### 8.1 Run the pipeline

```bash
cd $SHARED/precomputed
nextflow run $SHARED/neonatal-amr-nf/main.nf -profile singularity -c $SHARED/neonatal-amr-nf/conf/slurm.config \
    --input $SHARED/samplesheet.csv --outdir pipeline_results \
    --kraken2_db $SHARED/databases/kraken2 \
    --amrfinderplus_db $SHARED/databases/amrfinderplus/amrfinderdb.tar.gz

cd $SHARED/capstone
nextflow run $SHARED/neonatal-amr-nf/main.nf -profile singularity -c $SHARED/neonatal-amr-nf/conf/slurm.config \
    --input samplesheet.csv --outdir precomputed \
    --kraken2_db $SHARED/databases/kraken2 \
    --amrfinderplus_db $SHARED/databases/amrfinderplus/amrfinderdb.tar.gz
```

### 8.2 Arrange the assembly results for the manual steps

```bash
P=$SHARED/precomputed
R=$P/pipeline_results/02_assembly

for f in $R/shovill/*.contigs.fa; do
    s=$(basename $f .contigs.fa)
    mkdir -p $P/02_assembly/shovill/$s && cp $f $P/02_assembly/shovill/$s/contigs.fa
    [ -f $R/shovill/$s.contigs.gfa ] && cp $R/shovill/$s.contigs.gfa $P/02_assembly/shovill/$s/contigs.gfa
done
for f in $R/autocycler/*.autocycler.fasta; do
    s=$(basename $f .autocycler.fasta)
    mkdir -p $P/02_assembly/autocycler/$s/autocycler_out
    cp $f $P/02_assembly/autocycler/$s/autocycler_out/consensus_assembly.fasta
    cp $R/autocycler/$s.autocycler.gfa $P/02_assembly/autocycler/$s/autocycler_out/consensus_assembly.gfa
done
mkdir -p $P/02_assembly/hybrid && cp $R/hybrid/*.hybrid.fasta $P/02_assembly/hybrid/

# MOB-suite on every Shovill assembly (the pipeline runs it on the best assembly instead)
source $SHARED/neonatal-amr-nf/course/tools.sh && use_tools plasmids
for d in $P/02_assembly/shovill/*/; do
    s=$(basename $d)
    mob_recon --infile $d/contigs.fa --outdir $P/05_plasmids/shortread/$s -n 8
done
```

### 8.3 Time needed

Planning estimates, not measurements: run 3 isolates first and scale from Nextflow's execution report.

| Run | 32-core server | 128-core HPC |
| --- | --- | --- |
| Training set (29 isolates) | about 4–5 h | about 2.5–3 h |
| Capstone set (5 test isolates, 1 hybrid) | about 2–3 h (the hybrid assembly sets the pace) | about 1.5–2 h |
| MOB-suite on 29 Shovill assemblies | under 1 h | under 30 min |

## 9. Known gaps and decisions for you

- **No *E. coli* in the training data.** The *E. coli* steps (Day 3, section 3.5) are optional. If you add
  *E. coli* genomes, give them sample IDs starting `EC_` and species `ecoli` in the sample sheet.
- **AST is disc diffusion from the original study**, with a different panel per period. Ceftriaxone is
  represented by the cefpodoxime disc. Tell participants which breakpoint standard you assume.
- **Public ST39 genomes for Day 9** (option B) are not provided; download a set from Pathogenwatch or NCBI if you
  plan to use that option.
- **Dates are years only**, so the Day 7 dating exercise is about judging temporal signal, not producing precise
  dates.
- **Commands depend on tool versions.** Column names in Kleborate, AMRFinderPlus and Gubbins outputs have changed
  between versions. The course uses column-name matching where it can, but run every page once on your server.

---

[Course home](README.md) · **Trainer guide**
