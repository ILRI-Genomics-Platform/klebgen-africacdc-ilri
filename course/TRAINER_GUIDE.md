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
| 3 weeks before | Participant accounts `user1` to `user30` and a course Unix group (section 2); create `$SHARED`; set up the software (section 6) |
| 3 weeks before | Prepare and lock the pipeline's Singularity images (section 7) |
| 3 weeks before | Download data (section 4) and databases (section 5, after the images of section 7) |
| 2 weeks before | Run the pipeline on the training and capstone sets (section 8) |
| 2 weeks before | Pre-compute the manual steps for all isolates: `course/make_precomputed.sh` (section 8.2) |
| 1 week before | Dry run: one trainer works through every page, Day 1 to Day 10, on a participant account |
| 1 day before | Check every account can log in (`ssh -y userN@hpc.ilri.cgiar.org`), start a session on its node and run `use_tools qc && fastp --version` |
| 1 day before | As a participant: `bin/prepare_shared_containers.sh --check` shows all images present |

## 2. Hardware

| Item | Minimum | Comfortable for 25 participants |
| --- | --- | --- |
| CPU cores | 32 | 64–100 (about 4 per participant) |
| RAM | 128 GB | 256 GB or more |
| Disk | 500 GB | 1–2 TB (reads, databases, everyone's results) |

The Kraken2 standard-8 database needs about 8 GB of RAM per running job. If 25 people run Kraken2 at once, stagger
the session or have half the room use pre-computed reports.

**Participant accounts and nodes (ILRI HPC).** Participants log in as `user1` to `user30`
(`ssh -y userN@hpc.ilri.cgiar.org`) and work in interactive sessions with 1 CPU, split over the two course nodes
by account number:

| Accounts | Node | Command | Participants per node |
| --- | --- | --- | --- |
| Odd: `user1`, `user3`, … `user29` | `compute05` | `interactive -w compute05 -c 1` | 15 |
| Even: `user2`, `user4`, … `user30` | `compute06` | `interactive -w compute06 -c 1` | 15 |

Day 1 (section 0.1) adds a `start_session` command to each participant's `~/.bashrc` that picks the node from
the account number. Before the course, check that:

- `interactive -c 1` gives a session enough memory for Kraken2 Standard-8 (about 10 GB). If its default is lower,
  tell participants the option for more memory, or have them use the pre-computed Kraken2 reports.
- Each node has memory for 15 sessions at once (15 × 10 GB = 150 GB at the Kraken2 step).
- Participant accounts are allowed to request `compute05` and `compute06` with `-w`.
- Idle disconnects. Check whether the login shell logs out idle users (`echo $TMOUT` on the login node and in an
  interactive session; empty means no limit) and how long an idle interactive session is kept. Day 1 (section
  0.0) asks participants to turn on SSH keepalive, which prevents most network drops.

**When a participant is disconnected.** The interactive session ends with the SSH connection, and so does the
command that was running in it; files and `~/.bashrc` settings are kept. After `start_session`, the course tools
load automatically (`command_not_found_handle` in `course/tools.sh` runs a course tool from its image, or loads
its modules, when `use_tools` has not been run), so participants only re-run the interrupted step or copy its
pre-computed results. Nothing short of `screen`/`tmux` or a batch job keeps a running command alive through a
disconnect, so participants run slow steps only on their practice isolates and copy the rest with
`use_precomputed` (section 8).

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
│   ├── kraken2/                           Kraken2 Standard-8 (with Bracken files), VERSION (section 5)
│   ├── amrfinderplus/latest/              AMRFinderPlus database, plus amrfinderdb.tar.gz and VERSION
│   ├── checkm2/CheckM2_database/uniref100.KO.1.dmnd
│   └── bakta/db-light/                    Bakta light database v6
├── precomputed/                           section 8: mirrors a participant's $RESULTS
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
(`assets/13073_2024_1342_MOESM2_ESM.xlsx`, from [the paper](https://doi.org/10.1186/s13073-024-01342-3)), use the
`training-data` conda environment (Python with pandas and openpyxl; `bin/fetch_ena.py` below runs in it too):

```bash
cd $SHARED/neonatal-amr-nf
conda env create -f envs/training_data.yml        # once
conda activate training-data
for size in compact capstone full; do
    python3 bin/select_training_set.py assets/13073_2024_1342_MOESM2_ESM.xlsx training_data --size $size
done
```

The selection is deterministic: the same spreadsheet always gives the same isolates and tables.

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

| Database | Version | Download | On disk | Used on | Where |
| --- | --- | --- | --- | --- | --- |
| Kraken2 Standard-8, with Bracken files | 2026-06-26 (newest when written) | 5.5 GB | 7.5 GB; needs ~8 GB RAM | Day 1; pipeline (Days 5, 10) | `databases/kraken2/` |
| AMRFinderPlus | newest at download | small | small | Days 3–5; pipeline | `databases/amrfinderplus/latest/` and `amrfinderdb.tar.gz` |
| CheckM2 (DIAMOND) | database v2, for CheckM2 1.0.2 | 1.7 GB | a few GB | Day 2 | `databases/checkm2/CheckM2_database/uniref100.KO.1.dmnd` |
| Bakta light | database v6, for Bakta 1.12 | 1.3 GB | 3.9 GB | Day 5; pipeline `--run_bakta` | `databases/bakta/db-light/` |

In all, about **9 GB to download and 15 GB on disk**. Two other databases need nothing from you:

- **MOB-suite** (Days 4–5): its Bioconda package downloads its databases (about 0.5 GB) when it is installed, so
  they are inside its image. Check once: `use_tools plasmids && mob_typer --infile <any assembly FASTA> --out_file
  /tmp/mob_check.txt` should finish without trying to download anything.
- **Kleborate and Kaptive** (Days 3, 8): they ship with their own databases.

Standard-8 is the Kraken2 database for the course: it needs about 8 GB of memory, so it runs in a participant's
session. The full Standard database (80 GB download, 103 GB on disk, about 100 GB of memory) classifies a few
more reads but cannot run there. If your HPC already has one (on the ILRI HPC,
`/export/data/bio/kraken2/db/standard_20241228`), the pipeline runs that the trainers submit can use it with
`--kraken2_db`.

### 5.1 Download (once, as the trainer, on a node with internet access)

The tools come from the same images as the course (`course/tools.sh`), so each database matches the version of
the tool that reads it. **Prepare the images first** (section 7.1, `bin/prepare_shared_containers.sh`): it also
builds the course-only images (CheckM2, seqkit, FastQC, Bandage) that block 3 needs. Run the blocks one at a time;
each ends with a check.

```bash
export SHARED=/var/scratch/global/jjuma/neonatal-amr/shared
source $SHARED/neonatal-amr-nf/course/tools.sh
DB=$SHARED/databases
mkdir -p $DB

# 1. Kraken2 Standard-8 with Bracken files (5.5 GB download, 7.5 GB unpacked).
#    Use the newest file name from https://benlangmead.github.io/aws-indexes/k2
K2=k2_standard_08_GB_20260626
mkdir -p $DB/kraken2 && cd $DB/kraken2
wget -c https://genome-idx.s3.amazonaws.com/kraken/$K2.tar.gz
tar xzf $K2.tar.gz && rm $K2.tar.gz
echo "$K2" > VERSION
ls hash.k2d opts.k2d taxo.k2d database150mers.kmer_distrib     # all four must be there

# 2. AMRFinderPlus, with the pipeline's AMRFinderPlus 4.2.7; plus a .tar.gz copy for the pipeline (--amrfinderplus_db)
use_tools amr
amrfinder_update -d $DB/amrfinderplus
tar czf $DB/amrfinderplus/amrfinderdb.tar.gz -C $DB/amrfinderplus/latest .
amrfinder -d $DB/amrfinderplus/latest --database_version 2>&1 | tee $DB/amrfinderplus/VERSION

# 3. CheckM2 database v2, the one CheckM2 1.0.2 uses (1.7 GB download). Downloaded directly: inside the read-only
#    image, `checkm2 database --download` fails at its last step (it tries to save the path in its own install).
#    The course passes --database_path, so that saved path is not needed.
mkdir -p $DB/checkm2 && cd $DB/checkm2
wget -c -O checkm2_database.tar.gz "https://zenodo.org/records/5571251/files/checkm2_database.tar.gz?download=1"
echo "f35c40f58efaf112d29fc187a88af6f5  checkm2_database.tar.gz" | md5sum -c -
tar xzf checkm2_database.tar.gz && rm checkm2_database.tar.gz
ls -lh CheckM2_database/uniref100.KO.1.dmnd
echo "aff3badcd2bd389539e833b5a0d446c13fb257d71abcb6c73fac38e8ab070acb  CheckM2_database/uniref100.KO.1.dmnd" | sha256sum -c -
use_tools assembly
checkm2 testrun --database_path $DB/checkm2/CheckM2_database/uniref100.KO.1.dmnd --threads 1   # ends with test results

# 4. Bakta light database v6, with the pipeline's Bakta 1.12 (1.3 GB download, 3.9 GB unpacked);
#    bakta_db also sets up the AMRFinderPlus files Bakta needs inside it
use_tools plasmids
bakta_db download --output $DB/bakta --type light
ls $DB/bakta/db-light/version.json && cat $DB/bakta/db-light/version.json

# 5. Make the databases read-only for everyone
find $DB -type d -exec chmod 0555 {} +
find $DB -type f -exec chmod 0444 {} +
du -sh $DB/*
```

> [!NOTE]
> Record each database's version or download date (the `VERSION` files, Bakta's `version.json`) and give
> them to participants. AMRFinderPlus results in particular change between database releases. To update a
> database later, give yourself write permission again (`chmod -R u+w $DB/<name>`) and repeat its block.

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
| `plasmids` | MOB-suite (with Mash), Bakta | pipeline images (Bakta: the pipeline's `--run_bakta` image) |
| `phylo` | Snippy (with snp-sites), SKA2, Gubbins, snp-dists, IQ-TREE, seqkit | pipeline images (SKA2: the `--phylo_method ska2` image); seqkit: course image |
| `r` | R | module `R/4.5`; packages in each participant's own library |
| `nextflow` | Nextflow | module `nextflow/25.04` |

**Course-only images.** Four tools are used in the sessions but not by the pipeline: seqkit 2.8.2, FastQC
0.12.1, CheckM2 1.0.2 and Bandage 0.8.1 (`course/envs/<tool>/environment.yml`).
`bin/prepare_shared_containers.sh` (section 7) builds them into `$SHARED/containers/course/` with the pipeline's
images, and lists them in `MANIFEST.tsv`, where `use_tools` finds them. Check that everything is in place:

```bash
bin/prepare_shared_containers.sh --check          # every image, pipeline and course: present / MISSING
source course/tools.sh && list_tools              # as a participant: every tool and its image
```

Keep the databases in step with these versions: CheckM2 1.0.2 reads database v2; Bakta 1.12 (the pipeline's)
needs a version 6 database. Section 5.1 downloads both with these images, so they match.

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
course-built images (`containers/local/`). Participants set `SITE_CONFIG` to the one you choose in their Day 1 setup (it defaults
to `conf/slurm.config`).

| Folder | Images | Where they come from |
| --- | --- | --- |
| `containers/nf-core/` | 18, for the 21 nf-core module processes (with Bakta) | Downloaded from the nf-core / BioContainers and Seqera image stores, named as Nextflow expects |
| `containers/local/` | 5: Kleborate 3.2.4; Autocycler + Flye + Raven; Polypolish + BWA; EzClermont; SKA2 0.5.1 | Built from each local module's `environment.yml` |
| `containers/course/` | 4: seqkit, FastQC, CheckM2, Bandage | Hands-on sessions only (section 6); built from `course/envs/<tool>/environment.yml` |

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
  at the read-only folder. Before a run starts, the pipeline checks that the course-built images it will use exist (SKA2's only with
  `--phylo_method ska2`); a
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

By default the pipeline puts all isolates of a species in one group, whatever their ST (`Kpneumoniae`, `Ecoli`;
`--phylo_grouping kleborate` makes one group per species + ST instead, e.g. `Kpneumoniae_ST39`),
and maps each group to the reference genome listed for its species in `conf/phylo_references.csv`
(*K. pneumoniae*: NTUH-K2044; *K. quasipneumoniae*: ATCC 700603; *K. variicola*: DSM 15968; *E. coli*: K-12 MG1655).
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

Participants run each slow step themselves on three **practice isolates** (`$PRACTICE`, set in their
`~/.bashrc` on Day 1: `CHI11E BHB7VE CAA65Z`, all ST39; `CHI11E` is also a hybrid isolate), then copy your results
for the other isolates with `use_precomputed <folder>` (in `course/tools.sh`). `$SHARED/precomputed` therefore
mirrors a participant's `$RESULTS`, folder for folder, and `use_precomputed 01_qc/kraken2` copies
`$SHARED/precomputed/01_qc/kraken2` into `$RESULTS/01_qc/kraken2` without overwriting the participant's own files.

| Folder under `$SHARED/precomputed/` | Contents | Participants run themselves | Day |
| --- | --- | --- | --- |
| `01_qc/read_stats.tsv` | seqkit stats, every read file | practice isolates | 1 |
| `01_qc/fastp/` | fastp reports (`.json`, `.html`, `.log`) | practice isolates | 1 |
| `01_qc/nanoplot/` | NanoPlot before and after Filtlong, 3 hybrid isolates | `CHI11E` | 1 |
| `01_qc/kraken2/` | Kraken2 reports and Bracken tables | practice isolates | 1 |
| `02_assembly/shovill/<sample>/` | Shovill `contigs.fa` and `contigs.gfa` | `CHI11E` (more if time) | 2 |
| `02_assembly/checkm2/` | CheckM2 on every Shovill assembly | practice isolates | 2 |
| `02_assembly/autocycler/<sample>/autocycler_out/` | Autocycler consensus, 3 hybrid isolates | `CHI11E` | 2 |
| `02_assembly/hybrid/<sample>.hybrid.fasta` | polished hybrid assemblies | `CHI11E` | 2 |
| `03_typing/kleborate/` | Kleborate on every best assembly | practice isolates | 3 |
| `04_amr/amrfinderplus/<sample>.tsv` | AMRFinderPlus, one file per isolate | practice isolates | 3 |
| `05_plasmids/hybrid/<sample>/`, `05_plasmids/shortread/<sample>/` | MOB-suite on hybrid and Shovill assemblies | `CHI11E` | 4–5 |
| `06_phylogeny/ST39/snippy/<sample>/` | Snippy against NTUH-K2044 (the *K. pneumoniae* reference in `conf/phylo_references.csv`), 12 ST39 isolates | practice isolates | 6 |
| `06_phylogeny/ST39/` | core alignment, Gubbins, IQ-TREE, SNP distances, pair table; SKA2 (`ST39_ska.skf`, `ST39_ska_distances.tsv`) | all 12 (the pre-computed copy is a fallback) | 6 |
| `07_antigens/kaptive_k.tsv`, `kaptive_o.tsv` | Kaptive on every best assembly | practice isolates | 8 |
| `pipeline_results/` | full pipeline run, training set | | 5 |
| `capstone/precomputed/` (beside `precomputed/`) | full pipeline run, capstone set | | 10 |

Steps that are quick with 1 CPU stay on all isolates: QUAST, MultiQC, the plasmid and gene searches on the three
hybrid isolates, and, on Day 6, the core alignment, Gubbins, IQ-TREE and SKA2 for the ST39 group. Day 6 also
has a pre-computed copy of those (`use_precomputed 06_phylogeny/ST39`) for anyone who falls behind.

To change the practice set (for example to fit a shorter course), change `PRACTICE` in the Day 1 set-up block.
Choose ST39 isolates, so that every one can be used in the Day 6 phylogeny, and keep `CHI11E`, the hybrid
isolate used on Day 2. Day 6 maps to the reference that `conf/phylo_references.csv` lists for the isolates'
Kleborate species (NTUH-K2044, AP006725.1, for *K. pneumoniae*), as the pipeline does; `$SHARED/references` must
hold it (section 7.4). If you add an ST39-specific row to the sheet, change the file name on the Day 6 page too.

### 8.1 Run the pipeline

Each run takes hours (section 9 has estimates), so the pipeline runs as a Slurm job: it keeps going after you log
out, and Nextflow submits every step as its own job. `submit_neonatal_amr_reads.sbatch` has a preset for each
run: it sets the samplesheet, the results folder, a work folder in your own scratch space, and the course's
shared databases (Kraken2 Standard-8 and AMRFinderPlus from section 5). Both runs build each phylogeny with
Snippy and with SKA2 (`--phylo_method both`): Days 6, 7 and 10 compare the two with `bin/compare_snp_methods.py`.

```bash
echo $SHARED                                   # must print the shared folder, e.g. /var/scratch/global/jjuma/neonatal-amr/shared
cd $SHARED/neonatal-amr-nf
mkdir -p logs

sbatch --export=ALL,PRESET=training submit_neonatal_amr_reads.sbatch   # 29 isolates -> $SHARED/precomputed/pipeline_results
sbatch --export=ALL,PRESET=capstone submit_neonatal_amr_reads.sbatch   # capstone (section 4) -> $SHARED/capstone/precomputed

squeue -u $USER                                # the two head jobs, then their analysis jobs
tail -f logs/${USER}_reads_*.txt               # progress
```

The script creates the results folders, and stops before starting if the samplesheet, a read file or a
database is missing. If a run stops, submit the same command again: `-resume` reuses every step that already
finished. The work folders (`/var/scratch/global/$USER/neonatal-amr/work_precomputed` and `work_capstone`) can be
deleted once both runs are complete and checked.

### 8.2 Make the participants' pre-computed results

After the training pipeline run (8.1) has finished, `course/make_precomputed.sh` fills `$SHARED/precomputed` with
the folders in the table above. It reuses the pipeline run wherever the pipeline used the same tool, options and
database as the course page, and only renames files or adjusts columns to the course's format:

| From the pipeline run (reused) | Made by the script (the course page's own command) |
| --- | --- |
| fastp reports; Kraken2 reports; Shovill, Autocycler and hybrid assemblies; Kleborate and Kaptive tables (merged into one table each); AMRFinderPlus tables (with the `Name` column the course's `--name` adds); MOB-suite on each best assembly | seqkit stats; NanoPlot and Filtlong (3 hybrid isolates); Bracken (the course uses `-r 150`, the pipeline its defaults); CheckM2; MOB-suite on the Shovill assemblies of the 3 hybrid isolates; fastp and Snippy for the 12 ST39 isolates (same reference as the pipeline, but the pipeline does not keep the per-isolate alignments that `snippy-core` needs), then the Day 6 core alignment, Gubbins, IQ-TREE and SNP distances; SKA2 from the best assemblies (the pipeline's SKA2 used the reads) |

If a pipeline result is missing for a sample, the script runs the course's command for that sample instead.

```bash
cd $SHARED/neonatal-amr-nf && mkdir -p logs
sbatch course/make_precomputed.sh                                         # all steps, 16 CPUs, 64 GB
sbatch --export=ALL,STEPS="kraken2 checkm2" course/make_precomputed.sh    # only some steps
tail -f logs/${USER}_precomputed_*.txt
```

Finished steps are skipped, so after a failure fix the cause and submit it again. Trimmed reads and other
intermediate files go to `/var/scratch/global/$USER/neonatal-amr/precompute_work` (set `WORK` to change it),
not to the shared folder; delete it when you are done. The script ends with a summary of file counts per folder.
Then make the folder read-only for participants:

```bash
chmod -R o+rX,go-w $SHARED/precomputed
```

Before the course, log in as one participant account and check one copy: `use_precomputed 01_qc/kraken2`.

### 8.3 Time needed

Planning estimates, not measurements: run 3 isolates first and scale from Nextflow's execution report.

| Run | 32-core server | 128-core HPC |
| --- | --- | --- |
| Training set (29 isolates) | about 4–5 h | about 2.5–3 h |
| Capstone set (5 test isolates, 1 hybrid) | about 2–3 h (the hybrid assembly sets the pace) | about 1.5–2 h |
| `make_precomputed.sh` after the pipeline run (16 CPUs) | about 2–3 h (Snippy, CheckM2 and Gubbins take most of it) | the same: it runs as one 16-CPU job |

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
