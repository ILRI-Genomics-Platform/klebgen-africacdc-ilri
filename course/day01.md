[Course home](README.md) · **Day 1** · [Day 2 →](day02.md)

# Day 1: Foundations and read quality control

## At a glance

| Session | Time | Topic | Module |
| --- | --- | --- | --- |
| 1 | 09:00–10:30 | Welcome; the command line | 0 |
| 2 | 11:00–12:30 | Software environments, file formats, metadata | 0 |
| 3 | 13:30–15:00 | Illumina read QC and trimming | 1 |
| 4 | 15:30–17:00 | Nanopore read QC; contamination check | 1 |

**You will produce:** a QC table with a pass/fail decision for every sample.

## Learning objectives

By the end of today you can:

- Move around the file system and inspect FASTQ files from the command line.
- Activate a software environment and record a tool's version.
- Explain which metadata fields neonatal AMR surveillance needs, and why.
- Read FastQC, fastp and NanoPlot reports and decide whether a sample is fit for analysis.
- Detect a mixed culture or contamination with Kraken2 and Bracken.

---

## Module 0: Foundations

### 0.0 Log in and start your session

You have a course account, `user1` to `user30`; your trainer gives you yours and its password. Log in to the HPC
from a terminal (on Windows: PowerShell, MobaXterm or PuTTY), replacing `user7` with your account:

```bash
ssh -y user7@hpc.ilri.cgiar.org
```

You arrive on the **login node**, which everyone shares. Do the course's work in an **interactive session** on a
compute node, with 1 CPU. Half of the class works on each of the two course nodes, by account number:

| Your account | Node | Start your session with |
| --- | --- | --- |
| Odd number: `user1`, `user3`, … `user29` | `compute05` | `interactive -w compute05 -c 1` |
| Even number: `user2`, `user4`, … `user30` | `compute06` | `interactive -w compute06 -c 1` |

```bash
interactive -w compute05 -c 1        # user1, user3, user5, ...  (even numbers: compute06)
hostname                             # check: compute05 or compute06
```

From tomorrow, after section 0.1, the command `start_session` does this for you: it picks the node from your
account number. Type `exit` to end the session at the end of the day.

> [!IMPORTANT]
> Only one step in the course runs on the **login node** instead, because it needs the internet: installing
> your R packages (section 0.3). Everything else runs in your session.

**Keep the connection alive.** Many networks drop an SSH connection that has been quiet for a few minutes, and
your session ends with it. Ask your computer to send a small "keepalive" signal every minute:

| Your terminal | How |
| --- | --- |
| macOS, Linux, Windows PowerShell | `ssh -o ServerAliveInterval=60 -y user7@hpc.ilri.cgiar.org`, or once for all logins: add the two lines below to `~/.ssh/config` on your own computer |
| MobaXterm | Settings → Configuration → SSH → tick **SSH keepalive** |
| PuTTY | Connection → **Seconds between keepalives**: 60 (then save the session) |

```text
Host hpc.ilri.cgiar.org
    ServerAliveInterval 60
```

> [!WARNING]
> **If you are disconnected**, your interactive session ends, and so does any command that was running in it.
> Your files and your course settings are safe. To carry on:
>
> ```bash
> ssh -y user7@hpc.ilri.cgiar.org      # log in again
> start_session                        # a new session on your node
> cd $COURSE                           # or the folder you were working in
> ```
>
> You do not need to run `use_tools` again: a course tool you type is found and run from its image for you
> (section 0.3). Then check which outputs exist (`ls`), and run the unfinished step again for your practice
> isolates, or copy the trainers' results for it with `use_precomputed` (section 0.5).
> The long steps by hand (the Shovill and Autocycler assemblies on Day 2) offer pre-computed results. For a
> pipeline run (Days 5 and 10), go back to the same folder and run the same `nextflow run` command with
> `-resume` added: Nextflow reuses every step that had finished and only redoes the rest.

### 0.1 Set up your course folder

> [!IMPORTANT]
> Run this block at the start of the course. It creates your folders and the shell variables used on every
> page. It is safe to run again (for example if your trainer updates it): it first removes any earlier copy of
> the course lines from `~/.bashrc`, so they are only ever there once.

```bash
# Your trainer will tell you the value of SHARED and THREADS
[ -f ~/.bashrc.before_course ] || cp ~/.bashrc ~/.bashrc.before_course     # a copy of your .bashrc, made once

# Remove any earlier copy of the course lines, then add them (once)
sed -i '/>>> neonatal-amr course >>>/,/<<< neonatal-amr course <<</d' ~/.bashrc
cat >> ~/.bashrc <<'EOF'
# >>> neonatal-amr course >>>
export SHARED=/var/scratch/global/jjuma/neonatal-amr/shared
export COURSE=$HOME/neonatal-amr-course
export READS=$COURSE/reads
export RESULTS=$COURSE/results
export THREADS=1          # CPU threads per command: 1 in the course's interactive sessions
export PRACTICE="CHI11E BHB7VE CAA65Z"   # the isolates you run slow steps on yourself (use_precomputed for the rest)
# Site config for running the pipeline (Day 5 and Day 10): your trainer will tell you which
export SITE_CONFIG=$SHARED/neonatal-amr-nf/conf/slurm.config

# Course software (section 0.3): tools from the shared Singularity images.
# Your trainer will tell you if your server uses modules instead (TOOL_MODE=module).
export TOOL_MODE=singularity
source $SHARED/neonatal-amr-nf/course/tools.sh

# start_session: an interactive session with 1 CPU on your course node
# (odd account numbers, user1, user3, ...: compute05; even numbers: compute06)
start_session() {
    local n=${USER//[!0-9]/} node=compute05
    if [ -n "$n" ] && [ $(( 10#$n % 2 )) -eq 0 ]; then node=compute06; fi
    echo "Starting an interactive session on $node with 1 CPU (type 'exit' to end it)"
    interactive -w "$node" -c 1
}

# pick 'regex' file: print column 1 plus every column whose header matches regex
pick() {
    awk -F'\t' -v re="$1" 'NR==1{for(i=2;i<=NF;i++) if($i ~ re) c[++n]=i}
        {line=$1; for(j=1;j<=n;j++) line=line"\t"$c[j]; print line}' "$2"
}
# <<< neonatal-amr course <<<
EOF
source ~/.bashrc

mkdir -p $COURSE $READS $RESULTS
cp -r $SHARED/training_data $COURSE/
ls $COURSE/training_data
```

> [!TIP]
> **Switching the course set-up off** (after the course, or to use your own software): everything above sits
> between the `>>> neonatal-amr course >>>` and `<<< neonatal-amr course <<<` lines of `~/.bashrc`. This
> comments those lines out, and the second command puts them back:
>
> ```bash
> sed -i '/>>> neonatal-amr course >>>/,/<<< neonatal-amr course <<</ s/^/#OFF# /' ~/.bashrc   # off
> sed -i 's/^#OFF# //' ~/.bashrc                                                              # on again
> ```
>
> Then open a new terminal (or run `exec bash`). Your files in `$COURSE` are not touched.

<details>
<summary><b>The course lines are in my <code>~/.bashrc</code> more than once. How do I clean it up?</b></summary>

Count the copies, keep a backup, remove every copy, then run the block above again to add one:

```bash
grep -c ">>> neonatal-amr course >>>" ~/.bashrc       # copies with markers
grep -c "^export SHARED=" ~/.bashrc                    # all copies, including older ones without markers
cp ~/.bashrc ~/.bashrc.backup

# 1. Every copy between the >>> and <<< marker lines
sed -i '/>>> neonatal-amr course >>>/,/<<< neonatal-amr course <<</d' ~/.bashrc
# 2. Copies added before the markers existed: from 'export SHARED=...' to the closing '}' of pick()
sed -i '\#^export SHARED=/var/scratch/global/jjuma/neonatal-amr/shared$#,/^}$/d' ~/.bashrc

grep -c "^export SHARED=" ~/.bashrc                    # now 0; run the block above again, then: exec bash
```

Do not use `sort -u` or `uniq` on `~/.bashrc`: they remove repeated lines anywhere in the file, such as the `}`
that closes each function, and break it. If anything goes wrong, `cp ~/.bashrc.backup ~/.bashrc` puts it back.

</details>

Now make short, readable names for the read files. The sample sheet lists each sample's files:

```bash
head -3 $SHARED/samplesheet.csv

# Link each sample's reads as <sample>_R1/_R2/_ont.fastq.gz
tail -n +2 $SHARED/samplesheet.csv | while IFS=, read -r sample species r1 r2 ont group; do
    ln -sf "$r1" $READS/${sample}_R1.fastq.gz
    ln -sf "$r2" $READS/${sample}_R2.fastq.gz
    [ -n "$ont" ] && ln -sf "$ont" $READS/${sample}_ont.fastq.gz
done
ls $READS | head
ls $READS | wc -l
```

<details>
<summary><b>How many files should <code>$READS</code> contain?</b></summary>

29 samples × 2 Illumina files = 58, plus 3 Nanopore files for the hybrid isolates (`CHI11E`, `CNS1F1`,
`1013757`) = **61 files**.

</details>

### 0.2 The command line in ten commands

| Command | What it does | Example |
| --- | --- | --- |
| `pwd` | Show where you are | `pwd` |
| `ls -lh` | List files with sizes | `ls -lh $READS` |
| `cd` | Change folder | `cd $RESULTS` |
| `mkdir -p` | Make a folder (and its parents) | `mkdir -p $RESULTS/qc` |
| `less` | Page through a file (`q` to quit) | `less $COURSE/training_data/metadata_compact.tsv` |
| `head` / `tail` | First / last lines | `head -4 file.tsv` |
| `zcat` | Read a gzipped file | `zcat $READS/CAA65Z_R1.fastq.gz \| head -8` |
| `wc -l` | Count lines | `wc -l file.tsv` |
| `grep` | Find lines matching a pattern | `grep ST39 $COURSE/training_data/metadata_compact.tsv` |
| `cut -f` | Pick columns from a tab-separated file | `cut -f1,2,3 file.tsv` |

**Try it:** look at the first read of one sample.

```bash
zcat $READS/CAA65Z_R1.fastq.gz | head -4
```

A FASTQ record has four lines: a header starting with `@`, the sequence, a `+` separator, and one quality
character per base.

<details>
<summary><b>How many reads are in <code>CAA65Z_R1.fastq.gz</code>?</b></summary>

Each read is 4 lines, so divide the line count by 4:

```bash
echo $(( $(zcat $READS/CAA65Z_R1.fastq.gz | wc -l) / 4 ))
```

`seqkit stats` (next section) gives the same number faster.

</details>

### 0.3 Software

Nothing is installed in your account. The tools come from **Singularity images** that your trainers prepared in
`$SHARED/containers/`: one file per tool, holding the tool and everything it needs, at a fixed version. These
are the same images the pipeline uses on Day 5, so your results by hand and the pipeline's results come from
the same software. (On some servers your trainer uses **modules** instead; the commands are the same.)

Before using a group of tools, switch it on with `use_tools`:

| Group | Tools | Used on |
| --- | --- | --- |
| `qc` | seqkit, FastQC, fastp, NanoPlot, Filtlong, Kraken2, Bracken, MultiQC | Day 1 |
| `assembly` | Shovill, Autocycler, Flye, Raven, BWA, Polypolish, Pypolca, QUAST, CheckM2, Bandage | Day 2 |
| `medaka` | Medaka | Day 2 |
| `typing` | Kleborate, Kaptive, ECTyper, EzClermont | Days 3, 8 |
| `amr` | AMRFinderPlus | Days 3–4 |
| `plasmids` | MOB-suite, Mash, Bakta | Days 4–5 |
| `phylo` | Snippy, SKA2, Gubbins, snp-sites, snp-dists, IQ-TREE, seqkit | Days 6–7 |
| `r` | R with tidyverse, ggtree, BactDating (from your own R library) | Days 4, 7–9 |
| `nextflow` | Nextflow | Days 5, 10 |

```bash
use_tools qc
fastp --version
type fastp          # 'fastp is a function': it runs fastp inside its image for you
list_tools          # every tool and the image (or module) it comes from
```

After `use_tools qc` you type `fastp ...` as usual; behind the scenes it runs
`singularity exec $SHARED/containers/<fastp image> fastp ...`. Switching on another group keeps the earlier
groups' tools available.

**Forgot `use_tools`, or logged in again after a disconnect?** The course tools still work: when you type a
course tool that is not switched on, it is run from its image for you (in module mode, its modules are loaded
first); `nextflow`, `R` and `Rscript` load their modules the same way. `use_tools` remains the clear way to start each step, because it checks that every image is in place,
and `use_tools all` switches on every group at once. `type fastp` only says "function" after `use_tools`.

```bash
seqkit version            # works in a new session, without use_tools
export TOOLS_AUTO=0       # (optional) switch the automatic loading off for this session
```

**Once, today: install your R packages.** The R sessions (Days 4 and 7–9) use packages from your own R library
in your home folder. Installing them needs internet, so run this on the **login node**; it takes 20–40
minutes, so start it before a break:

```bash
$SHARED/neonatal-amr-nf/bin/run_R.sh install_r_packages.R      # ends with "All packages installed."
```

> [!TIP]
> Record the version of every tool you use. A result from Kleborate 2 and one from Kleborate 3 are not always
> comparable, and AMR databases change every few weeks.

### 0.4 Metadata

Open the course metadata:

```bash
less -S $COURSE/training_data/metadata_compact.tsv
```

| Field | Why it matters |
| --- | --- |
| Sample ID | Joins every result table together; must be unique |
| Collection date | Needed to tell transmission from separate introductions, and to date trees |
| Facility and ward | Places isolates in space; clusters within one ward suggest transmission |
| Patient age (days) | Early-onset (≤72 h) vs late-onset sepsis point to different sources |
| Specimen (blood, CSF) | Invasive disease vs colonisation |
| AST method and breakpoints | CLSI and EUCAST breakpoints differ; needed to compare genotype with phenotype |

> [!WARNING]
> Never put patient names, hospital numbers or full dates of birth in a sequencing metadata sheet. Use coded IDs
> and keep the link file separately, under your data-sharing agreement.

### Check your understanding

<details>
<summary><b>1. The metadata gives collection year only, not full dates. What analysis does that limit?</b></summary>

Dating a cluster (Day 7). With years only, a time-scaled tree can say a cluster arose in roughly a given year,
not in which month. For transmission you also cannot tell whether two infections were days or months apart.

</details>

<details>
<summary><b>2. Two isolates from the same ward have the same sample ID. What goes wrong?</b></summary>

One result file overwrites the other, or the join in R silently duplicates rows. The pipeline you use on Day 5
stops with an error if a sample ID is repeated.

</details>

### 0.5 Practice isolates and pre-computed results

Some steps take minutes per sample, and with 1 CPU each, running them on all 29 isolates would take most of a
day. So for each slow step you **run the command yourself on three practice isolates**, `$PRACTICE`
(`CHI11E`, `BHB7VE`, `CAA65Z`: three ST39 isolates; `CHI11E` also has Nanopore reads), then **copy the
trainers' results for the other isolates**. The trainers ran exactly the same commands on every sample, and their
results sit in the same folders as yours, so the next step finds everything where it expects it:

```bash
echo $PRACTICE                    # CHI11E BHB7VE CAA65Z
use_precomputed 01_qc/kraken2     # copies $SHARED/precomputed/01_qc/kraken2 into $RESULTS/01_qc/kraken2
use_precomputed                   # lists what is available
```

`use_precomputed` never overwrites a file you made yourself, so it is safe to run at any time, including after a
disconnect. Where a tool writes one table for all samples (Kleborate, CheckM2, Kaptive), you write your
practice run to a separate `..._practice` file, compare your rows with the trainers' table, and then use the
trainers' table for the rest of the exercise.

> [!TIP]
> Time to spare? Run a step on more isolates: put their names in place of `$PRACTICE`, e.g.
> `for s in $PRACTICE CES13U D38119; do ...`.

---

## Module 1: Read quality control

### 1.1 Read statistics

```bash
use_tools qc
mkdir -p $RESULTS/01_qc
# Your practice isolates
seqkit stats -a -T -j $THREADS $(for s in $PRACTICE; do ls $READS/${s}_*.fastq.gz; done) \
    > $RESULTS/01_qc/read_stats_practice.tsv
cut -f1,4,5,7,13 $RESULTS/01_qc/read_stats_practice.tsv

# The same table for every sample, from the trainers: check that your rows match it
use_precomputed 01_qc/read_stats.tsv
grep -E "$(echo $PRACTICE | tr ' ' '|')" $RESULTS/01_qc/read_stats.tsv | cut -f1,4,5,7,13
```

Estimate **sequencing depth** for each Illumina sample: total bases in R1 and R2 divided by genome size
(about 5.5 Mb for *K. pneumoniae*).

```bash
awk -F'\t' 'NR>1 && $1 ~ /_R[12]\.fastq\.gz$/ {
    s=$1; sub(/.*\//,"",s); sub(/_R[12]\.fastq\.gz$/,"",s); bases[s]+=$5
} END { for (s in bases) printf "%s\t%.0fx\n", s, bases[s]/5.5e6 }' \
    $RESULTS/01_qc/read_stats.tsv | sort > $RESULTS/01_qc/depth.tsv
cat $RESULTS/01_qc/depth.tsv
```

> [!NOTE]
> Working targets: **about 40× or more** Illumina depth for typing and SNP calling; **about 50× or more** Nanopore
> depth, with read N50 above 10 kb, for complete plasmids. These are guides, not hard rules.

### 1.2 Illumina: FastQC, then trimming with fastp

```bash
mkdir -p $RESULTS/01_qc/fastqc $RESULTS/01_qc/fastp $COURSE/trimmed

# FastQC on two samples first, to learn the report
fastqc -t $THREADS -o $RESULTS/01_qc/fastqc $READS/CAA65Z_R{1,2}.fastq.gz $READS/CES13U_R{1,2}.fastq.gz

# Trim your practice isolates
for s in $PRACTICE; do
    fastp -i $READS/${s}_R1.fastq.gz -I $READS/${s}_R2.fastq.gz \
          -o $COURSE/trimmed/${s}_R1.fastq.gz -O $COURSE/trimmed/${s}_R2.fastq.gz \
          --qualified_quality_phred 20 --length_required 50 \
          -j $RESULTS/01_qc/fastp/${s}.fastp.json -h $RESULTS/01_qc/fastp/${s}.fastp.html \
          -w $THREADS 2>&1 | tee $RESULTS/01_qc/fastp/${s}.log
done

# fastp reports for the other samples, from the trainers (their trimmed reads are not needed later)
use_precomputed 01_qc/fastp

# One report for all samples
multiqc $RESULTS/01_qc -o $RESULTS/01_qc/multiqc
```

Open `$RESULTS/01_qc/multiqc/multiqc_report.html` in your browser and look at:

1. **Per-base quality**: does quality fall towards the end of the reads?
2. **Adapter content**: are adapters present before trimming? After?
3. **Duplication**: very high duplication suggests a low-input library.
4. **Reads passing filters** in the fastp section.

### 1.3 Nanopore: NanoPlot and Filtlong

Only the three hybrid isolates have Nanopore reads. You process `CHI11E` (the hybrid isolate in your practice
set, assembled on Day 2); the trainers' reports cover `CNS1F1` and `1013757`.

```bash
mkdir -p $RESULTS/01_qc/nanoplot
for s in CHI11E; do
    NanoPlot --fastq $READS/${s}_ont.fastq.gz -o $RESULTS/01_qc/nanoplot/${s}_raw -t $THREADS

    filtlong --min_length 1000 --keep_percent 90 $READS/${s}_ont.fastq.gz \
        | gzip > $COURSE/trimmed/${s}_ont.fastq.gz

    NanoPlot --fastq $COURSE/trimmed/${s}_ont.fastq.gz -o $RESULTS/01_qc/nanoplot/${s}_filtered -t $THREADS
done
use_precomputed 01_qc/nanoplot          # CNS1F1 and 1013757
grep -E "Read length N50|Mean read quality|Total bases" $RESULTS/01_qc/nanoplot/*/NanoStats.txt
```

<details>
<summary><b>Filtlong removed 10% of the reads. Which 10%, and why does that help assembly?</b></summary>

`--keep_percent 90` drops the worst 10% of bases, scored on read length and quality, after first removing reads
under 1 kb. Long, accurate reads are what let an assembler span repeats and close plasmids; short, poor reads
mostly add noise.

</details>

### 1.4 Contamination and species check: Kraken2 and Bracken

```bash
mkdir -p $RESULTS/01_qc/kraken2
KDB=$SHARED/databases/kraken2          # Kraken2 Standard-8, with Bracken files (needs ~8 GB of memory)
cat $KDB/VERSION                       # record the database version with your results

for s in $PRACTICE; do                 # each takes 5-15 minutes with 1 thread
    kraken2 --db $KDB --threads $THREADS --paired \
        --report $RESULTS/01_qc/kraken2/${s}.k2report --output /dev/null \
        $COURSE/trimmed/${s}_R1.fastq.gz $COURSE/trimmed/${s}_R2.fastq.gz
    bracken -d $KDB -i $RESULTS/01_qc/kraken2/${s}.k2report \
        -o $RESULTS/01_qc/kraken2/${s}.bracken.tsv -l S -r 150
done
use_precomputed 01_qc/kraken2          # every other sample, including the two QC failures

# Top three species per sample
for f in $RESULTS/01_qc/kraken2/*.bracken.tsv; do
    echo "== $(basename $f .bracken.tsv)"
    sort -t$'\t' -k7,7gr $f | head -3 | cut -f1,7
done
```

> [!WARNING]
> Kraken2 often splits reads between *K. pneumoniae*, *K. quasipneumoniae* and *K. variicola*, because these
> species are close relatives. A few percent of "other *Klebsiella*" is normal. Two different genera, each with a
> large share, is not. Confirm the species on Day 3 with Kleborate.

### 1.5 Make your QC decision

Copy this table into a file (`$RESULTS/01_qc/qc_decisions.tsv`) and fill it in for every sample:

| Sample | Depth (×) | Reads passing fastp (%) | Top species (%) | Second species (%) | Decision |
| --- | --- | --- | --- | --- | --- |
| CAA65Z | | | | | pass / re-sequence / re-culture |

Use these rules:

- **Pass**: depth ≥ 40× and the top species ≥ 90%.
- **Re-sequence**: depth too low, but the culture looks pure.
- **Re-culture**: a second species with a large share (a mixed culture).

### Check your understanding

<details>
<summary><b>1. Two samples were included as QC failures: <code>D38119</code> and <code>CES13U</code>. What did the study find for them?</b></summary>

From the published tables (`training_data/training_set_compact.tsv`):

- `D38119`: assembly of **9.46 Mb** in 141 contigs. That is almost twice the expected *K. pneumoniae* genome size,
  which suggests two genomes in one culture. Check whether Kraken2 shows a second organism.
- `CES13U`: **1,468 contigs** (expected: under about 100). A very fragmented assembly usually means low depth,
  short or poor reads, or contamination. Check its depth and fastp report.

You will see both again when you assemble on Day 2.

</details>

<details>
<summary><b>2. A sample has 25× depth and is 99% <i>K. pneumoniae</i>. Pass, re-sequence or re-culture?</b></summary>

**Re-sequence** by the rules above: the culture is pure but depth is below 40×. In practice it may still be usable
for MLST and gene detection, but SNP calling for transmission analysis needs good depth across the genome.

</details>

<details>
<summary><b>3. Why trim before Kraken2 rather than after assembly?</b></summary>

Kraken2 classifies reads. Adapters and low-quality ends add k-mers that match nothing, or the wrong thing, and
lower the share classified. Checking reads also catches contamination before you spend hours assembling.

</details>

---

## End-of-day checklist

- [ ] `$COURSE`, `$READS`, `$RESULTS`, `$SHARED` and `$THREADS` are set in a new terminal (`echo $COURSE`).
- [ ] `$READS` has 61 linked read files.
- [ ] `$COURSE/trimmed/` has trimmed reads for all 29 samples and filtered Nanopore reads for 3.
- [ ] You have opened the MultiQC report.
- [ ] `$RESULTS/01_qc/qc_decisions.tsv` has a decision for every sample.

## Further reading

- [fastp](https://github.com/OpenGene/fastp), [NanoPlot](https://github.com/wdecoster/NanoPlot),
  [Filtlong](https://github.com/rrwick/Filtlong), [Kraken2](https://github.com/DerrickWood/kraken2/wiki),
  [Bracken](https://github.com/jenniferlu717/Bracken)

---

[Course home](README.md) · **Day 1** · [Day 2 →](day02.md)
