#!/usr/bin/env python3
"""
aggregate_results.py - collect neonatal-amr-nf results into one table per unit of observation.

Most people want "one row per isolate". That works only for results that ARE one row per
isolate, so this script writes several tables instead of forcing everything into one:

  master_table.tsv        one row per ISOLATE   (the wide table: joins + chosen summaries)
  <tool>_all.tsv          one row per TOOL ROW  (--tool-tables; the tool's own columns, verbatim)
  amr_genes_long.tsv      one row per GENE HIT  (--long)
  plasmids_long.tsv       one row per PLASMID   (--long)
  snp_distances_long.tsv  one row per PAIR      (--long)
  aggregate_summary.txt   what was found, what was missing

Run 'aggregate_results.py --describe' for what can and cannot be combined, and why.
Run 'aggregate_results.py --list-sources' for the source names and where they are looked for.

Input can be the pipeline's --outdir, or individual directories, or both: --dir overrides
one source's location, which is how you aggregate results that are not under one outdir
(a precomputed folder, several runs, one tool re-run elsewhere).

Only the standard library is used, so this runs with a bare python3 on the cluster
(unlike bin/compare_snp_methods.py, which needs pandas).

Examples
--------
    bin/aggregate_results.py --results results --long --tool-tables
    bin/aggregate_results.py --dir kleborate=/path/to/kleborate --dir amrfinderplus=/path/to/amr \\
                             --outdir combined --tool-tables
    bin/aggregate_results.py --results results --dir kleborate=$SHARED/precomputed/kleborate
    bin/aggregate_results.py --results results --tool-tables --beside-source
    bin/aggregate_results.py --results results --samplesheet samplesheet.csv --strict
"""

import argparse
import csv
import json
import re
import sys
from collections import Counter, OrderedDict
from pathlib import Path

DESCRIBE = """
What can and cannot be combined
===============================

The question is always: what is one ROW? Two tables can be joined only when they share
a unit of observation (or when one is a lookup for a key the other holds).

1. ONE ROW PER ISOLATE - joined directly into master_table.tsv
   fastp               reads before/after, Q30, duplication
   Kraken2 / Bracken   top species and its percentage
   QUAST               contigs, total length, N50, GC, largest contig
   CheckM2             completeness, contamination
   Kleborate           species, ST, K locus, O type, resistance and virulence scores
   Kaptive             K and O: best match locus, type, confidence, problems
   ECTyper/EzClermont  E. coli serotype and phylogroup
   These are safe: one file (or row) per isolate, so the join cannot change the row count.

2. MANY ROWS PER ISOLATE - can only enter the wide table as a SUMMARY you choose
   AMRFinderPlus       many gene hits per isolate
   MOB-suite           many plasmids per isolate
   Bakta / Prokka      many features per isolate
   Summarising loses information BY DEFINITION. 'n_plasmids = 3' and 'amr_genes = a;b;c'
   cannot tell you WHICH plasmid carries WHICH gene. If you need that, use --long; the
   wide columns are for counting and filtering, not for analysing the relationship.

3. CANNOT BE COMBINED - a different unit of observation entirely
   SNP distances       the unit is a PAIR of isolates, not an isolate. An n-isolate group
                       has n(n-1)/2 pairs. Joining these onto isolates is a category error.
   Trees               the unit is a GROUP, and the content is Newick text, not a table.
   Gubbins predictions the unit is a BRANCH or a recombinant BLOCK.
   Alignments          sequence data.
   Plasmid clusters    the unit is a PLASMID, which can occur in many isolates: a
                       many-to-many (bipartite) relation between plasmids and isolates.
   Per-contig reports  the unit is a CONTIG.
   MultiQC             an HTML report.

4. LEGITIMATE DERIVED COLUMNS - pair or group data reduced to per-isolate facts
   These ARE in master_table.tsv, because the reduction is explicit and stated:
   phylo_group, phylo_group_n          group broadcast onto its members (many-to-one)
   nearest_sample, nearest_snps        from the pairwise matrix, per isolate
   n_within_<threshold>_snps           ditto, with --snp-threshold
   Note the difference: 'the SNP distance' is not an isolate-level fact, but 'the distance
   to this isolate's nearest neighbour' is.

5. PER-TOOL TABLES (--tool-tables) - one tool, every row, its own columns
   Concatenation, not a join: no unit of observation is crossed, nothing is summarised,
   and the tool's columns are kept verbatim with a 'sample' column added. This is the
   lossless form, and the right input for your own analysis in R or pandas.
   kleborate_all.tsv is one row per isolate; amrfinderplus_all.tsv is one row per gene
   hit; mobtyper_all.tsv is one row per plasmid. Same tool, different units - which is
   exactly why they are separate files.

A rule of thumb: if adding a column would change the number of rows, it does not belong
in the wide table. Put it in a long or per-tool table and join when you need it.
"""

# Where each source lives under the pipeline's --outdir. Override any of them with
# --dir NAME=PATH; see --list-sources.
SOURCE_DIRS = OrderedDict((
    ("fastp", "01_qc/fastp"),
    ("kraken2", "01_qc/kraken2"),
    ("bracken", "01_qc/bracken"),
    ("quast", "02_assembly/quast"),
    ("checkm2", "02_assembly/checkm2"),
    ("bakta", "02_assembly/bakta"),
    ("prokka", "02_assembly/prokka"),
    ("kleborate", "03_typing/kleborate"),
    ("kaptive", "03_typing/kaptive"),
    ("ectyper", "03_typing/ectyper"),
    ("ezclermont", "03_typing/ezclermont"),
    ("amrfinderplus", "04_amr/amrfinderplus"),
    ("mobsuite", "05_plasmids/mobsuite"),
    ("phylogeny", "06_phylogeny"),
))

# Which --sections each source belongs to.
SECTIONS = OrderedDict((
    ("qc", ["fastp", "kraken2", "bracken"]),
    ("assembly", ["quast", "checkm2"]),
    ("typing", ["kleborate", "kaptive", "ectyper", "ezclermont"]),
    ("amr", ["amrfinderplus"]),
    ("plasmids", ["mobsuite"]),
    ("annotation", ["bakta", "prokka"]),
    ("phylogeny", ["phylogeny"]),
))

# Per-tool concatenated tables: raw key -> file name. 'beside' is the directory the file
# goes next to with --beside-source (the parent of the source directory).
TOOL_TABLES = OrderedDict((
    ("fastp", ("fastp_all.tsv", "fastp")),
    ("kraken2", ("kraken2_species_all.tsv", "kraken2")),
    ("bracken", ("bracken_all.tsv", "bracken")),
    ("quast", ("quast_all.tsv", "quast")),
    ("checkm2", ("checkm2_all.tsv", "checkm2")),
    ("bakta", ("bakta_all.tsv", "bakta")),
    ("prokka", ("prokka_all.tsv", "prokka")),
    ("kleborate", ("kleborate_all.tsv", "kleborate")),
    ("kaptive_K", ("kaptive_K_all.tsv", "kaptive")),
    ("kaptive_O", ("kaptive_O_all.tsv", "kaptive")),
    ("ectyper", ("ectyper_all.tsv", "ectyper")),
    ("ezclermont", ("ezclermont_all.tsv", "ezclermont")),
    ("amrfinderplus", ("amrfinderplus_all.tsv", "amrfinderplus")),
    ("mobsuite", ("mobtyper_all.tsv", "mobsuite")),
))

# AMRFinderPlus column names changed between versions; try each in order.
AMR_GENE_COLS = ("Element symbol", "Gene symbol")
AMR_TYPE_COLS = ("Element type", "Type")

# Default --amr-flag rules: NAME=REGEX, matched case-insensitively against Class + Subclass.
DEFAULT_AMR_FLAGS = OrderedDict((
    ("esbl_or_ampc", "CEPHALOSPORIN"),
    ("carbapenemase", "CARBAPENEM"),
    ("gentamicin_r", "GENTAMICIN"),
    ("colistin_r", "COLISTIN"),
))


# --------------------------------------------------------------------------- plumbing

class Dirs:
    """Resolves a source name to a directory: --dir override, else <results>/<default>."""

    def __init__(self, results, overrides, report):
        self.results = Path(results) if results else None
        self.overrides = {}
        self.report = report
        for item in overrides or []:
            if "=" not in item:
                report.problem("ignoring --dir '%s' (expected NAME=PATH)" % item)
                continue
            name, path = item.split("=", 1)
            name, path = name.strip(), path.strip()
            if name not in SOURCE_DIRS:
                report.problem("ignoring --dir '%s': unknown source (see --list-sources)" % name)
                continue
            p = Path(path)
            if not p.is_dir():
                report.problem("--dir %s=%s is not a directory; skipping that source" % (name, path))
                continue
            self.overrides[name] = p

    def get(self, name):
        if name in self.overrides:
            return self.overrides[name]
        if self.results is not None:
            p = self.results / SOURCE_DIRS[name]
            if p.is_dir():
                return p
        return None

    def explicit(self, name):
        return name in self.overrides


class Report:
    def __init__(self):
        self.lines = []
        self.problems = []

    def source(self, name, files, samples, note="", where=None):
        self.lines.append((name, len(files), len(samples), note, str(where) if where else "-"))

    def skipped(self, name):
        self.lines.append((name, 0, 0, "directory not found", "-"))

    def problem(self, text):
        self.problems.append(text)
        print("WARNING: %s" % text, file=sys.stderr)

    def write(self, path, nsamples, tables):
        with open(path, "w", encoding="utf-8") as fh:
            fh.write("aggregate_results.py summary\n" + "=" * 70 + "\n\n")
            fh.write("Isolates in master_table.tsv: %d\n\n" % nsamples)
            fh.write("%-16s %6s %8s  %-26s %s\n" % ("source", "files", "samples", "note", "directory"))
            fh.write("%-16s %6s %8s  %-26s %s\n" % ("-" * 16, "-" * 6, "-" * 8, "-" * 26, "-" * 30))
            for name, nfiles, ns, note, where in self.lines:
                fh.write("%-16s %6d %8d  %-26s %s\n" % (name, nfiles, ns, note, where))
            fh.write("\nTables written:\n")
            for label, p, nrows in tables:
                fh.write("  %-32s %7d rows  %s\n" % (label, nrows, p))
            if self.problems:
                fh.write("\nProblems (%d):\n" % len(self.problems))
                for text in self.problems:
                    fh.write("  - %s\n" % text)
            else:
                fh.write("\nNo problems reported.\n")


RAW = OrderedDict()           # tool-table key -> [(sample, OrderedDict row), ...]


def raw_add(key, sample, row):
    RAW.setdefault(key, []).append((sample, row))


def read_tsv(path, delimiter="\t"):
    try:
        with open(path, newline="", encoding="utf-8-sig", errors="replace") as fh:
            return list(csv.DictReader(fh, delimiter=delimiter))
    except OSError:
        return []


def first_col(row, names):
    for name in names:
        if name in row and row[name] not in (None, ""):
            return row[name]
    return ""


def sample_from(path, *suffixes):
    name = Path(path).name
    for suffix in suffixes:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return Path(name).stem


def clean_id(value):
    """Normalise a sample ID taken from inside a file (tools differ in what they print)."""
    value = Path(str(value).strip()).name
    for suffix in (".fasta.gz", ".fa.gz", ".fasta", ".fna", ".fa", ".gfa"):
        if value.endswith(suffix):
            value = value[: -len(suffix)]
            break
    for suffix in ("_contigs", ".contigs", ".hybrid", ".scaffolds", ".shovill"):
        if value.endswith(suffix):
            value = value[: -len(suffix)]
    return value


def num(value, cast=float):
    try:
        return cast(str(value).replace(",", "").strip())
    except (TypeError, ValueError):
        return None


def joined(values):
    return ";".join(sorted({v for v in values if v}))


def put(table, sample, **cols):
    table.setdefault(sample, {}).update(cols)


# --------------------------------------------------------------------------- readers

def read_fastp(dirs, wide, report):
    d = dirs.get("fastp")
    if d is None:
        return report.skipped("fastp")
    files = sorted(d.glob("*.json"))
    seen = set()
    for path in files:
        sample = sample_from(path, ".fastp.json", ".json")
        try:
            with open(path, encoding="utf-8") as fh:
                data = json.load(fh)
        except (OSError, ValueError):
            report.problem("could not parse %s" % path)
            continue
        before = data.get("summary", {}).get("before_filtering", {})
        after = data.get("summary", {}).get("after_filtering", {})
        put(wide, sample,
            reads_before=before.get("total_reads"),
            reads_after=after.get("total_reads"),
            q30_rate_after=after.get("q30_rate"),
            duplication_rate=data.get("duplication", {}).get("rate"))
        row = OrderedDict()
        for stage, block in (("before", before), ("after", after)):
            for key, value in block.items():
                if not isinstance(value, (dict, list)):
                    row["%s_%s" % (stage, key)] = value
        row["duplication_rate"] = data.get("duplication", {}).get("rate")
        row["insert_size_peak"] = data.get("insert_size", {}).get("peak")
        raw_add("fastp", sample, row)
        seen.add(sample)
    report.source("fastp", files, seen, where=d)


def _kraken_rows(path):
    """Yield (percent, clade_reads, taxon_reads, rank, taxid, name) tuples."""
    try:
        with open(path, newline="", encoding="utf-8", errors="replace") as fh:
            for row in csv.reader(fh, delimiter="\t"):
                if len(row) >= 6:
                    yield [c.strip() for c in row[:6]]
    except OSError:
        return


def read_kraken2(dirs, wide, report):
    d = dirs.get("kraken2")
    if d is None:
        return report.skipped("kraken2")
    files = sorted(d.glob("*report.txt"))
    seen = set()
    for path in files:
        sample = sample_from(path, ".kraken2.report.txt", ".report.txt", ".txt")
        best, unclassified = None, None
        for percent, clade, taxon, rank, taxid, name in _kraken_rows(path):
            if rank == "U":
                unclassified = num(percent)
            if rank == "S":
                reads = num(taxon, int) or 0
                if best is None or reads > best[0]:
                    best = (reads, name, num(percent))
                raw_add("kraken2", sample, OrderedDict((
                    ("rank", rank), ("taxid", taxid), ("name", name),
                    ("percent", num(percent)), ("clade_reads", num(clade, int)),
                    ("taxon_reads", num(taxon, int)))))
        put(wide, sample,
            kraken2_top_species=best[1] if best else None,
            kraken2_top_percent=best[2] if best else None,
            kraken2_unclassified_percent=unclassified)
        seen.add(sample)
    report.source("kraken2", files, seen, "species rows kept", where=d)


def read_bracken(dirs, wide, report):
    d = dirs.get("bracken")
    if d is None:
        return report.skipped("bracken")
    files = [p for p in sorted(d.glob("*")) if p.is_file() and p.suffix in (".txt", ".tsv")]
    seen = set()
    for path in files:
        sample = sample_from(path, ".bracken.report.txt", "_bracken.txt", ".tsv", ".txt")
        rows = read_tsv(path)
        if rows and "name" in rows[0] and "fraction_total_reads" in rows[0]:
            top = max(rows, key=lambda r: num(r.get("fraction_total_reads")) or 0)
            put(wide, sample,
                bracken_top_species=top.get("name"),
                bracken_top_fraction=num(top.get("fraction_total_reads")))
            for row in rows:
                raw_add("bracken", sample, OrderedDict(row))
        else:
            best = None
            for percent, _clade, taxon, rank, _taxid, name in _kraken_rows(path):
                if rank == "S":
                    reads = num(taxon, int) or 0
                    if best is None or reads > best[0]:
                        best = (reads, name, num(percent))
            put(wide, sample,
                bracken_top_species=best[1] if best else None,
                bracken_top_fraction=best[2] if best else None)
        seen.add(sample)
    report.source("bracken", files, seen, where=d)


QUAST_METRICS = OrderedDict((
    ("# contigs", "contigs"),
    ("Total length", "assembly_length"),
    ("N50", "n50"),
    ("GC (%)", "gc_percent"),
    ("Largest contig", "largest_contig"),
))


def read_quast(dirs, wide, report):
    d = dirs.get("quast")
    if d is None:
        return report.skipped("quast")
    files = sorted(d.glob("*.tsv"))
    seen = set()
    for path in files:
        sample = sample_from(path, ".tsv")
        cols, raw = {}, OrderedDict()
        try:
            with open(path, newline="", encoding="utf-8", errors="replace") as fh:
                for row in csv.reader(fh, delimiter="\t"):
                    if len(row) < 2:
                        continue
                    metric, value = row[0].strip(), row[-1]
                    raw[metric] = value
                    if metric in QUAST_METRICS:
                        cols[QUAST_METRICS[metric]] = num(value)
        except OSError:
            report.problem("could not read %s" % path)
            continue
        put(wide, sample, **cols)
        raw_add("quast", sample, raw)
        seen.add(sample)
    report.source("quast", files, seen, "all metrics in tool table", where=d)


def read_checkm2(dirs, wide, report):
    d = dirs.get("checkm2")
    if d is None:
        return report.skipped("checkm2")
    files = sorted(d.glob("*quality_report.tsv")) or sorted(d.glob("*.tsv"))
    seen = set()
    for path in files:
        for row in read_tsv(path):
            sample = clean_id(row.get("Name", ""))
            if not sample:
                continue
            put(wide, sample,
                completeness=num(row.get("Completeness")),
                contamination=num(row.get("Contamination")))
            raw_add("checkm2", sample, OrderedDict(row))
            seen.add(sample)
    report.source("checkm2", files, seen, where=d)


# Kleborate 3 prefixes its columns with the module name, e.g.
# 'klebsiella_pneumo_complex__mlst__ST'. Match on the suffix instead.
KLEBORATE_COLS = OrderedDict((
    ("kleborate_species", r"__species$"),
    ("ST", r"mlst__ST$"),
    ("kleborate_K_locus", r"K_locus$"),
    ("kleborate_O_type", r"O_type$"),
    ("resistance_score", r"resistance_score$"),
    ("virulence_score", r"virulence_score$"),
))


def read_kleborate(dirs, wide, report, extra_patterns=()):
    d = dirs.get("kleborate")
    if d is None:
        return report.skipped("kleborate")
    files = sorted(d.glob("*.txt")) + sorted(d.glob("*.tsv"))
    wanted = list(KLEBORATE_COLS.items()) + [(p, p) for p in extra_patterns]
    seen = set()
    for path in files:
        rows = read_tsv(path)
        if not rows:
            continue
        file_sample = sample_from(path, ".kleborate.txt", ".txt", ".tsv")
        multi = len(rows) > 1          # a combined Kleborate table, not one file per isolate
        for row in rows:
            inner = clean_id(next(iter(row.values())) or "")
            sample = inner if multi else (file_sample or inner)
            if not sample:
                continue
            cols = {}
            for out_name, pattern in wanted:
                for key in row:
                    if re.search(pattern, key):
                        cols[out_name] = row[key]
                        break
            put(wide, sample, **cols)
            raw_add("kleborate", sample, OrderedDict(row))
            seen.add(sample)
    report.source("kleborate", files, seen, where=d)


KAPTIVE_COLS = OrderedDict((
    ("locus", "Best match locus"),
    ("type", "Best match type"),
    ("confidence", "Match confidence"),
    ("problems", "Problems"),
))


def read_kaptive(dirs, wide, report):
    d = dirs.get("kaptive")
    if d is None:
        return report.skipped("kaptive")
    allfiles, seen = [], set()
    for letter in ("K", "O"):
        files = sorted(d.glob("*.kaptive_%s.tsv" % letter)) or sorted(d.glob("*_%s.tsv" % letter))
        allfiles += files
        for path in files:
            file_sample = sample_from(path, ".kaptive_%s.tsv" % letter, "_%s.tsv" % letter)
            rows = read_tsv(path)
            multi = len(rows) > 1
            for row in rows:
                inner = clean_id(first_col(row, ("Assembly", "assembly", "Name")))
                sample = inner if multi else (file_sample or inner)
                if not sample:
                    continue
                put(wide, sample, **{
                    "kaptive_%s_%s" % (letter, out): row.get(col, "")
                    for out, col in KAPTIVE_COLS.items()})
                raw_add("kaptive_%s" % letter, sample, OrderedDict(row))
                seen.add(sample)
    report.source("kaptive", allfiles, seen, "K and O read separately", where=d)


def read_ectyper(dirs, wide, report):
    d = dirs.get("ectyper")
    if d is None:
        return report.skipped("ectyper")
    files = sorted(d.rglob("*.tsv"))
    seen = set()
    for path in files:
        for row in read_tsv(path):
            sample = clean_id(first_col(row, ("Name", "name", "Sample")))
            if not sample:
                continue
            put(wide, sample,
                ectyper_serotype=first_col(row, ("Serotype", "serotype")),
                ectyper_o_type=first_col(row, ("O-type", "O_type")),
                ectyper_h_type=first_col(row, ("H-type", "H_type")))
            raw_add("ectyper", sample, OrderedDict(row))
            seen.add(sample)
    report.source("ectyper", files, seen, where=d)


def read_ezclermont(dirs, wide, report):
    d = dirs.get("ezclermont")
    if d is None:
        return report.skipped("ezclermont")
    files = sorted(p for p in d.rglob("*") if p.is_file() and p.suffix in (".txt", ".tsv"))
    seen = set()
    for path in files:
        sample = sample_from(path, ".ezclermont.txt", ".tsv", ".txt")
        try:
            text = open(path, encoding="utf-8", errors="replace").read().strip()
        except OSError:
            continue
        if not text:
            continue
        last = text.splitlines()[-1]
        phylogroup = last.split("\t")[-1].strip()
        put(wide, sample, ezclermont_phylogroup=phylogroup)
        raw_add("ezclermont", sample,
                OrderedDict((("phylogroup", phylogroup), ("raw_line", last))))
        seen.add(sample)
    report.source("ezclermont", files, seen, where=d)


def read_amrfinder(dirs, wide, report, flags):
    d = dirs.get("amrfinderplus")
    if d is None:
        report.skipped("amrfinderplus")
        return []
    files = [p for p in sorted(d.glob("*.tsv")) if not p.name.endswith("-mutations.tsv")]
    long_rows, seen = [], set()
    for path in files:
        file_sample = sample_from(path, ".tsv")
        rows = read_tsv(path)
        by_sample = {}
        for row in rows:
            inner = clean_id(first_col(row, ("Name", "Protein id", "Protein identifier")))
            # A per-isolate file is named after the isolate; a combined file has the ID inside.
            sample = file_sample if file_sample else inner
            if len({clean_id(first_col(r, ("Name",))) for r in rows if first_col(r, ("Name",))}) > 1:
                sample = inner or file_sample
            by_sample.setdefault(sample, []).append(row)
        for sample, srows in by_sample.items():
            if not sample:
                continue
            genes, classes, subclasses = [], [], []
            for row in srows:
                gene = first_col(row, AMR_GENE_COLS)
                etype = first_col(row, AMR_TYPE_COLS)
                klass, subclass = row.get("Class", ""), row.get("Subclass", "")
                if etype.upper() == "AMR":
                    genes.append(gene)
                    classes.append(klass)
                subclasses.append("%s %s" % (klass, subclass))
                raw_add("amrfinderplus", sample, OrderedDict(row))
                long_rows.append(OrderedDict((
                    ("sample", sample), ("gene", gene), ("element_type", etype),
                    ("element_subtype", first_col(row, ("Element subtype", "Subtype"))),
                    ("class", klass), ("subclass", subclass),
                    ("contig", first_col(row, ("Contig id", "Contig"))),
                    ("coverage", first_col(row, ("% Coverage of reference sequence",
                                                 "% Coverage of reference"))),
                    ("identity", first_col(row, ("% Identity to reference sequence",
                                                 "% Identity to reference"))))))
            haystack = " | ".join(subclasses).upper()
            cols = {"amr_n_genes": len(genes), "amr_genes": joined(genes),
                    "amr_classes": joined(classes)}
            for name, pattern in flags.items():
                cols[name] = int(bool(re.search(pattern, haystack, re.IGNORECASE)))
            put(wide, sample, **cols)
            seen.add(sample)
    report.source("amrfinderplus", files, seen, "%d gene hits" % len(long_rows), where=d)
    return long_rows


MOBILITY_KEYS = ("conjugative", "mobilizable", "non-mobilizable")


def read_mobsuite(dirs, wide, report):
    d = dirs.get("mobsuite")
    if d is None:
        report.skipped("mobsuite")
        return []
    files = sorted(d.glob("*/mobtyper_results.txt")) or sorted(d.glob("*mobtyper*"))
    long_rows, seen = [], set()
    for path in files:
        if path.name == "mobtyper_results.txt":
            sample = path.parent.name          # mobsuite/<isolate>/mobtyper_results.txt
        else:                                  # or one renamed file per isolate
            sample = sample_from(path, ".mobtyper_results.txt", ".mobtyper.txt",
                                 "_mobtyper.txt", ".txt", ".tsv")
        rows = read_tsv(path)
        sizes, reps, clusters = [], [], []
        mobility = Counter()
        for row in rows:
            size = num(row.get("size"), int)
            rep = row.get("rep_type(s)", "")
            cluster = row.get("primary_cluster_id", "")
            mob = (row.get("predicted_mobility", "") or "").strip().lower()
            sizes.append(size or 0)
            reps += [r.strip() for r in rep.split(",") if r.strip()]
            clusters.append(cluster)
            mobility[mob] += 1
            raw_add("mobsuite", sample, OrderedDict(row))
            long_rows.append(OrderedDict((
                ("sample", sample), ("plasmid_id", row.get("sample_id", "")),
                ("primary_cluster_id", cluster),
                ("secondary_cluster_id", row.get("secondary_cluster_id", "")),
                ("size", size), ("gc", row.get("gc", "")), ("rep_types", rep),
                ("relaxase_types", row.get("relaxase_type(s)", "")),
                ("predicted_mobility", mob))))
        cols = {"n_plasmids": len(rows), "plasmid_total_bp": sum(sizes),
                "plasmid_rep_types": joined(reps), "plasmid_clusters": joined(clusters)}
        for key in MOBILITY_KEYS:
            cols["n_plasmids_%s" % key.replace("-", "_")] = mobility.get(key, 0)
        put(wide, sample, **cols)
        seen.add(sample)
    report.source("mobsuite", files, seen, "%d plasmids" % len(long_rows), where=d)
    return long_rows


ANNOT_KEYS = ("cds", "cdss", "trna", "trnas", "rrna", "rrnas", "contigs", "sequences", "count")


def read_annotation(dirs, wide, report):
    for tool in ("bakta", "prokka"):
        d = dirs.get(tool)
        if d is None:
            report.skipped(tool)
            continue
        files = sorted(d.rglob("*.txt"))
        seen = set()
        for path in files:
            sample = path.parent.name if path.parent != d else sample_from(path, ".txt")
            cols, raw = {}, OrderedDict()
            try:
                for line in open(path, encoding="utf-8", errors="replace"):
                    if ":" not in line:
                        continue
                    key, value = line.split(":", 1)
                    key = key.strip().lower().replace(" ", "_").replace("(s)", "s")
                    raw[key] = value.strip()
                    if key in ANNOT_KEYS:
                        cols["%s_%s" % (tool, key)] = num(value, int)
            except OSError:
                continue
            if cols:
                put(wide, sample, **cols)
                raw_add(tool, sample, raw)
                seen.add(sample)
        report.source(tool, files, seen, "feature counts only", where=d)


# ------------------------------------------------------- phylogeny (pairs -> per isolate)

def _read_snpdists_matrix(path):
    """snp-dists matrix -> (samples, {(a, b): snps}). Drops the Snippy 'Reference' row."""
    try:
        with open(path, newline="", encoding="utf-8", errors="replace") as fh:
            rows = [r for r in csv.reader(fh, delimiter="\t") if r]
    except OSError:
        return [], {}
    if len(rows) < 2:
        return [], {}
    header = [c.strip() for c in rows[0][1:]]
    dist = {}
    for row in rows[1:]:
        a = row[0].strip()
        for b, value in zip(header, row[1:]):
            n = num(value, int)
            if n is not None and a != b:
                dist[(a, b)] = n
    return [s for s in header if s and s.lower() != "reference"], dist


def read_phylogeny(dirs, wide, report, threshold):
    d = dirs.get("phylogeny")
    if d is None:
        report.skipped("phylogeny")
        return []
    matrices = sorted(d.glob("*/*.snp_distances.tsv")) or sorted(d.glob("*.snp_distances.tsv"))
    long_rows, seen = [], set()
    for path in matrices:
        group = path.parent.name if path.parent != d else sample_from(path, ".snp_distances.tsv")
        samples, dist = _read_snpdists_matrix(path)
        if not samples:
            report.problem("no usable SNP matrix in %s" % path)
            continue
        for a in samples:
            others = [(dist[(a, b)], b) for b in samples if b != a and (a, b) in dist]
            cols = {"phylo_group": group, "phylo_group_n": len(samples)}
            if others:
                nearest = min(others)
                cols["nearest_sample"] = nearest[1]
                cols["nearest_snps"] = nearest[0]
                cols["n_within_%d_snps" % threshold] = sum(1 for x, _ in others if x <= threshold)
            put(wide, a, **cols)
            seen.add(a)
        for (a, b), snps in sorted(dist.items()):
            if a < b and a in samples and b in samples:
                long_rows.append(OrderedDict((("group", group), ("sample1", a), ("sample2", b),
                                              ("snps", snps), ("method", "snp-dists"))))
    for path in sorted(d.glob("*/*.ska_distances.tsv")) + sorted(d.glob("*.ska_distances.tsv")):
        group = path.parent.name if path.parent != d else sample_from(path, ".ska_distances.tsv")
        for row in read_tsv(path):
            a, b = row.get("Sample1", ""), row.get("Sample2", "")
            if a and b:
                long_rows.append(OrderedDict((("group", group), ("sample1", a), ("sample2", b),
                                              ("snps", row.get("Mismatches", row.get("Distance", ""))),
                                              ("method", "ska distance"))))
    trees = {p.parent.name: p for p in sorted(d.glob("*/*.treefile"))}
    for cols in wide.values():
        group = cols.get("phylo_group")
        if group and group in trees:
            cols["tree_file"] = str(trees[group])
    report.source("phylogeny", matrices, seen,
                  "%d pairs, %d tree(s)" % (len(long_rows), len(trees)), where=d)
    return long_rows


# --------------------------------------------------------------------------- extras

def read_metadata(path, wide, report):
    delim = "\t" if str(path).endswith((".tsv", ".txt")) else ","
    rows = read_tsv(path, delimiter=delim)
    if not rows:
        report.problem("no rows read from --metadata %s" % path)
        return
    header = list(rows[0])
    id_col = next((c for c in ("sample", "id", "Name", "strain") if c in header), header[0])
    seen = set()
    for row in rows:
        sample = clean_id(row.get(id_col, ""))
        if sample:
            put(wide, sample, **{k: v for k, v in row.items() if k != id_col})
            seen.add(sample)
    report.source("metadata", [path], seen, "ID column '%s'" % id_col, where=Path(path).parent)


def samplesheet_ids(path, report):
    rows = read_tsv(path, delimiter="," if str(path).endswith(".csv") else "\t")
    if not rows:
        report.problem("no rows read from --samplesheet %s" % path)
        return []
    header = list(rows[0])
    id_col = next((c for c in ("sample", "id", "Sample") if c in header), header[0])
    out, seen = [], set()
    for row in rows:
        sample = clean_id(row.get(id_col, ""))
        if sample and sample not in seen:      # a samplesheet may list a sample twice (R1/R2 rows)
            out.append(sample)
            seen.add(sample)
    return out


# --------------------------------------------------------------------------- output

def fmt(value, na):
    """Missing -> na; a float with no fractional part -> int (QUAST counts, plasmid sizes)."""
    if value is None or value == "":
        return na
    if isinstance(value, float) and value.is_integer():
        return int(value)
    return value


def write_table(path, rows, columns, na):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        for row in rows:
            writer.writerow([fmt(row.get(c), na) for c in columns])
    return len(rows)


def write_tool_tables(outdir, dirs, beside, na, report):
    """One file per tool: every row, the tool's own columns, plus a leading 'sample'."""
    tables = []
    for key, (name, source) in TOOL_TABLES.items():
        entries = RAW.get(key)
        if not entries:
            continue
        columns = ["sample"]
        for _sample, row in entries:                 # union of headers, first-seen order
            for col in row:
                if col not in columns:
                    columns.append(col)
        rows = []
        for sample, row in entries:
            merged = OrderedDict(sample=sample)
            merged.update(row)
            rows.append(merged)
        target = outdir / name
        if beside:
            d = dirs.get(source)
            if d is not None:
                target = d.parent / name
            else:
                report.problem("--beside-source: no directory for %s, writing %s to --outdir"
                               % (source, name))
        tables.append(("%s (tool table)" % name[:-4], target,
                       write_table(target, rows, columns, na)))
    return tables


def order_columns(wide, leading):
    seen = list(leading)
    for cols in wide.values():
        for key in cols:
            if key not in seen:
                seen.append(key)
    present = {k for cols in wide.values() for k in cols}
    return [c for c in seen if c in present or c == "sample"]


LEADING = ["sample", "phylo_group", "phylo_group_n", "kleborate_species", "ST",
           "kleborate_K_locus", "kleborate_O_type", "kaptive_K_locus", "kaptive_K_type",
           "kaptive_K_confidence", "kaptive_O_locus", "kaptive_O_type", "kaptive_O_confidence",
           "resistance_score", "virulence_score", "amr_n_genes", "n_plasmids",
           "completeness", "contamination", "contigs", "assembly_length", "n50"]


# --------------------------------------------------------------------------- CLI

def parse_amr_flags(values, report):
    flags = OrderedDict(DEFAULT_AMR_FLAGS)
    for item in values or []:
        if "=" not in item:
            report.problem("ignoring --amr-flag '%s' (expected NAME=REGEX)" % item)
            continue
        name, pattern = item.split("=", 1)
        name = name.strip()
        if pattern.strip():
            flags[name] = pattern.strip()
        else:
            flags.pop(name, None)              # NAME= removes a default flag
    return flags


def build_parser():
    p = argparse.ArgumentParser(
        prog="aggregate_results.py",
        description="Aggregate neonatal-amr-nf results into one table per unit of observation.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Source names for --dir: " + ", ".join(SOURCE_DIRS) +
               "\nSee --list-sources for where each one is looked for, and --describe for "
               "what can and cannot be combined.")
    p.add_argument("--results", "-r", metavar="DIR",
                   help="the pipeline's --outdir (the folder holding 01_qc, 02_assembly, ...). "
                        "Optional if you give --dir for every source you want.")
    p.add_argument("--dir", action="append", metavar="SOURCE=PATH", dest="dirs",
                   help="read one source from this directory instead of under --results. "
                        "Repeatable, e.g. --dir kleborate=/scratch/kleborate")
    p.add_argument("--outdir", "-o", metavar="DIR",
                   help="where to write the tables [<results>/09_aggregate; required without --results]")
    p.add_argument("--samplesheet", metavar="CSV",
                   help="defines the expected isolates; extras and missing ones are reported")
    p.add_argument("--metadata", metavar="TSV",
                   help="extra per-isolate columns to merge (ward, year, AST, ...)")
    p.add_argument("--sections", metavar="LIST", default="all",
                   help="comma-separated subset of: " + ",".join(SECTIONS) + ",all [all]")
    p.add_argument("--tool-tables", action="store_true",
                   help="also write one lossless table per tool (kleborate_all.tsv, "
                        "amrfinderplus_all.tsv, mobtyper_all.tsv, ...)")
    p.add_argument("--beside-source", action="store_true",
                   help="with --tool-tables, write each one next to its source directory "
                        "(04_amr/amrfinderplus_all.tsv) instead of into --outdir")
    p.add_argument("--long", action="store_true",
                   help="also write the long tables (gene hits, plasmids, SNP pairs)")
    p.add_argument("--snp-threshold", type=int, default=20, metavar="N",
                   help="for the n_within_<N>_snps column [20]")
    p.add_argument("--amr-flag", action="append", metavar="NAME=REGEX",
                   help="add or replace an AMR flag column (NAME= with no regex removes one). "
                        "Defaults: " + ", ".join("%s=%s" % kv for kv in DEFAULT_AMR_FLAGS.items()))
    p.add_argument("--kleborate-extra", action="append", metavar="REGEX", default=[],
                   help="extra Kleborate columns to keep, matched by regex on the column name")
    p.add_argument("--na", default="NA", metavar="STR", help="text for missing values [NA]")
    p.add_argument("--prefix", default="master_table", metavar="NAME",
                   help="base name of the wide table [master_table]")
    p.add_argument("--strict", action="store_true",
                   help="exit non-zero if a samplesheet isolate is missing from the table")
    p.add_argument("--list-sources", action="store_true",
                   help="print the source names, their default locations and tool tables, then exit")
    p.add_argument("--describe", action="store_true",
                   help="print what can and cannot be combined, then exit")
    return p


def list_sources():
    tool_by_source = {}
    for key, (name, source) in TOOL_TABLES.items():
        tool_by_source.setdefault(source, []).append(name)
    section_of = {s: sec for sec, names in SECTIONS.items() for s in names}
    print("%-16s %-24s %-12s %s" % ("source", "default path under --results", "--sections", "tool table"))
    print("%-16s %-24s %-12s %s" % ("-" * 16, "-" * 24, "-" * 12, "-" * 28))
    for name, rel in SOURCE_DIRS.items():
        print("%-16s %-24s %-12s %s" % (name, rel, section_of.get(name, "-"),
                                        ", ".join(tool_by_source.get(name, ["-"]))))
    print("\nOverride any of them with --dir SOURCE=PATH (repeatable).")


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.describe:
        print(DESCRIBE.strip())
        return 0
    if args.list_sources:
        list_sources()
        return 0
    if not args.results and not args.dirs:
        parser.error("give --results, or --dir SOURCE=PATH for at least one source "
                     "(--list-sources shows the names)")
    if args.results and not Path(args.results).is_dir():
        parser.error("--results %s is not a directory" % args.results)
    if not args.results and not args.outdir:
        parser.error("--outdir is required when --results is not given")
    if args.beside_source and not args.tool_tables:
        parser.error("--beside-source only applies with --tool-tables")

    report = Report()
    dirs = Dirs(args.results, args.dirs, report)
    outdir = Path(args.outdir) if args.outdir else Path(args.results) / "09_aggregate"

    sections = {s.strip() for s in args.sections.split(",") if s.strip()}
    if "all" in sections:
        sections = set(SECTIONS)
    unknown = sections - set(SECTIONS)
    if unknown:
        parser.error("unknown --sections: %s (choose from %s,all)"
                     % (",".join(sorted(unknown)), ",".join(SECTIONS)))
    # A source named with --dir is wanted, whatever --sections says.
    for sec, names in SECTIONS.items():
        if any(dirs.explicit(n) for n in names):
            sections.add(sec)

    wide = {}
    amr_long = plasmid_long = snp_long = []
    if "qc" in sections:
        read_fastp(dirs, wide, report)
        read_kraken2(dirs, wide, report)
        read_bracken(dirs, wide, report)
    if "assembly" in sections:
        read_quast(dirs, wide, report)
        read_checkm2(dirs, wide, report)
    if "typing" in sections:
        read_kleborate(dirs, wide, report, args.kleborate_extra)
        read_kaptive(dirs, wide, report)
        read_ectyper(dirs, wide, report)
        read_ezclermont(dirs, wide, report)
    if "amr" in sections:
        amr_long = read_amrfinder(dirs, wide, report, parse_amr_flags(args.amr_flag, report))
    if "plasmids" in sections:
        plasmid_long = read_mobsuite(dirs, wide, report)
    if "annotation" in sections:
        read_annotation(dirs, wide, report)
    if "phylogeny" in sections:
        snp_long = read_phylogeny(dirs, wide, report, args.snp_threshold)
    if args.metadata:
        read_metadata(Path(args.metadata), wide, report)

    expected = samplesheet_ids(Path(args.samplesheet), report) if args.samplesheet else []
    missing = [s for s in expected if s not in wide]
    extra = [s for s in sorted(wide) if expected and s not in expected]
    if missing:
        report.problem("%d samplesheet isolate(s) have no results at all: %s"
                       % (len(missing), ", ".join(missing[:10])))
    if extra:
        report.problem("%d isolate(s) in the results are not in the samplesheet: %s"
                       % (len(extra), ", ".join(extra[:10])))
    for sample in expected:
        wide.setdefault(sample, {})
    if not wide:
        report.problem("no results found - check --results / --dir (see --list-sources)")

    rows = []
    for sample in (expected if expected else sorted(wide)):
        row = OrderedDict(sample=sample)
        row.update(wide[sample])
        rows.append(row)
    columns = order_columns(wide, LEADING)

    outdir.mkdir(parents=True, exist_ok=True)
    master = outdir / ("%s.tsv" % args.prefix)
    tables = [("master_table (isolates)", master, write_table(master, rows, columns, args.na))]

    if args.tool_tables:
        tables += write_tool_tables(outdir, dirs, args.beside_source, args.na, report)
    if args.long:
        for label, name, data in (("amr_genes_long (gene hits)", "amr_genes_long.tsv", amr_long),
                                  ("plasmids_long (plasmids)", "plasmids_long.tsv", plasmid_long),
                                  ("snp_distances_long (pairs)", "snp_distances_long.tsv", snp_long)):
            if data:
                path = outdir / name
                tables.append((label, path, write_table(path, data, list(data[0]), args.na)))

    summary = outdir / "aggregate_summary.txt"
    report.write(summary, len(rows), tables)

    print("Isolates: %d   Columns: %d" % (len(rows), len(columns)))
    for label, path, nrows in tables:
        print("  %-30s %7d rows  ->  %s" % (label, nrows, path))
    print("  %-30s %7s        ->  %s" % ("summary", "", summary))

    if args.strict and missing:
        print("ERROR: --strict and %d isolate(s) missing" % len(missing), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
