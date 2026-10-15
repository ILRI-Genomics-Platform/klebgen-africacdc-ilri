#!/usr/bin/env python3
"""
fetch_references.py - download the phylogeny reference genomes once, into a shared read-only folder.

The pipeline (--phylo_grouping kleborate, the default) maps each species/ST group to a reference listed in
conf/phylo_references.csv and reads it from <dir>/<accession>.fasta. Run this once as the trainer, with
network access (e.g. on the login node); participants and compute nodes only read the files.

    bin/fetch_references.py                         # download what is missing, then lock the folder
    bin/fetch_references.py --check                 # only report present / missing (changes nothing)
    bin/fetch_references.py --group amr-course      # readable only by this Unix group
    bin/fetch_references.py --dir /other/folder     # then pass the same --phylo_reference_dir to the pipeline

Sources: NCBI E-utilities (efetch), then ENA if NCBI fails. Writes MANIFEST.tsv (accession, name, species,
length, contigs, md5, source, date) next to the FASTA files.
"""

import argparse
import csv
import datetime
import grp
import hashlib
import io
import os
import stat
import sys
import time
import urllib.request

PIPELINE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_SHEET = os.path.join(PIPELINE_DIR, "conf", "phylo_references.csv")
DEFAULT_DIR = "/var/scratch/global/jjuma/neonatal-amr/shared/references"

SOURCES = [
    ("NCBI", "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={acc}&rettype=fasta&retmode=text"),
    ("ENA", "https://www.ebi.ac.uk/ena/browser/api/fasta/{acc}?download=true"),
]


def read_sheet(path):
    """Rows of the reference sheet; lines starting with '#' are comments."""
    with open(path) as handle:
        lines = [line for line in handle if line.strip() and not line.lstrip().startswith("#")]
    rows = list(csv.DictReader(io.StringIO("".join(lines))))
    if not rows:
        sys.exit(f"ERROR: no rows in {path}")
    for col in ("species", "st", "name", "accession"):
        if col not in rows[0]:
            sys.exit(f"ERROR: {path} has no '{col}' column (expected species,st,name,accession)")
    return [{k: (v or "").strip() for k, v in row.items()} for row in rows]


def check_fasta(data, acc):
    """Return (length, contigs) or raise ValueError if this is not a nucleotide FASTA."""
    text = data.decode("ascii", errors="replace")
    if not text.startswith(">"):
        raise ValueError(f"not a FASTA file (starts with {text[:60]!r})")
    contigs, length = 0, 0
    for line in text.splitlines():
        if line.startswith(">"):
            contigs += 1
        else:
            seq = line.strip().upper()
            bad = set(seq) - set("ACGTNRYKMSWBDHV")
            if bad:
                raise ValueError(f"unexpected characters {''.join(sorted(bad))} in the sequence")
            length += len(seq)
    if length < 100_000:
        raise ValueError(f"only {length} bp: is {acc} a complete genome or chromosome?")
    return length, contigs


def download(acc):
    errors = []
    for name, url in SOURCES:
        for attempt in (1, 2, 3):
            try:
                req = urllib.request.Request(url.format(acc=acc), headers={"User-Agent": "neonatal-amr-nf"})
                with urllib.request.urlopen(req, timeout=120) as resp:
                    data = resp.read()
                length, contigs = check_fasta(data, acc)
                return data, name, length, contigs
            except Exception as exc:  # network error or bad content: retry, then next source
                errors.append(f"{name} (attempt {attempt}): {exc}")
                if isinstance(exc, ValueError):
                    break
                time.sleep(3 * attempt)
    raise RuntimeError("; ".join(errors))


def lock(folder, group):
    """Read-only for everyone (0444/0555) or for one group (0440/0550); nobody keeps write permission."""
    gid = grp.getgrnam(group).gr_gid if group else -1
    for root, dirs, files in os.walk(folder):
        for f in files:
            p = os.path.join(root, f)
            if group:
                os.chown(p, -1, gid)
            os.chmod(p, 0o440 if group else 0o444)
    if group:
        os.chown(folder, -1, gid)
    os.chmod(folder, 0o550 if group else 0o555)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sheet", default=DEFAULT_SHEET, help=f"reference sheet (default: {DEFAULT_SHEET})")
    ap.add_argument("--dir", default=DEFAULT_DIR, help=f"reference folder (default: {DEFAULT_DIR})")
    ap.add_argument("--check", action="store_true", help="only report present / missing; change nothing")
    ap.add_argument("--force", action="store_true", help="download again even if the file exists")
    ap.add_argument("--group", help="make the files readable only by this Unix group")
    ap.add_argument("--no-lock", action="store_true", help="leave permissions unchanged")
    args = ap.parse_args()

    rows = read_sheet(args.sheet)
    accessions = {}
    for row in rows:
        accessions.setdefault(row["accession"], row)
    folder = os.path.abspath(args.dir)
    print(f"Reference sheet : {args.sheet} ({len(rows)} rows, {len(accessions)} genomes)")
    print(f"Reference folder: {folder}")

    missing = [a for a in accessions if args.force or not os.path.isfile(os.path.join(folder, f"{a}.fasta"))]
    for acc, row in accessions.items():
        state = "missing" if acc in missing else "present"
        print(f"  {state:8} {acc:14} {row['name']:16} {row['species']}{' ' + row['st'] if row['st'] else ''}")
    if args.check:
        sys.exit(1 if missing else 0)
    if not missing:
        print("All reference genomes present.")
    else:
        try:
            os.makedirs(folder, exist_ok=True)
            os.chmod(folder, os.stat(folder).st_mode | stat.S_IWUSR)   # unlock after an earlier run
        except PermissionError:
            sys.exit(f"ERROR: cannot create or write to {folder}; use --dir with a folder you own")
        manifest_path = os.path.join(folder, "MANIFEST.tsv")
        manifest = {}
        if os.path.isfile(manifest_path):
            with open(manifest_path) as handle:
                for rec in csv.DictReader(handle, delimiter="\t"):
                    manifest[rec["accession"]] = rec
            os.chmod(manifest_path, 0o644)

        failed = []
        for acc in missing:
            row = accessions[acc]
            print(f"Downloading {acc} ({row['name']})...", flush=True)
            try:
                data, source, length, contigs = download(acc)
            except RuntimeError as exc:
                print(f"  FAILED: {exc}", file=sys.stderr)
                failed.append(acc)
                continue
            target = os.path.join(folder, f"{acc}.fasta")
            if os.path.exists(target):
                os.chmod(target, 0o644)
            tmp = target + ".part"
            with open(tmp, "wb") as handle:
                handle.write(data)
            os.replace(tmp, target)
            manifest[acc] = {
                "accession": acc, "name": row["name"], "species": row["species"], "length_bp": length,
                "contigs": contigs, "md5": hashlib.md5(data).hexdigest(), "source": source,
                "downloaded": datetime.date.today().isoformat(),
            }
            print(f"  {length:,} bp in {contigs} sequence(s) from {source}")

        with open(manifest_path, "w", newline="") as handle:
            cols = ["accession", "name", "species", "length_bp", "contigs", "md5", "source", "downloaded"]
            writer = csv.DictWriter(handle, fieldnames=cols, delimiter="\t")
            writer.writeheader()
            for acc in sorted(manifest):
                writer.writerow({c: manifest[acc].get(c, "") for c in cols})
        if failed:
            sys.exit(f"ERROR: could not download {', '.join(failed)}. Check the network (NCBI and ENA) and "
                     f"the accession(s) in {args.sheet}, then run this script again.")

    if not args.no_lock and os.path.isdir(folder):
        lock(folder, args.group)
        print(f"Locked: readable by {'group ' + args.group if args.group else 'everyone'}; nobody can write. "
              "Run this script again to add references.")
    print(f"Pipeline option: --phylo_reference_dir {folder}  (already set in conf/slurm.config and conf/local.config)")


if __name__ == "__main__":
    main()
