#!/usr/bin/env python3
"""
make_assembly_samplesheet.py - build the pipeline's assembly samplesheet from a metadata table and a folder of
contigs, and a matching metadata table for the R scripts (transmission pairs, figures, BactDating).

    bin/make_assembly_samplesheet.py \\
        --metadata ethiopia_metadata.tsv \\
        --contigs-dir /path/to/contigs \\
        --out samplesheet_assemblies.csv

    nextflow run main.nf -profile singularity -c conf/slurm.config \\
        --input samplesheet_assemblies.csv --outdir results_ethiopia      # --input_type assembly is detected

One samplesheet per species (each with its own metadata table), to run K. pneumoniae and E. coli separately:

    bin/make_assembly_samplesheet.py --metadata ethiopia_metadata.tsv --contigs-dir /path/to/contigs \\
        --out ethiopia/samplesheet_assemblies.csv --split-by-species
    #   ethiopia/samplesheet_assemblies_kpneumoniae.csv   ethiopia/metadata_assemblies_kpneumoniae.tsv
    #   ethiopia/samplesheet_assemblies_ecoli.csv         ethiopia/metadata_assemblies_ecoli.tsv
    #   ethiopia/samplesheet_assemblies_other_species.tsv (organisms the pipeline does not type; for the record)

Only one species:  --species kpneumoniae   or   --species ecoli

Defaults match a tab-delimited table with these columns (change them with the options below):
    entity:Ethiopia_sample_id  -> sample ID; the contigs file is <ID>_contigs.fasta (also <ID>.fasta, .fa, .fna)
    Organism_detected_phenotypically -> species: 'K. pneumoniae', 'Klebsiella pneumoniae' -> kpneumoniae,
                                        'E. coli', 'Escherichia coli' -> ecoli; anything else is left out
    Year, Month, Day, Collection_date, Site -> collection_year, collection_date, ward in the metadata output

Samples are reported and left out when their contigs file is missing or compressed, or their organism is not
K. pneumoniae or E. coli (those are listed in <out>_other_species.tsv); contigs files without a metadata row are
reported too.
"""

import argparse
import csv
import datetime
import os
import re
import sys

MONTHS = {m: i for i, m in enumerate(
    ["jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"], start=1)}


def species_code(value):
    """kpneumoniae / ecoli, or None for organisms the pipeline does not type."""
    s = re.sub(r"[^a-z]", "", (value or "").lower())
    if s in ("kpneumoniae", "klebsiellapneumoniae", "kpn", "kp") or "pneumoniae" in s or "variicola" in s:
        return "kpneumoniae"
    if s in ("ecoli", "escherichiacoli", "eco"):
        return "ecoli"
    return None


def find_contigs(folder, sample, suffix):
    """Path of the sample's contigs file: <ID><suffix>, then <ID>.fasta / .fa / .fna / _contigs.fa."""
    for name in [sample + suffix] + [sample + s for s in (".fasta", ".fa", ".fna", "_contigs.fa", "_contigs.fna")]:
        path = os.path.join(folder, name)
        if os.path.isfile(path):
            return path
        if os.path.isfile(path + ".gz"):
            return path + ".gz"
    return None


def iso_date(row, year_col, month_col, day_col, date_col):
    """YYYY-MM-DD, YYYY-MM or YYYY from Year/Month/Day (or an ISO Collection_date); '' if unknown."""
    y = (row.get(year_col) or "").strip()
    m = (row.get(month_col) or "").strip()
    d = (row.get(day_col) or "").strip()
    raw = (row.get(date_col) or "").strip()
    if re.fullmatch(r"\d{4}-\d{2}(-\d{2})?", raw):
        return raw
    if not re.fullmatch(r"\d{4}", y):
        found = re.search(r"(19|20)\d{2}", raw)
        return found.group(0) if found else ""
    month = int(m) if m.isdigit() else MONTHS.get(m[:3].lower())
    if not month or not 1 <= month <= 12:
        return y
    if d.isdigit():
        try:
            return datetime.date(int(y), month, int(d)).isoformat()
        except ValueError:
            pass
    return f"{y}-{month:02d}"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--metadata", required=True, help="metadata table (tab-delimited; .csv is read as comma-separated)")
    ap.add_argument("--contigs-dir", required=True, help="folder with the contigs FASTA files")
    ap.add_argument("--out", default="samplesheet_assemblies.csv", help="samplesheet to write [%(default)s]")
    ap.add_argument("--metadata-out", default=None,
                    help="metadata for the R scripts [<out folder>/metadata_assemblies.tsv]")
    ap.add_argument("--id-col", default="entity:Ethiopia_sample_id", help="sample ID column [%(default)s]")
    ap.add_argument("--species-col", default="Organism_detected_phenotypically", help="organism column [%(default)s]")
    ap.add_argument("--suffix", default="_contigs.fasta", help="contigs file = <ID><suffix> [%(default)s]")
    ap.add_argument("--group-col", default=None, help="optional column to use as phylo_group")
    ap.add_argument("--split-by-species", action="store_true",
                    help="write one samplesheet and one metadata table per species (<out>_kpneumoniae.csv, <out>_ecoli.csv)")
    ap.add_argument("--species", choices=["all", "kpneumoniae", "ecoli"], default="all",
                    help="keep only this species [%(default)s]")
    ap.add_argument("--year-col", default="Year")
    ap.add_argument("--month-col", default="Month")
    ap.add_argument("--day-col", default="Day")
    ap.add_argument("--date-col", default="Collection_date")
    ap.add_argument("--ward-col", default="Site", help="column used as 'ward' in the R scripts [%(default)s]")
    args = ap.parse_args()

    delim = "," if args.metadata.lower().endswith(".csv") else "\t"
    with open(args.metadata, newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle, delimiter=delim))
    if not rows:
        sys.exit(f"ERROR: no rows in {args.metadata}")
    cols = list(rows[0].keys())
    for col in (args.id_col, args.species_col):
        if col not in cols:
            sys.exit(f"ERROR: column '{col}' not in {args.metadata}. Columns: {', '.join(cols)}")

    folder = os.path.abspath(args.contigs_dir)
    sheet, meta_out, skipped, seen, others, not_selected = [], [], [], set(), [], set()
    for row in rows:
        sample = (row.get(args.id_col) or "").strip()
        if not sample:
            continue
        if sample in seen:
            skipped.append((sample, "duplicate sample ID in the metadata"))
            continue
        seen.add(sample)
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", sample):
            skipped.append((sample, "ID may only contain letters, digits, '_', '.' and '-'"))
            continue
        organism = (row.get(args.species_col) or "").strip()
        species = species_code(organism)
        contigs = find_contigs(folder, sample, args.suffix)
        if species is None:
            others.append({"sample": sample, "organism": organism or "(empty)", "fasta": contigs or "(not found)"})
            continue
        if args.species != "all" and species != args.species:
            not_selected.add(sample)
            continue
        if contigs is None:
            skipped.append((sample, f"no contigs file ({sample}{args.suffix}) in {folder}"))
            continue
        if contigs.endswith(".gz"):
            skipped.append((sample, f"{os.path.basename(contigs)} is compressed: gunzip -k it first"))
            continue
        group = (row.get(args.group_col) or "").strip() if args.group_col else ""
        sheet.append({"sample": sample, "species": species, "fasta": contigs, "phylo_group": group})

        out = {"sample": sample,
               "collection_year": (row.get(args.year_col) or "").strip() or iso_date(row, args.year_col, args.month_col, args.day_col, args.date_col)[:4],
               "collection_date": iso_date(row, args.year_col, args.month_col, args.day_col, args.date_col),
               "ward": (row.get(args.ward_col) or "").strip() or "unknown",
               "species_phenotypic": organism, "species": species}
        out.update({k: v for k, v in row.items() if k != args.id_col})
        meta_out.append(out)

    # Contigs files with no metadata row
    ids = {r["sample"] for r in sheet} | {s for s, _ in skipped} | {o["sample"] for o in others} | not_selected
    orphans = sorted(f for f in os.listdir(folder)
                     if re.search(r"\.(fa|fasta|fna)(\.gz)?$", f)
                     and not any(f.startswith(i + suffix) or f == i + suffix
                                 for i in ids for suffix in (args.suffix, ".fasta", ".fa", ".fna", "_contigs.fa", "_contigs.fna")))

    out_dir = os.path.dirname(os.path.abspath(args.out))
    os.makedirs(out_dir, exist_ok=True)
    stem, ext = os.path.splitext(args.out)
    meta_path = args.metadata_out or os.path.join(out_dir, "metadata_assemblies.tsv")
    meta_stem, meta_ext = os.path.splitext(meta_path)

    def write_sheet(path, rows_):
        with open(path, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["sample", "species", "fasta", "phylo_group"])
            writer.writeheader()
            writer.writerows(rows_)

    def write_meta(path, rows_):
        with open(path, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows_[0].keys()), delimiter="\t", extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows_)

    written = []
    if args.split_by_species:
        for sp in ("kpneumoniae", "ecoli"):
            rows_sp = [r for r in sheet if r["species"] == sp]
            if not rows_sp:
                continue
            ids_sp = {r["sample"] for r in rows_sp}
            write_sheet(f"{stem}_{sp}{ext or '.csv'}", rows_sp)
            write_meta(f"{meta_stem}_{sp}{meta_ext or '.tsv'}", [m for m in meta_out if m["sample"] in ids_sp])
            written.append((f"{stem}_{sp}{ext or '.csv'}", f"{meta_stem}_{sp}{meta_ext or '.tsv'}", len(rows_sp), sp))
    elif sheet:
        write_sheet(args.out, sheet)
        write_meta(meta_path, meta_out)
        by_species = {}
        for r in sheet:
            by_species[r["species"]] = by_species.get(r["species"], 0) + 1
        written.append((args.out, meta_path, len(sheet), ", ".join(f"{n} {s}" for s, n in sorted(by_species.items()))))
    if others:
        other_path = f"{stem}_other_species.tsv"
        with open(other_path, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["sample", "organism", "fasta"], delimiter="\t")
            writer.writeheader()
            writer.writerows(others)

    for sheet_path, meta_file, n, label in written:
        print(f"Samplesheet: {sheet_path}  ({n} samples: {label})")
        print(f"  metadata for the R scripts: {meta_file}  (sample, collection_year, collection_date, ward = {args.ward_col}, ...)")
    if others:
        counts = {}
        for o in others:
            counts[o["organism"]] = counts.get(o["organism"], 0) + 1
        print(f"Other species (not typed by the pipeline): {stem}_other_species.tsv  ("
              + ", ".join(f"{n} {o}" for o, n in sorted(counts.items())) + ")")
    if skipped:
        print(f"\nLeft out ({len(skipped)}):")
        for sample, why in skipped:
            print(f"  {sample}: {why}")
    if orphans:
        print(f"\nContigs files without a metadata row ({len(orphans)}): {', '.join(orphans[:20])}"
              + (" ..." if len(orphans) > 20 else ""))
    if not written:
        sys.exit("ERROR: no K. pneumoniae or E. coli samples written" + (f" for --species {args.species}" if args.species != "all" else ""))


if __name__ == "__main__":
    main()
