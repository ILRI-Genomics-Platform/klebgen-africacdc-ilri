#!/usr/bin/env python3
"""
filter_terra_table.py - keep only the rows of a Terra data table (TSV) that match a specimen type, a
maximum age and, optionally, an organism, e.g. neonatal K. pneumoniae bloodstream isolates. The output is a
TSV with the same columns, ready for bin/download_terra_outputs.py.

    python3 bin/filter_terra_table.py --terra-output KlebGenNov2025.tsv --out KlebGen_blood_neonates.tsv
        # defaults: Specimen_type is 'blood' (any case) and Age <= 30 days; every organism

    python3 bin/filter_terra_table.py --terra-output KlebGenNov2025.tsv --out KlebGen_kpn_blood_neonates.tsv \
        --organism "K. pneumoniae"
        # only Organism_detected_phenotypically = K. pneumoniae (also 'E. coli', or both: --organism kpn ecoli)

    python3 bin/download_terra_outputs.py --terra-output KlebGen_blood_neonates.tsv --outdir ... --by-country

Ages are read with their unit: '30D', '30 d', '30 days', '2W', '3M' (x 30 days), '1Y' (x 365.25), '12H'
(hours). A '<' or '>' before the number is ignored. A number without a unit is ambiguous (days? months?
years?), so it is NOT kept unless --unitless-unit says what it means. Every row left out is counted by reason,
and the unreadable ages and the other specimen types and organisms are listed, so nothing is dropped silently.

Organisms are matched as in bin/make_assembly_samplesheet.py: 'K. pneumoniae', 'Klebsiella pneumoniae',
'kpneumoniae' (and K. variicola / quasipneumoniae, the rest of the K. pneumoniae species complex) count as
K. pneumoniae; 'E. coli', 'Escherichia coli', 'ecoli' as E. coli.
Standard library only.
"""

import argparse
import csv
import re
import sys
from collections import Counter

UNIT_DAYS = {"h": 1 / 24, "hr": 1 / 24, "hrs": 1 / 24, "hour": 1 / 24, "hours": 1 / 24,
             "d": 1, "day": 1, "days": 1, "w": 7, "wk": 7, "wks": 7, "week": 7, "weeks": 7,
             "m": 30, "mo": 30, "mon": 30, "month": 30, "months": 30,
             "y": 365.25, "yr": 365.25, "yrs": 365.25, "year": 365.25, "years": 365.25}
# Ontology accessions that turn up in Terra country columns (same rules as download_terra_outputs.py)
ONTOLOGY_ID = re.compile(r"(?<![A-Za-z0-9])(?:GAZ|NCBITaxon|UBERON|OBI|DOID|MONDO|PATO|EFO|BTO|ENVO)"
                         r"[:_]\d+(?![A-Za-z0-9])", re.IGNORECASE)


def clean_country(value):
    """'Ethiopia GAZ_00000567' -> 'Ethiopia'; a value that is only an ID is returned unchanged."""
    raw = (value or "").strip()
    text = ONTOLOGY_ID.sub(" ", raw)
    text = re.sub(r"[\[\](){}<>,;:|/]+", " ", text)
    text = re.sub(r"\s+", " ", text).strip(" _-.")
    return text or raw


def species_code(value):
    """kpneumoniae, ecoli, or None for other organisms (same rules as make_assembly_samplesheet.py)."""
    s = re.sub(r"[^a-z]", "", (value or "").lower())
    if s in ("kpneumoniae", "klebsiellapneumoniae", "kpn", "kp") or "pneumoniae" in s or "variicola" in s:
        return "kpneumoniae"
    if s in ("ecoli", "escherichiacoli", "eco"):
        return "ecoli"
    return None


AGE = re.compile(r"^[<>~]?\s*(\d+(?:\.\d+)?)\s*([a-z]*)\.?$")


def age_in_days(value, unitless_unit):
    """Age in days, or None if it cannot be read."""
    m = AGE.match(value.strip().lower())
    if not m:
        return None
    number, unit = float(m.group(1)), m.group(2) or unitless_unit
    return number * UNIT_DAYS[unit] if unit in UNIT_DAYS else None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--terra-output", required=True, help="Terra data table (TSV)")
    ap.add_argument("--out", required=True, help="filtered table to write (TSV)")
    ap.add_argument("--specimen-col", default="Specimen_type", help="specimen column [Specimen_type]")
    ap.add_argument("--specimen", nargs="+", default=["blood"],
                    help="specimen types to keep, any case [blood]")
    ap.add_argument("--specimen-contains", action="store_true",
                    help="keep specimens that CONTAIN the word (e.g. 'Blood culture'), not only exact matches")
    ap.add_argument("--organism-col", default="Organism_detected_phenotypically",
                    help="organism column [Organism_detected_phenotypically]")
    ap.add_argument("--organism", nargs="+",
                    help="organisms to keep: 'K. pneumoniae' and/or 'E. coli' (any spelling) [all]")
    ap.add_argument("--country-col", default="Country", help="country column, with --countries [Country]")
    ap.add_argument("--countries", nargs="+",
                    help="keep only these countries; ontology IDs in the column are ignored, so "
                         "'--countries Ethiopia' also matches 'Ethiopia GAZ_00000567'")
    ap.add_argument("--age-col", default="Age", help="age column [Age]")
    ap.add_argument("--max-age-days", type=float, default=30, help="keep ages <= this many days [30]")
    ap.add_argument("--unitless-unit", choices=["d", "w", "m", "y"],
                    help="unit for ages given as a bare number; without it those rows are left out")
    args = ap.parse_args()

    wanted = [s.strip().lower() for s in args.specimen]
    organisms = None
    if args.organism:
        organisms = {species_code(o) for o in args.organism}
        if None in organisms:
            sys.exit(f"--organism: use 'K. pneumoniae' or 'E. coli' (got {', '.join(args.organism)})")
    with open(args.terra_output, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        needed = [args.specimen_col, args.age_col] + ([args.organism_col] if organisms else []) \
            + ([args.country_col] if args.countries else [])
        for col in needed:
            if col not in header:
                sys.exit(f"No '{col}' column in {args.terra_output}. Columns: {', '.join(header)}")
        s_idx, a_idx = header.index(args.specimen_col), header.index(args.age_col)
        o_idx = header.index(args.organism_col) if organisms else None
        c_idx = header.index(args.country_col) if args.countries else None
        countries = {c.strip().lower() for c in args.countries} if args.countries else None
        kept, reasons, other_specimens, bad_ages, other_orgs = [], Counter(), Counter(), Counter(), Counter()
        other_countries = Counter()
        for row in reader:
            if not row or not row[0].strip():
                continue
            row += [""] * (len(header) - len(row))
            if countries is not None:
                raw = row[c_idx].strip()
                if clean_country(raw).lower() not in countries and raw.lower() not in countries:
                    reasons["other country"] += 1
                    other_countries[raw or "(empty)"] += 1
                    continue
            specimen = row[s_idx].strip().lower()
            ok = any(w in specimen for w in wanted) if args.specimen_contains else specimen in wanted
            if not ok:
                reasons["other specimen type"] += 1
                other_specimens[row[s_idx].strip() or "(empty)"] += 1
                continue
            if organisms and species_code(row[o_idx]) not in organisms:
                reasons["other organism"] += 1
                other_orgs[row[o_idx].strip() or "(empty)"] += 1
                continue
            days = age_in_days(row[a_idx], args.unitless_unit)
            if days is None:
                reasons["age missing or unreadable"] += 1
                bad_ages[row[a_idx].strip() or "(empty)"] += 1
                continue
            if days > args.max_age_days:
                reasons[f"older than {args.max_age_days:g} days"] += 1
                continue
            kept.append(row)

    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(kept)

    total = len(kept) + sum(reasons.values())
    print(f"{len(kept)} of {total} rows kept -> {args.out}")
    for reason, n in reasons.most_common():
        print(f"  left out, {reason}: {n}")
    if other_countries:
        print("  other countries: " + ", ".join(f"{k} ({v})" for k, v in other_countries.most_common(15)))
    if other_specimens:
        print("  other specimen types: " + ", ".join(f"{k} ({v})" for k, v in other_specimens.most_common(15)))
    if other_orgs:
        print("  other organisms: " + ", ".join(f"{k} ({v})" for k, v in other_orgs.most_common(15)))
    if bad_ages:
        print("  ages not read: " + ", ".join(f"{k!r} ({v})" for k, v in bad_ages.most_common(15)))
        if any(re.fullmatch(r"\d+(\.\d+)?", k) for k in bad_ages):
            print("  -> bare numbers have no unit: add --unitless-unit d (or w, m, y) if you know it")


if __name__ == "__main__":
    main()
