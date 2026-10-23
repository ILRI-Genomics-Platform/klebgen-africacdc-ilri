#!/usr/bin/env python3
"""
download_terra_outputs.py - download selected output files (assemblies, Kleborate reports, MLST) from a
Terra data table, in parallel, ready for the pipeline's assembly route.

    python3 bin/download_terra_outputs.py --terra-output KlebGenNov2025.tsv --outdir results/africa
        # assemblies only (the default): results/africa/assemblies/<sample>_contigs.fasta

    python3 bin/download_terra_outputs.py --terra-output KlebGenNov2025.tsv --outdir results/africa \\
        --files assembly_fasta kleborate_klebsiella_output_file kleborate_ecoli_output_file ts_mlst_results \\
        --by-country --jobs 16

The Terra table is the TSV downloaded from the workspace's Data tab: its first column holds the sample ID
(e.g. 'entity:sample_id'), and each output column holds a gs:// path. Standard library only (no pandas);
needs the Google Cloud CLI ('gcloud', logged in with 'gcloud auth login').

Faster than one 'gcloud storage cp' per file:
  - one pass over the table for every file type requested;
  - files already downloaded (non-empty) are skipped, so a rerun only fetches what is missing;
  - when the files of one type have distinct names in the bucket (the usual case), ONE 'gcloud storage cp -I'
    call downloads them all, with gcloud's own parallel transfers; otherwise --jobs downloads at a time;
  - each file is written to <name>.part and renamed when complete, so an interrupted run never leaves a
    truncated file that a rerun would skip.

Output names match bin/make_assembly_samplesheet.py: assemblies/<sample>_contigs.fasta (gzipped assemblies are
decompressed, as the pipeline needs plain FASTA). With --by-country, files go to <outdir>/<Country>/<type>/.

Country folders are named after the country alone: ontology IDs in the country column ('Ethiopia
GAZ_00000567', 'Kenya [GAZ:00001101]') are dropped, so the folder is 'Ethiopia', not 'Ethiopia_GAZ_00000567_'.
A value that is nothing but an ID has no name to keep, so it is used as it is and reported at the end;
give those a name with --country-map (a TSV of 'raw<TAB>country'), or keep the IDs with --keep-country-ids.

One log per country is written as <folder>/terra_downloads_<Country>.tsv (in that country's folder with
--by-country, otherwise in --outdir), with the columns sample, country, country_raw, type, source,
destination and status. --log-prefix renames the file.
"""

import argparse
import csv
import gzip
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

# column in the Terra table -> (sub-folder, file name suffix)
KINDS = {
    "assembly_fasta": ("assemblies", "_contigs.fasta"),
    "kleborate_klebsiella_output_file": ("kleborate-klebsiella", "_kleborate_kpsc_out.tsv"),
    "kleborate_ecoli_output_file": ("kleborate-ecoli", "_escherichia_output.txt"),
    "ts_mlst_results": ("ts-mlst", "_ts_mlst.tsv"),
}
SAFE_ID = re.compile(r"^[A-Za-z0-9_.-]+$")
# Ontology accessions that turn up in Terra metadata columns, e.g. GAZ_00001093, GAZ:00000567
# '_' is a word character, so \b would miss 'Ethiopia_GAZ_00000567': match on non-alphanumeric instead
ONTOLOGY_ID = re.compile(r"(?<![A-Za-z0-9])(?:GAZ|NCBITaxon|UBERON|OBI|DOID|MONDO|PATO|EFO|BTO|ENVO)"
                         r"[:_]\d+(?![A-Za-z0-9])", re.IGNORECASE)


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def clean_country(value, mapping=None, keep_ids=False):
    """Country name with ontology IDs removed: 'Ethiopia GAZ_00000567' -> 'Ethiopia'.

    A value that is only an ID ('GAZ_00001093') has no name to keep, so it is returned unchanged rather
    than dropped: nothing is silently lost, and main() reports those values. --country-map names them.
    """
    raw = (value or "").strip()
    if mapping and raw in mapping:
        return mapping[raw]
    if keep_ids:
        return raw
    text = ONTOLOGY_ID.sub(" ", raw)
    text = re.sub(r"[\[\](){}<>,;:|/]+", " ", text)      # brackets or separators the ID sat in
    text = re.sub(r"\s+", " ", text).strip(" _-.")
    return text or raw


def country_dir_name(country):
    """Folder name for a country: letters, digits, _ . - only, no leading or trailing punctuation."""
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", country).strip("_.-") or "unknown"


def read_country_map(path):
    """TSV (or CSV) of 'raw value<TAB>country name'; lines starting with '#' are comments."""
    mapping = {}
    with open(path, newline="") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            parts = line.split("\t") if "\t" in line else line.split(",")
            if len(parts) < 2:
                sys.exit(f"{path}: need two columns 'raw<TAB>country', got {line!r}")
            mapping[parts[0].strip()] = parts[1].strip()
    if not mapping:
        sys.exit(f"{path} has no rows")
    return mapping


def read_table(path, kinds, country_col, countries, samples, by_country, cmap=None, keep_ids=False):
    """Rows to download: dicts with sample, country, kind, src. Also counts of what was skipped."""
    with open(path, newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        missing = [k for k in kinds if k not in header]
        if missing:
            sys.exit(f"Column(s) not in {path}: {', '.join(missing)}\nColumns: {', '.join(header)}")
        idx = {name: i for i, name in enumerate(header)}
        c_idx = idx.get(country_col)
        if c_idx is None and (countries or by_country):
            sys.exit(f"No '{country_col}' column in {path} (set --country-col)")
        jobs, empty, bad_ids, seen = [], Counter(), set(), set()
        for row in reader:
            if not row or not row[0].strip():
                continue
            sample = row[0].strip()
            raw = (row[c_idx].strip() if c_idx is not None and c_idx < len(row) else "")
            country = clean_country(raw, cmap, keep_ids) or "unknown"
            # --countries matches the cleaned name or the raw value, so 'Ethiopia' also selects
            # rows whose country column reads 'Ethiopia GAZ_00000567'
            if countries and country.lower() not in countries and raw.lower() not in countries:
                continue
            if samples and sample not in samples:
                continue
            if sample in seen:
                log(f"WARNING: sample {sample} appears more than once in the table; using its first row")
                continue
            seen.add(sample)
            if not SAFE_ID.match(sample):
                bad_ids.add(sample)
            for kind in kinds:
                src = row[idx[kind]].strip() if idx[kind] < len(row) else ""
                if not src.startswith("gs://"):
                    empty[kind] += 1
                    continue
                jobs.append({"sample": sample, "country": country, "country_raw": raw or "unknown",
                             "kind": kind, "src": src})
    return jobs, empty, bad_ids


def destination(job, outdir, by_country):
    folder, suffix = KINDS[job["kind"]]
    country_dir = country_dir_name(job["country"]) if by_country else ""
    return os.path.join(outdir, country_dir, folder, job["sample"] + suffix)


def gcloud_base(args):
    cmd = ["gcloud", "storage", "cp"]
    if args.billing_project:
        cmd.append(f"--billing-project={args.billing_project}")
    return cmd


def finish(tmp_path, job):
    """Move a downloaded file into place (decompressing gzipped assemblies); returns a status."""
    dest = job["dest"]
    if not os.path.isfile(tmp_path) or os.path.getsize(tmp_path) == 0:
        return "failed: empty or missing after download"
    part = dest + ".part"
    if job["kind"] == "assembly_fasta" and job["src"].endswith(".gz"):
        with gzip.open(tmp_path, "rb") as fin, open(part, "wb") as fout:
            shutil.copyfileobj(fin, fout, 1 << 20)
        os.remove(tmp_path)
    else:
        shutil.move(tmp_path, part)
    if job["kind"] == "assembly_fasta":
        with open(part, "rb") as fh:
            if fh.read(1) != b">":
                os.remove(part)
                return "failed: not a FASTA file"
    os.replace(part, dest)
    return "downloaded"


def gcloud_error(stderr):
    """The informative part of a gcloud error (gcloud wraps long messages over several lines)."""
    text = " ".join(stderr.split())
    m = re.search(r"ERROR:.*", text)
    return (m.group(0) if m else text)[:300] or "gcloud error"


def check_login(args):
    """Stop early, with a clear message, if gcloud has no active account."""
    r = subprocess.run(["gcloud", "auth", "list", "--filter=status:ACTIVE", "--format=value(account)"],
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    account = r.stdout.strip()
    if r.returncode != 0 or not account:
        sys.exit("gcloud has no active account on this machine.\n"
                 "  Log in (on a login node):   gcloud auth login --no-launch-browser\n"
                 "  or pick an account:         gcloud auth list ; gcloud config set account <ACCOUNT>")
    log(f"gcloud account: {account}")


def download_one(job, args):
    """Fallback: one gcloud call for one file."""
    tmp = job["dest"] + ".download"
    cmd = gcloud_base(args) + ["--no-user-output-enabled", job["src"], tmp]
    for attempt in range(1, args.retries + 2):
        r = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
        if r.returncode == 0:
            return finish(tmp, job)
        if re.search(r"not found|No URLs matched|404|403|AccessDenied|active account|auth login|"
                     r"reauthenticat", r.stderr, re.I):
            break                                  # permanent: retrying will not help
        time.sleep(2 * attempt)
    if os.path.exists(tmp):
        os.remove(tmp)
    return "failed: " + gcloud_error(r.stderr)


def download_batch(kind_jobs, args, outdir):
    """One 'gcloud storage cp -I' for many files with distinct bucket names; gcloud runs the transfers in
    parallel. Returns {dest: status}; files it did not deliver are retried one by one."""
    staging = tempfile.mkdtemp(prefix=".terra_staging_", dir=outdir)
    cmd = gcloud_base(args) + ["-I", staging + "/"]
    r = subprocess.run(cmd, input="\n".join(j["src"] for j in kind_jobs) + "\n",
                       stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
    status = {}
    for job in kind_jobs:
        got = os.path.join(staging, os.path.basename(job["src"]))
        if os.path.isfile(got):
            status[job["dest"]] = finish(got, job)
    shutil.rmtree(staging, ignore_errors=True)
    if r.returncode != 0:
        n = len(kind_jobs) - len(status)
        log(f"  gcloud reported errors ({n} file(s) not delivered): {gcloud_error(r.stderr)}")
        log("  retrying those one at a time")
    return status


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--terra-output", required=True, help="Terra data table (TSV, first column = sample ID)")
    ap.add_argument("--outdir", "--outDir", default=".", help="output folder [.]")
    ap.add_argument("--files", "--files-to-download", nargs="+", default=["assembly_fasta"], choices=list(KINDS),
                    help="output columns to download [assembly_fasta]")
    ap.add_argument("--jobs", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 8)),
                    help="parallel downloads in the one-file-at-a-time fallback [$SLURM_CPUS_PER_TASK or 8]")
    ap.add_argument("--by-country", action="store_true", help="put files in <outdir>/<Country>/<type>/")
    ap.add_argument("--country-col", default="Country", help="country column [Country]")
    ap.add_argument("--countries", nargs="+",
                    help="only these countries (case-insensitive; matches the cleaned name or the raw value)")
    ap.add_argument("--country-map", help="TSV 'raw value<TAB>country name', for country values that are "
                                          "only an ontology ID (e.g. GAZ_00001093)")
    ap.add_argument("--keep-country-ids", action="store_true",
                    help="keep ontology IDs in country folder names (the old behaviour)")
    ap.add_argument("--log-prefix", default="terra_downloads",
                    help="stem of the per-country log file [terra_downloads]")
    ap.add_argument("--samples", help="only the sample IDs listed in this file (one per line)")
    ap.add_argument("--billing-project", help="Google Cloud project to bill (requester-pays buckets)")
    ap.add_argument("--retries", type=int, default=2, help="retries per file in the fallback [2]")
    ap.add_argument("--dry-run", action="store_true", help="list what would be downloaded, download nothing")
    args = ap.parse_args()

    if not args.dry_run:
        if shutil.which("gcloud") is None:
            sys.exit("gcloud not found: activate the environment with the Google Cloud CLI first")
        check_login(args)
    outdir = os.path.abspath(args.outdir)
    samples = None
    if args.samples:
        with open(args.samples) as fh:
            samples = {line.strip() for line in fh if line.strip()}
    countries = {c.lower() for c in args.countries} if args.countries else None
    cmap = read_country_map(args.country_map) if args.country_map else None

    jobs, empty, bad_ids = read_table(args.terra_output, args.files, args.country_col, countries, samples,
                                      args.by_country, cmap, args.keep_country_ids)
    for kind, n in empty.items():
        log(f"{kind}: {n} sample(s) have no file in the table (skipped)")
    # country values that are only an ontology ID: no name to use, so they keep the ID
    id_only = sorted({j["country"] for j in jobs if ONTOLOGY_ID.fullmatch(j["country"].strip())})
    if id_only:
        log(f"WARNING: {len(id_only)} country value(s) are only an ontology ID, so the folder keeps the ID: "
            f"{', '.join(id_only[:5])}{' ...' if len(id_only) > 5 else ''}. "
            f"Give them names with --country-map.")
    if bad_ids:
        log(f"WARNING: {len(bad_ids)} sample ID(s) contain characters the pipeline rejects "
            f"(allowed: letters, digits, _ . -), e.g. {sorted(bad_ids)[0]!r}")

    todo, skipped = [], 0
    for job in jobs:
        job["dest"] = destination(job, outdir, args.by_country)
        if os.path.isfile(job["dest"]) and os.path.getsize(job["dest"]) > 0:
            job["status"] = "already present"
            skipped += 1
        else:
            todo.append(job)
    log(f"{len(jobs)} file(s) listed: {skipped} already present, {len(todo)} to download")
    if args.dry_run:
        for job in todo:
            print(f"{job['src']}\t{job['dest']}")
        return
    for job in todo:
        os.makedirs(os.path.dirname(job["dest"]), exist_ok=True)

    start = time.time()
    by_kind = defaultdict(list)
    for job in todo:
        by_kind[job["kind"]].append(job)
    single = []
    for kind, kind_jobs in by_kind.items():
        names = Counter(os.path.basename(j["src"]) for j in kind_jobs)
        if len(kind_jobs) > 1 and max(names.values()) == 1:
            log(f"{kind}: {len(kind_jobs)} file(s) in one batch")
            status = download_batch(kind_jobs, args, outdir)
            for job in kind_jobs:
                if job["dest"] in status and status[job["dest"]] == "downloaded":
                    job["status"] = "downloaded"
                else:
                    single.append(job)
        else:
            single.extend(kind_jobs)
    if single:
        log(f"{len(single)} file(s) one at a time, {args.jobs} in parallel")
        with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
            futures = {pool.submit(download_one, job, args): job for job in single}
            for i, fut in enumerate(as_completed(futures), 1):
                job = futures[fut]
                job["status"] = fut.result()
                if not job["status"].startswith("downloaded") or i % 50 == 0 or i == len(single):
                    log(f"  {i}/{len(single)} {job['sample']} {job['kind']}: {job['status']}")

    # One log per country, beside that country's files (with --by-country) or in --outdir.
    # A table with no country column puts everything in one '<prefix>.tsv'.
    per_country = defaultdict(list)
    for job in jobs:
        per_country[job["country"]].append(job)
    reports = []
    for country, country_jobs in sorted(per_country.items()):
        cdir = country_dir_name(country)
        folder = os.path.join(outdir, cdir) if args.by_country else outdir
        stem = args.log_prefix if country == "unknown" else f"{args.log_prefix}_{cdir}"
        report = os.path.join(folder, stem + ".tsv")
        os.makedirs(folder, exist_ok=True)
        with open(report, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["sample", "country", "country_raw", "type", "source", "destination", "status"])
            for job in country_jobs:
                w.writerow([job["sample"], job["country"], job["country_raw"], job["kind"],
                            job["src"], job["dest"], job.get("status", "")])
        reports.append(report)
    totals = Counter(j.get("status", "").split(":")[0] for j in jobs)
    log(f"Done in {time.time() - start:.0f} s: " + ", ".join(f"{v} {k}" for k, v in totals.items()) +
        f". Log(s): {', '.join(reports)}")
    if totals.get("failed"):
        sys.exit(1)


if __name__ == "__main__":
    main()
