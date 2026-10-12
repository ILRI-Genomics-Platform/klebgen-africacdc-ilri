#!/usr/bin/env python3
"""
fetch_ena.py - pull a chosen list of isolates from ENA and build the
neonatal-amr-nf samplesheet
(sample,species,fastq_1,fastq_2,long_fastq,phylo_group).

Standard library only (Python 3.8+), so it runs on any training server.

Two subcommands
---------------
1) list   - inspect a BioProject before choosing isolates
     python3 bin/fetch_ena.py list PRJEB42462 -o PRJEB42462_runs.tsv

   Writes one row per run with sample metadata (alias, date, country,
   platform, layout) so trainers can pick isolates in a spreadsheet.

2) fetch  - download the chosen isolates and write samplesheet.csv
     python3 bin/fetch_ena.py fetch accessions.tsv -d data -o samplesheet.csv

   accessions.tsv (tab-separated, header required; phylo_group optional):
     id        species       accession            phylo_group
     MW_ST39_01 kpneumoniae  ERS1234567
     MW_ST39_02 kpneumoniae  ERR111111;ERR222222
     BN_EC_01   ecoli        SAMEA7654321

   - accession: any ENA run (ERR/SRR/DRR), experiment (ERX/SRX), sample
     (ERS/SRS/SAMEA/SAMN/SAMD) accession. Separate several with ';' - e.g.
     an isolate whose Illumina and ONT runs are listed separately.
   - species: kpneumoniae | ecoli (as the workflow expects).
   Illumina paired-end runs become fastq_1/fastq_2; Oxford Nanopore runs
   become long_fastq; phylo_group is copied through (isolates sharing a
   label are analysed together in the phylogeny step).
   Several runs of the same type for one isolate are concatenated.

Useful options
   --dry-run      resolve accessions and write the sheet + a URL list,
                  but download nothing (check before pulling ~GBs)
   --no-ont       ignore ONT runs (short-read only training)
   --threads N    parallel downloads (default 4)

Every download is checked against ENA's MD5. Re-running skips files that
are already present and verified. A manifest (runs_manifest.tsv) records
which run accessions were used, for citation in course materials.
"""

import argparse
import concurrent.futures as cf
import csv
import hashlib
import os
import shutil
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

ENA_PORTAL = "https://www.ebi.ac.uk/ena/portal/api/filereport"
RUN_FIELDS = [
    "run_accession", "sample_accession", "secondary_sample_accession",
    "sample_alias", "scientific_name", "instrument_platform",
    "instrument_model", "library_layout", "library_strategy",
    "read_count", "base_count", "collection_date", "country",
    "fastq_ftp", "fastq_md5", "fastq_bytes",
]
VALID_SPECIES = {"kpneumoniae", "ecoli"}


# ---------------------------------------------------------------- ENA API --

def http_get(url, retries=5, timeout=120):
    """GET with polite retries/backoff (ENA rate-limits bursts)."""
    delay = 5
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "neonatal-amr-nf/fetch_ena"})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read().decode("utf-8")
        except urllib.error.HTTPError as e:
            if e.code == 204:          # ENA: no content for this accession
                return ""
            if e.code in (429, 500, 502, 503, 504) and attempt < retries:
                time.sleep(delay)
                delay *= 2
                continue
            raise
        except urllib.error.URLError:
            if attempt < retries:
                time.sleep(delay)
                delay *= 2
                continue
            raise
    return ""


def ena_runs(accession):
    """Return a list of run dicts for any ENA accession (project, sample, run...)."""
    query = urllib.parse.urlencode({
        "accession": accession,
        "result": "read_run",
        "fields": ",".join(RUN_FIELDS),
        "format": "tsv",
        "limit": 0,
    })
    text = http_get(f"{ENA_PORTAL}?{query}")
    if not text.strip():
        return []
    return list(csv.DictReader(text.splitlines(), delimiter="\t"))


# ------------------------------------------------------- run classification --

def to_https(path):
    """ENA returns 'ftp.sra.ebi.ac.uk/vol1/...' - use HTTPS (firewall-friendly)."""
    return "https://" + path if not path.startswith(("http://", "https://")) else path


def run_files(run):
    """Pair each FASTQ path with its MD5 and size."""
    paths = [p for p in run.get("fastq_ftp", "").split(";") if p]
    md5s = run.get("fastq_md5", "").split(";")
    sizes = run.get("fastq_bytes", "").split(";")
    out = []
    for i, p in enumerate(paths):
        out.append({
            "url": to_https(p),
            "name": os.path.basename(p),
            "md5": md5s[i] if i < len(md5s) else "",
            "bytes": int(sizes[i]) if i < len(sizes) and sizes[i].isdigit() else 0,
        })
    return out


def classify(run):
    """Return ('illumina', (f1, f2)) | ('ont', [files]) | (None, reason)."""
    platform = run.get("instrument_platform", "").upper()
    files = run_files(run)
    if not files:
        return None, "no FASTQ files"
    if platform == "ILLUMINA":
        r1 = [f for f in files if f["name"].endswith(("_1.fastq.gz", "_1.fq.gz"))]
        r2 = [f for f in files if f["name"].endswith(("_2.fastq.gz", "_2.fq.gz"))]
        if len(r1) == 1 and len(r2) == 1:
            return "illumina", (r1[0], r2[0])
        return None, f"Illumina run is not paired-end ({len(files)} file(s))"
    if platform == "OXFORD_NANOPORE":
        return "ont", files
    return None, f"platform {platform or 'unknown'} not used by the workflow"


# --------------------------------------------------------------- download --

def md5sum(path, chunk=1 << 20):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def download(fileinfo, dest_dir):
    """Download one file with MD5 check; skip if already present and valid."""
    dest = os.path.join(dest_dir, fileinfo["name"])
    if os.path.exists(dest) and (not fileinfo["md5"] or md5sum(dest) == fileinfo["md5"]):
        return dest, "cached"
    tmp = dest + ".part"
    for attempt in range(1, 4):
        try:
            req = urllib.request.Request(fileinfo["url"], headers={"User-Agent": "neonatal-amr-nf/fetch_ena"})
            with urllib.request.urlopen(req, timeout=600) as resp, open(tmp, "wb") as out:
                shutil.copyfileobj(resp, out, length=1 << 20)
            if fileinfo["md5"] and md5sum(tmp) != fileinfo["md5"]:
                raise IOError("MD5 mismatch")
            os.replace(tmp, dest)
            return dest, "downloaded"
        except Exception as e:  # noqa: BLE001 - report and retry any failure
            if attempt == 3:
                if os.path.exists(tmp):
                    os.remove(tmp)
                raise RuntimeError(f"{fileinfo['url']}: {e}") from e
            time.sleep(10 * attempt)
    return dest, "failed"


def concat_gz(parts, dest):
    """Concatenated gzip members form a valid gzip file - no recompression needed."""
    if len(parts) == 1:
        return parts[0]
    with open(dest, "wb") as out:
        for p in parts:
            with open(p, "rb") as fh:
                shutil.copyfileobj(fh, out, length=1 << 20)
    return dest


# --------------------------------------------------------------- commands --

def ensure_output_dir(path):
    """Create the folder an output file will be written to, if it is missing.
    Returns that folder's absolute path."""
    out_dir = os.path.dirname(os.path.abspath(path))
    if os.path.exists(out_dir) and not os.path.isdir(out_dir):
        sys.exit(f"Cannot write {path}: {out_dir} exists but is a file, not a folder")
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir, exist_ok=True)
        print(f"Created output folder {out_dir}", file=sys.stderr)
    return out_dir

def cmd_list(args):
    runs = ena_runs(args.project)
    if not runs:
        sys.exit(f"No runs found for {args.project}")
    cols = ["run_accession", "sample_accession", "sample_alias", "scientific_name",
            "instrument_platform", "instrument_model", "library_layout",
            "read_count", "base_count", "collection_date", "country"]
    if args.output:
        ensure_output_dir(args.output)
    out = open(args.output, "w", newline="") if args.output else sys.stdout
    w = csv.DictWriter(out, fieldnames=cols, delimiter="\t", extrasaction="ignore")
    w.writeheader()
    for r in runs:
        w.writerow(r)
    if args.output:
        out.close()
        platforms = {}
        for r in runs:
            platforms[r.get("instrument_platform", "?")] = platforms.get(r.get("instrument_platform", "?"), 0) + 1
        summary = ", ".join(f"{k}: {v}" for k, v in sorted(platforms.items()))
        print(f"{len(runs)} runs written to {args.output} ({summary})", file=sys.stderr)


def read_accession_table(path):
    with open(path, newline="") as fh:
        lines = [ln for ln in fh if ln.strip() and not ln.lstrip().startswith("#")]
    rows = list(csv.DictReader(lines, delimiter="\t"))
    placeholders = [r["id"] for r in rows if "<" in (r.get("accession") or "")]
    if placeholders:
        sys.exit(f"{path}: replace the <placeholder> accessions for: {', '.join(placeholders)}")
    if not rows:
        sys.exit(f"{path}: no rows")
    missing = {"id", "species", "accession"} - set(rows[0])
    if missing:
        sys.exit(f"{path}: missing column(s) {', '.join(sorted(missing))} (file must be tab-separated)")
    seen = set()
    for i, r in enumerate(rows, start=2):
        r["id"] = r["id"].strip()
        r["species"] = r["species"].strip().lower()
        if not r["id"]:
            sys.exit(f"{path}: line {i}: empty id")
        if r["id"] in seen:
            sys.exit(f"{path}: line {i}: duplicate id {r['id']}")
        if r["species"] not in VALID_SPECIES:
            sys.exit(f"{path}: line {i}: species must be one of {sorted(VALID_SPECIES)}, got '{r['species']}'")
        seen.add(r["id"])
    return rows


def resolve(rows, use_ont, query=None):
    """Map each isolate to its Illumina pair(s) and ONT file(s)."""
    query = query or ena_runs
    plan, problems = [], []
    for r in rows:
        illumina, ont, runs_used = [], [], []
        for acc in [a.strip() for a in r["accession"].split(";") if a.strip()]:
            runs = query(acc)
            if not runs:
                problems.append(f"{r['id']}: no runs found for {acc}")
                continue
            for run in runs:
                kind, payload = classify(run)
                if kind == "illumina":
                    illumina.append(payload)
                    runs_used.append((run, "illumina"))
                elif kind == "ont" and use_ont:
                    ont.extend(payload)
                    runs_used.append((run, "ont"))
                elif kind is None:
                    problems.append(f"{r['id']}: skipped {run.get('run_accession')} ({payload})")
        if not illumina:
            problems.append(f"{r['id']}: EXCLUDED - no paired-end Illumina run (the workflow needs r1/r2)")
            continue
        plan.append({"id": r["id"], "species": r["species"], "illumina": illumina,
                     "phylo_group": (r.get("phylo_group") or "").strip(),
                     "ont": ont, "runs": runs_used})
    return plan, problems


def cmd_fetch(args):
    rows = read_accession_table(args.accessions)
    print(f"Resolving {len(rows)} isolate(s) against ENA...", file=sys.stderr)
    plan, problems = resolve(rows, use_ont=not args.no_ont)
    for p in problems:
        print("WARNING:", p, file=sys.stderr)
    if not plan:
        sys.exit("Nothing to fetch.")

    data_dir = os.path.abspath(args.data_dir)
    raw_dir = os.path.join(data_dir, "raw")
    os.makedirs(raw_dir, exist_ok=True)

    all_files = {}
    for iso in plan:
        for f1, f2 in iso["illumina"]:
            all_files[f1["url"]] = f1
            all_files[f2["url"]] = f2
        for f in iso["ont"]:
            all_files[f["url"]] = f
    total_gb = sum(f["bytes"] for f in all_files.values()) / 1e9
    n_ont = sum(1 for i in plan if i["ont"])
    print(f"{len(plan)} isolate(s), {n_ont} with ONT; {len(all_files)} file(s), {total_gb:.1f} GB",
          file=sys.stderr)

    # Manifest of runs used - for citing ENA accessions in course materials.
    # The sample sheet, manifest and URL list go in the sample sheet's folder.
    out_dir = ensure_output_dir(args.output)
    manifest = os.path.join(out_dir, "runs_manifest.tsv")
    with open(manifest, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["id", "species", "type", "run_accession", "sample_accession",
                    "sample_alias", "instrument_model", "collection_date", "country"])
        for iso in plan:
            for run, kind in iso["runs"]:
                w.writerow([iso["id"], iso["species"], kind, run.get("run_accession", ""),
                            run.get("sample_accession", ""), run.get("sample_alias", ""),
                            run.get("instrument_model", ""), run.get("collection_date", ""),
                            run.get("country", "")])

    if args.dry_run:
        urls = os.path.join(out_dir, "download_urls.txt")
        with open(urls, "w") as fh:
            fh.write("\n".join(sorted(all_files)) + "\n")
        local = {u: os.path.join(raw_dir, f["name"]) for u, f in all_files.items()}
    else:
        local = {}
        with cf.ThreadPoolExecutor(max_workers=args.threads) as pool:
            futures = {pool.submit(download, f, raw_dir): u for u, f in all_files.items()}
            for i, fut in enumerate(cf.as_completed(futures), start=1):
                url = futures[fut]
                path, status = fut.result()
                local[url] = path
                print(f"[{i}/{len(futures)}] {status}: {os.path.basename(path)}", file=sys.stderr)

    # Build per-isolate files (concatenate multiple runs) and the sample sheet
    sheet_rows = []
    for iso in plan:
        r1_parts = [local[f1["url"]] for f1, _ in iso["illumina"]]
        r2_parts = [local[f2["url"]] for _, f2 in iso["illumina"]]
        ont_parts = [local[f["url"]] for f in iso["ont"]]
        if args.dry_run or len(r1_parts) == 1:
            r1 = r1_parts[0] if len(r1_parts) == 1 else os.path.join(data_dir, f"{iso['id']}_R1.fastq.gz")
            r2 = r2_parts[0] if len(r2_parts) == 1 else os.path.join(data_dir, f"{iso['id']}_R2.fastq.gz")
        else:
            r1 = concat_gz(r1_parts, os.path.join(data_dir, f"{iso['id']}_R1.fastq.gz"))
            r2 = concat_gz(r2_parts, os.path.join(data_dir, f"{iso['id']}_R2.fastq.gz"))
        ont = ""
        if ont_parts:
            if len(ont_parts) == 1:
                ont = ont_parts[0]
            else:
                ont = os.path.join(data_dir, f"{iso['id']}_ont.fastq.gz")
                if not args.dry_run:
                    concat_gz(ont_parts, ont)
        sheet_rows.append({"sample": iso["id"], "species": iso["species"], "fastq_1": r1, "fastq_2": r2,
                           "long_fastq": ont, "phylo_group": iso["phylo_group"]})

    with open(args.output, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["sample", "species", "fastq_1", "fastq_2", "long_fastq", "phylo_group"])
        w.writeheader()
        w.writerows(sheet_rows)

    print(f"Sample sheet: {args.output} ({len(sheet_rows)} isolates)", file=sys.stderr)
    print(f"Run manifest: {manifest}", file=sys.stderr)
    if args.dry_run:
        print(f"Dry run: nothing downloaded; URL list in {os.path.join(out_dir, 'download_urls.txt')}", file=sys.stderr)
    if problems:
        print(f"{len(problems)} warning(s) above - review before running the workflow.", file=sys.stderr)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_list = sub.add_parser("list", help="list runs and metadata for a BioProject")
    p_list.add_argument("project", help="e.g. PRJEB42462")
    p_list.add_argument("-o", "--output", help="TSV file (default: stdout); its folder is created if missing")
    p_list.set_defaults(func=cmd_list)

    p_fetch = sub.add_parser("fetch", help="download chosen isolates and write the sample sheet")
    p_fetch.add_argument("accessions", help="TSV with columns id, species, accession")
    p_fetch.add_argument("-d", "--data-dir", default="data", help="download directory (default: data)")
    p_fetch.add_argument("-o", "--output", default="samplesheet.csv",
                         help="sample sheet path; its folder is created if missing and also receives "
                              "runs_manifest.tsv (and download_urls.txt with --dry-run)")
    p_fetch.add_argument("--threads", type=int, default=4, help="parallel downloads (default 4)")
    p_fetch.add_argument("--no-ont", action="store_true", help="ignore Oxford Nanopore runs")
    p_fetch.add_argument("--dry-run", action="store_true", help="resolve and write sheet, no downloads")
    p_fetch.set_defaults(func=cmd_fetch)

    args = ap.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
