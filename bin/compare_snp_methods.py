#!/usr/bin/env python3
"""
compare_snp_methods.py - put the SNP distances from Snippy and SKA2 side by side, for a pipeline run with
--phylo_method both (results in 06_phylogeny/<group>/ and 06_phylogeny/<group>_ska2/).

    bin/compare_snp_methods.py --phylo results/06_phylogeny
    bin/compare_snp_methods.py --phylo results/06_phylogeny --threshold 10 --out ST39_comparison

For every group with both results, each pair of isolates gets:
  snippy        SNPs between the two from Snippy's core alignment (after Gubbins, for groups of 3 or more)
  ska2          the same from SKA2's alignment (ska map onto the same reference; after Gubbins likewise)
  ska2_ref_free SKA2's reference-free distance (ska distance: all split k-mers the two share, before Gubbins)

Writes <out>_pairs.tsv (one row per pair) and <out>_summary.tsv (one row per group), and prints the summary:
how well the two methods agree (Spearman correlation of the distances, typical difference) and whether they
flag the same pairs as close (at or below --threshold SNPs, a common cut-off for possible recent transmission).
Standard library only: runs with any Python 3.
"""

import argparse
import csv
import os
import statistics
import sys


def read_matrix(path):
    """snp-dists matrix -> {(a, b): distance} for a < b, without the reference."""
    with open(path, newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    if not rows:
        return {}
    names = rows[0][1:]
    dist = {}
    for row in rows[1:]:
        a = row[0]
        for b, value in zip(names, row[1:]):
            if a == b or "Reference" in (a, b):
                continue
            key = tuple(sorted((a, b)))
            dist[key] = float(value)
    return dist


def read_ska_distance(path):
    """ska distance output (Sample1, Sample2, Distance, Mismatches) -> {(a, b): distance}."""
    dist = {}
    if not os.path.exists(path):
        return dist
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            a, b = row.get("Sample1"), row.get("Sample2")
            if a and b and row.get("Distance") not in (None, ""):
                dist[tuple(sorted((a, b)))] = float(row["Distance"])
    return dist


def ranks(values):
    order = sorted(range(len(values)), key=lambda i: values[i])
    r = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        for k in range(i, j + 1):
            r[order[k]] = (i + j) / 2 + 1
        i = j + 1
    return r


def spearman(x, y):
    if len(x) < 3:
        return None
    rx, ry = ranks(x), ranks(y)
    mx, my = statistics.fmean(rx), statistics.fmean(ry)
    sxy = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    sxx = sum((a - mx) ** 2 for a in rx)
    syy = sum((b - my) ** 2 for b in ry)
    if sxx == 0 or syy == 0:
        return None
    return sxy / (sxx * syy) ** 0.5


def fmt(value):
    if value is None:
        return "NA"
    return str(int(value)) if float(value).is_integer() else f"{value:.2f}"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--phylo", required=True, help="the run's 06_phylogeny folder")
    ap.add_argument("--threshold", type=float, default=20, help="SNPs at or below which a pair counts as close [%(default)s]")
    ap.add_argument("--out", default=None, help="output prefix [<phylo>/snippy_vs_ska2]")
    args = ap.parse_args()

    out = args.out or os.path.join(args.phylo, "snippy_vs_ska2")
    groups = sorted(d for d in os.listdir(args.phylo)
                    if os.path.isdir(os.path.join(args.phylo, d)) and not d.endswith("_ska2")
                    and os.path.isdir(os.path.join(args.phylo, d + "_ska2")))
    if not groups:
        sys.exit(f"No group in {args.phylo} has both <group>/ and <group>_ska2/: was the pipeline run with --phylo_method both?")

    pair_rows, summary_rows = [], []
    t = args.threshold
    for g in groups:
        snippy_file = os.path.join(args.phylo, g, f"{g}.snp_distances.tsv")
        ska_file = os.path.join(args.phylo, f"{g}_ska2", f"{g}_ska2.snp_distances.tsv")
        if not (os.path.exists(snippy_file) and os.path.exists(ska_file)):
            print(f"{g}: missing {snippy_file if not os.path.exists(snippy_file) else ska_file}; skipped", file=sys.stderr)
            continue
        snippy, ska = read_matrix(snippy_file), read_matrix(ska_file)
        ref_free = read_ska_distance(os.path.join(args.phylo, f"{g}_ska2", f"{g}_ska2.ska_distances.tsv"))
        pairs = sorted(set(snippy) & set(ska))
        if not pairs:
            print(f"{g}: no isolate pairs in common; skipped", file=sys.stderr)
            continue
        for a, b in pairs:
            s, k = snippy[(a, b)], ska[(a, b)]
            pair_rows.append({
                "group": g, "isolate1": a, "isolate2": b,
                "snippy": fmt(s), "ska2": fmt(k), "ska2_ref_free": fmt(ref_free.get((a, b))),
                "difference": fmt(k - s),
                "close_snippy": "yes" if s <= t else "no", "close_ska2": "yes" if k <= t else "no",
            })
        xs = [snippy[p] for p in pairs]
        ys = [ska[p] for p in pairs]
        diffs = [abs(y - x) for x, y in zip(xs, ys)]
        close_s = {p for p in pairs if snippy[p] <= t}
        close_k = {p for p in pairs if ska[p] <= t}
        summary_rows.append({
            "group": g, "pairs": len(pairs),
            "spearman_rho": fmt(spearman(xs, ys)),
            "median_abs_difference": fmt(statistics.median(diffs)),
            "max_abs_difference": fmt(max(diffs)),
            f"close_snippy_le{fmt(t)}": len(close_s), f"close_ska2_le{fmt(t)}": len(close_k),
            "close_in_both": len(close_s & close_k), "close_in_one_only": len(close_s ^ close_k),
        })

    if not summary_rows:
        sys.exit("Nothing to compare.")
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    for suffix, rows in (("_pairs.tsv", pair_rows), ("_summary.tsv", summary_rows)):
        with open(out + suffix, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(rows)

    print(f"Snippy vs SKA2 SNP distances (close = {fmt(t)} SNPs or fewer)")
    cols = list(summary_rows[0].keys())
    widths = [max(len(c), *(len(str(r[c])) for r in summary_rows)) for c in cols]
    print("  ".join(c.ljust(w) for c, w in zip(cols, widths)))
    for r in summary_rows:
        print("  ".join(str(r[c]).ljust(w) for c, w in zip(cols, widths)))
    print(f"\nPairs: {out}_pairs.tsv\nSummary: {out}_summary.tsv")


if __name__ == "__main__":
    main()
