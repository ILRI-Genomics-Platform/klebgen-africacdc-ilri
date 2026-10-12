#!/usr/bin/env python3
"""
select_training_set.py - choose a training set from the supplementary tables
of Heinz et al., Genome Medicine 2024 (doi:10.1186/s13073-024-01342-3;
Blantyre, Malawi, 1998-2020).

Usage
  python3 select_training_set.py MOESM2_ESM.xlsx OUT_DIR [--size full|compact]

  --size full     65 isolates: three candidate transmission clusters
  --size compact  29 isolates: two clusters, fewer diversity examples
                  (a strict subset of 'full', so pre-computed full results
                  cover it)
  --size capstone 36 isolates: 'full' minus 'compact' (unseen data for the
                  Day 10 capstone, including the ST14 PICU/HDU 2019 cluster)

Outputs (in OUT_DIR, suffixed with the size)
  training_set_<size>.tsv         one row per isolate, with group + reason
  accessions_training_<size>.tsv  ready for bin/fetch_ena.py fetch
  metadata_<size>.tsv             sample, collection year, ward, neonate,
                                  published ST / K / O (for Modules 6 and 8)
  ast_<size>.tsv                  disc-diffusion AST (S/I/R) from Table S4,
                                  drugs tested in the set only (Module 4)
  published_K_O_all.tsv           K/O types, year, ST and neonate flag for all
                                  1,373 study isolates (Module 7 coverage)

Rules (deterministic: same input -> same selection)
  * Good quality = species_match 'strong', <=100 contigs, total size
    5.1-6.0 Mb, no QC warning other than ambiguous_bases.
  * Prefer isolates with AST results (Table S4) for the geno-pheno module.
  * Transmission clusters: within a neonatal ST/ward/year group, take the
    most common genotype profile (K locus + YbST + acquired AMR genes +
    plasmid replicons), optionally plus isolates with a minor profile.
    Background isolates = same ST from the most distant and nearest years.
  * Long-read isolates come from Table S6 (Illumina run + ONT experiment).
  * Requires pandas and openpyxl.
"""

import argparse
import pandas as pd

# How many isolates each group gets, per size. Compact counts are always <=
# full counts and use the same ordering, so compact is a subset of full.
SIZES = {
    "full": {
        "clusters": [  # ST, year, ward, n_dominant, n_minor, n_background
            ("ST39", 2017, "CHATINKHA NURSERY", 14, 0, 3),
            ("ST35", 2018, "CHATINKHA NURSERY", 6, 2, 3),
            ("ST14", 2019, "PICU/HDU", 8, 0, 3),
        ],
        "diversity": ["ST15", "ST307", "ST17", "ST101", "ST340", "ST25", "ST1552", "ST86", "ST48", "ST607"],
        "virulence": 3, "per_species": 2, "qc_fragmented": 2, "qc_size": 2,
    },
    "compact": {
        "clusters": [
            ("ST39", 2017, "CHATINKHA NURSERY", 8, 0, 2),
            ("ST35", 2018, "CHATINKHA NURSERY", 4, 1, 1),
        ],
        "diversity": ["ST15", "ST307", "ST17"],
        "virulence": 1, "per_species": 1, "qc_fragmented": 1, "qc_size": 1,
    },
}


def load(xlsx):
    s1 = pd.read_excel(xlsx, "TableS1")
    s3 = pd.read_excel(xlsx, "TableS3")
    s4 = pd.read_excel(xlsx, "TableS4")
    s5 = pd.read_excel(xlsx, "TableS5")
    s6 = pd.read_excel(xlsx, "TableS6")

    # Table S3 has no ID column; its rows follow the order of Table S2's 'name'
    # (checked: replicons match the closed plasmids of the Table S6 isolates).
    s3.index = pd.read_excel(xlsx, "TableS2")["name"].values
    rep_cols = list(s3.columns)
    s1["replicons"] = s1["strain"].map(
        s3.apply(lambda r: ";".join(c for c in rep_cols if r[c] == 1), axis=1))

    ast = s4.set_index("strain").replace("NA", pd.NA)
    s1["n_ast"] = s1["strain"].map(ast.notna().sum(axis=1)).fillna(0).astype(int)
    s1 = s1.merge(s5[["Assembly", "Otype", "Ktype"]], left_on="strain", right_on="Assembly", how="left")

    amr_cols = [c for c in s1.columns if c.endswith("_acquired")]
    s1["profile"] = (s1["K_locus"].astype(str) + "/" + s1["YbST"].astype(str) + "/"
                     + s1[amr_cols].astype(str).agg("|".join, axis=1) + "/" + s1["replicons"].fillna(""))
    s1["good_qc"] = ((s1["species_match"] == "strong") & (s1["contig_count"] <= 100)
                     & s1["total_size"].between(5.1e6, 6.0e6)
                     & s1["QC_warnings"].isin(["ambiguous_bases", "-"]))
    return s1, s6, ast


def select(s1, s6, cfg):
    picked = {}  # strain -> dict(group, reason[, ont])

    def take(df, n, group, reason, prefer_ast=True):
        """Pick n not-yet-chosen isolates. prefer_ast: most AST results, then
        fewest contigs; otherwise keep the caller's order."""
        if n <= 0:
            return
        df = df[~df["strain"].isin(picked)]
        if prefer_ast:
            df = df.sort_values(["n_ast", "contig_count", "strain"], ascending=[False, True, True])
        for _, r in df.head(n).iterrows():
            picked[r["strain"]] = {"group": group, "reason": reason}

    good = s1[s1["good_qc"]]

    # 1. Long-read (ONT) isolates from Table S6 -> Modules 2, 5
    for sample, srx in s6.groupby("sample")["long-read accession"].first().items():
        row = s1[s1["SupplierName"].astype(str) == str(sample)].iloc[0]
        picked[row["strain"]] = {"group": "A_hybrid", "ont": srx,
                                 "reason": f"Illumina + ONT ({srx}); closed plasmids in Table S6"}

    # 2. Candidate transmission clusters + background -> Module 6
    for st, year, ward, n_dom, n_minor, n_bg in cfg["clusters"]:
        tag = f"B_cluster_{st}_{year}"
        grp = good[(good["main_ST"] == st) & (good["DATE"] == year) & (good["major_wards"] == ward)]
        top = grp["profile"].value_counts()
        dom = grp[grp["profile"] == top.index[0]]
        take(dom, n_dom, tag, f"neonatal {ward.title()} {year}; dominant {st} genotype "
             f"({len(dom)}/{len(grp)} share it) - candidate transmission")
        take(grp[grp["profile"] != top.index[0]], n_minor, tag,
             f"{ward.title()} {year}, {st}, different genotype - test whether linked")
        others = good[(good["main_ST"] == st) & (good["DATE"] != year)]
        years = sorted(others["DATE"].dropna().unique(), key=lambda y: (abs(y - year), y))
        for yr in ([years[-1]] + years[:2])[:n_bg]:   # most distant first, then nearest
            take(others[others["DATE"] == yr], 1, tag + "_background",
                 f"{st} from {int(yr)} - background to test against the cluster")

    # 3. Carbapenemase isolates (NDM-1) -> Module 4, kept regardless of contig count
    for _, r in s1[s1["Bla_Carb_acquired"] != "-"].iterrows():
        qc = "" if r["good_qc"] else f" ({int(r['contig_count'])} contigs - borderline QC)"
        picked[r["strain"]] = {"group": "C_carbapenemase", "reason": f"{r['Bla_Carb_acquired']} carbapenemase{qc}"}

    # 4. Lineage / AMR / K-O diversity -> Modules 3, 4, 7
    for st in cfg["diversity"]:
        take(good[good["main_ST"] == st], 1, "D_diversity", f"common / high-risk lineage {st}")

    # 5. Virulence markers -> Module 4
    take(good[good["virulence_score"] >= 3], cfg["virulence"], "E_virulence",
         "virulence score >=3 (aerobactin and/or rmp) - hypervirulence markers")

    # 6. Species complex -> Module 3
    for sp in ["Klebsiella quasipneumoniae subsp. similipneumoniae", "Klebsiella variicola subsp. variicola"]:
        take(good[good["species"] == sp], cfg["per_species"], "F_species",
             f"{sp} - often misidentified as K. pneumoniae")

    # 7. Natural QC failures -> Module 1
    bad = s1[~s1["good_qc"]]
    take(bad.sort_values("contig_count", ascending=False), cfg["qc_fragmented"], "G_qc_fail",
         "very fragmented assembly - QC failure exercise", prefer_ast=False)
    take(bad[bad["QC_warnings"].str.contains("total_size", na=False)].sort_values("total_size", ascending=False),
         cfg["qc_size"], "G_qc_fail",
         "assembly ~9 Mb, far above expected (likely mixed culture) - QC failure exercise", prefer_ast=False)

    return picked


def write(s1, ast, picked, out_dir, size):
    sel = s1[s1["strain"].isin(picked)].copy()
    sel["group"] = sel["strain"].map(lambda s: picked[s]["group"])
    sel["reason"] = sel["strain"].map(lambda s: picked[s]["reason"])
    sel["ont_accession"] = sel["strain"].map(lambda s: picked[s].get("ont", ""))
    sel = sel.sort_values(["group", "DATE", "strain"])

    cols = ["group", "strain", "SupplierName", "LaneAcc", "SampleAcc", "ont_accession", "species", "ST",
            "K_locus", "Otype", "DATE", "neonate", "major_wards", "contig_count", "total_size",
            "virulence_score", "resistance_score", "Bla_ESBL_acquired", "Bla_Carb_acquired",
            "replicons", "n_ast", "reason"]
    sel[cols].to_csv(f"{out_dir}/training_set_{size}.tsv", sep="\t", index=False)

    # phylo_group: isolates of a cluster ST (clusters, backgrounds and any other
    # good-quality isolate of that ST, e.g. a hybrid) share one phylogeny
    cluster_sts = {c[0] for c in SIZES.get(size, SIZES["full"])["clusters"]}
    phylo_group = [st if (st in cluster_sts and g != "G_qc_fail") else ""
                   for st, g in zip(sel["main_ST"], sel["group"])]
    pd.DataFrame({
        "id": sel["SupplierName"].astype(str).str.replace(r"[^A-Za-z0-9_.-]", "_", regex=True),
        "species": "kpneumoniae",   # all KpSC use the workflow's kpsc preset
        "accession": [f"{l};{o}" if o else l for l, o in zip(sel["LaneAcc"], sel["ont_accession"])],
        "phylo_group": phylo_group,
    }).to_csv(f"{out_dir}/accessions_training_{size}.tsv", sep="\t", index=False)

    sample_ids = sel["SupplierName"].astype(str).str.replace(r"[^A-Za-z0-9_.-]", "_", regex=True)

    # Metadata for trainees (Modules 6 and 8). Published typing is kept so
    # trainees can compare their own results with the paper's.
    meta = pd.DataFrame({
        "sample": sample_ids.values,
        "collection_year": sel["DATE"].astype("Int64").values,
        "ward": sel["major_wards"].values,
        "neonate": sel["neonate"].map({"y": "yes", "n": "no", "u": "unknown"}).values,
        "training_group": sel["group"].values,
        "published_species": sel["species"].values,
        "published_ST": sel["ST"].values,
        "published_K_locus": sel["K_locus"].values,
        "published_O_type": sel["Otype"].values,
        "ena_run": sel["LaneAcc"].values,
        "ena_ont_experiment": sel["ont_accession"].values,
    })
    meta.to_csv(f"{out_dir}/metadata_{size}.tsv", sep="\t", index=False)

    # AST results (Table S4): keep drugs tested in at least one selected isolate
    a = ast.reindex(sel["strain"])
    a = a.loc[:, a.notna().any()]
    a.insert(0, "sample", sample_ids.values)
    a.to_csv(f"{out_dir}/ast_{size}.tsv", sep="\t", index=False, na_rep="NA")

    print(f"[{size}] {len(sel)} isolates; {(sel.n_ast > 0).sum()} with AST; "
          f"{(sel.neonate == 'y').sum()} neonatal")
    print(sel.groupby("group").size().to_string())
    return sel


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("xlsx", help="13073_2024_1342_MOESM2_ESM.xlsx")
    ap.add_argument("out_dir", nargs="?", default=".")
    ap.add_argument("--size", choices=sorted(SIZES) + ["capstone"], default="full")
    args = ap.parse_args()
    s1, s6, ast = load(args.xlsx)
    if args.size == "capstone":
        full = select(s1, s6, SIZES["full"])
        compact = select(s1, s6, SIZES["compact"])
        picked = {k: v for k, v in full.items() if k not in compact}
    else:
        picked = select(s1, s6, SIZES[args.size])
    write(s1, ast, picked, args.out_dir, args.size)
    write_ko_table(args.xlsx, s1, args.out_dir)


def write_ko_table(xlsx, s1, out_dir):
    """All study isolates' published K/O types, for vaccine-coverage exercises."""
    s5 = pd.read_excel(xlsx, "TableS5")
    ko = s5.merge(s1[["strain", "SupplierName", "neonate", "ST"]], left_on="Assembly", right_on="strain", how="left")
    ko = pd.DataFrame({
        "sample": ko["SupplierName"].astype(str),
        "collection_year": ko["DATE"].astype("Int64"),
        "neonate": ko["neonate"].map({"y": "yes", "n": "no", "u": "unknown"}),
        "ST": ko["ST"],
        "K_locus": ko["Klocus"],
        "K_type": ko["Ktype"],
        "O_locus": ko["Olocus"],
        "O_type": ko["Otype"],
    })
    ko.to_csv(f"{out_dir}/published_K_O_all.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
