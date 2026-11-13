#!/usr/bin/env python3
"""
54_supplement_tables.py

Supplemental Tables S1-S6 of the MAP manuscript, one Excel file each (a
README sheet holds the legend), written to 11_supplement/.

  S1  cohort: accession, study, platform, metadata, QC, lineage, typing,
      clade D, read depth, mixed-site count/flag; second sheet: excluded
      genomes with reasons
  S2  assembly and annotation statistics (QUAST, CheckM2, Prokka counts
      from each genome's GFF)
  S3  pangenome summary, Panaroo vs Roary (summary_statistics.txt)
  S4  lineage- or clade-specific genes confirmed by reads (49 PART B)
  S5  gene-level convergence test for all K-10 CDS (49)
  S6  origin prediction: balanced accuracy per label/scheme/model and per
      held-out study (32)

  python 54_supplement_tables.py
Environment: ml_env (+ openpyxl; pip install openpyxl if missing; without it
TSV files are written instead).
"""

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/cfs/earth/scratch/tabasran/MAP")

LEGENDS = {
    "S1": "TABLE S1. Genomes analysed. One row per genome of the final cohort (n = 1,792): accession, "
          "study and sequencing platform (ENA), curated metadata, assembly quality, lineage assignment, "
          "typing level, clade D membership, read depth on K-10 and number of mixed core SNP sites "
          "(possible mixed-strain sample if >= 20). Sheet 'excluded' lists genomes removed after "
          "annotation, with the reason.",
    "S2": "TABLE S2. Assembly and annotation statistics: QUAST metrics, CheckM2 completeness and "
          "contamination, and numbers of CDS, tRNA and rRNA features annotated by Prokka.",
    "S3": "TABLE S3. Pangenome summary from Panaroo and Roary run on the same 1,792 annotations.",
    "S4": "TABLE S4. Lineage-specific genes. Genes present in at least 95% of the genomes with a three-state "
          "read call on one side of a lineage contrast and in at most 5% on the other. For each gene: the "
          "contrast, presence frequency and number of called genomes per side, position on K-10 or on the "
          "S-type reference Telford, the block of genes less than 5 kb apart that it belongs to, the "
          "overlapping K-10 coding sequence, and the eggNOG-mapper annotation. Sheet 'blocks': one row per "
          "block.",
    "S5": "TABLE S5. Gene-level convergence test. For every K-10 coding sequence outside masked repeats: "
          "number of SNP sites, number of distinct tree branches with a change (read-verified calls), "
          "expected number from gene length, Poisson P value and Benjamini-Hochberg q value.",
    "S6": "TABLE S6. Prediction of geographic origin from core SNPs: balanced accuracy by label, validation "
          "scheme and model (sheet 'metrics') and accuracy for each held-out study (sheet 'per_study').",
}


def write(out, name, sheets):
    try:
        import openpyxl  # noqa: F401
        p = out / "Table{}.xlsx".format(name)
        with pd.ExcelWriter(p) as xw:
            pd.DataFrame({"legend": [LEGENDS[name]]}).to_excel(xw, sheet_name="README", index=False)
            for sn, df in sheets.items():
                df.to_excel(xw, sheet_name=sn[:31], index=False)
    except ImportError:
        for sn, df in sheets.items():
            df.to_csv(out / "Table{}_{}.tsv".format(name, sn), sep="\t", index=False)
        (out / "Table{}_README.txt".format(name)).write_text(LEGENDS[name] + "\n")
        p = out
    print("  Table {}: {}".format(name, ", ".join("{} {}x{}".format(k, *v.shape) for k, v in sheets.items())))
    return p


def read(p, **kw):
    p = Path(p)
    if not p.exists():
        print("  [i] missing", p)
        return None
    return pd.read_csv(p, **kw)


def keycol(df, ids):
    """column of df whose values overlap the cohort ids most."""
    best, n = None, 0
    for c in df.columns:
        k = df[c].astype(str).str.strip().isin(ids).sum()
        if k > n:
            best, n = c, k
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=str(ROOT))
    args = ap.parse_args()
    root = Path(args.root)
    out = root / "11_supplement"
    out.mkdir(parents=True, exist_ok=True)
    ids = [l.strip() for l in open(root / "01_metadata/cohort_final.txt") if l.strip()]
    T = pd.DataFrame({"genome": ids})

    # ---------------------------------------------------------------- S1
    ena = read(root / "01_metadata/ena_raw.tsv", sep="\t", dtype=str)
    if ena is not None:
        ena = ena.drop_duplicates("run_accession")
        cols = [c for c in ("run_accession", "study_accession", "sample_accession", "instrument_platform",
                            "instrument_model", "library_layout", "read_count", "first_public") if c in ena.columns]
        T = T.merge(ena[cols].rename(columns={"run_accession": "genome"}), on="genome", how="left")
    md = read(root / "01_metadata/metadata_clean.tsv", sep="\t", dtype=str)
    if md is not None:
        k = keycol(md, set(ids))
        md = md.drop_duplicates(k).rename(columns={k: "genome"})
        md.columns = ["genome"] + ["meta_" + c for c in md.columns if c != "genome"]
        T = T.merge(md, on="genome", how="left")
    dates = read(root / "09_dating/dates.tsv", sep="\t", dtype=str)
    if dates is not None:
        T = T.merge(dates[["genome", "raw", "precision"]].rename(columns={"raw": "collection_date",
                                                                         "precision": "date_precision"}),
                    on="genome", how="left")
    lin = read(root / "05_phylogeny/typing/lineage_assignments.tsv", sep="\t", dtype=str)
    if lin is not None:
        T = T.merge(lin, on="genome", how="left")
    sp = read(root / "05_phylogeny/typing/split_assignments.tsv", sep="\t", dtype=str)
    if sp is not None:
        T = T.merge(sp[["genome", "tree_side", "two_snp_call"]].rename(
            columns={"tree_side": "C_split_group", "two_snp_call": "clade_D_two_snp_call"}), on="genome", how="left")
    q = read(root / "02_qc/genome_qc_master.tsv", sep="\t")
    rq = read(root / "08_reads/summary/genome_qc_reads.tsv", sep="\t")
    if rq is not None:
        rq = rq.rename(columns={rq.columns[0]: "genome"})
        keep = [c for c in ("genome", "depth", "breadth", "mixed_sites") if c in rq.columns]
        T = T.merge(rq[keep].rename(columns={"depth": "K10_read_depth", "breadth": "K10_breadth_pct"}),
                    on="genome", how="left")
        T["possible_mixed_sample"] = np.where(T["mixed_sites"].isna(), "no reads",
                                              np.where(T["mixed_sites"] >= 20, "yes", "no"))
    exc = []
    e3 = read(root / "02_qc/excluded_step3_assembly_qc.tsv", sep="\t", dtype=str)
    if e3 is not None:
        e3 = e3.copy()
        e3["step"] = "CheckM2 (completeness < 95% or contamination > 5%)"
        exc.append(e3)
    e4 = read(root / "01_metadata/manual_exclusions.tsv", sep="\t", header=None, names=["genome", "reason"], dtype=str)
    if e4 is not None:
        e4["step"] = "manual check"
        exc.append(e4)
    sheets = {"genomes": T}
    if exc:
        sheets["excluded"] = pd.concat(exc, ignore_index=True, sort=False)
    write(out, "S1", sheets)

    # ---------------------------------------------------------------- S2
    S2 = pd.DataFrame({"genome": ids})
    if q is not None:
        kq = "genome" if "genome" in q.columns else keycol(q.astype(str), set(ids))
        S2 = S2.merge(q.rename(columns={kq: "genome"}), on="genome", how="left")
    counts = []
    for g in ids:
        f = root / "03_genomes/gff" / (g + ".gff")
        c = {"genome": g, "CDS": 0, "tRNA": 0, "rRNA": 0, "tmRNA": 0}
        if f.exists():
            for line in open(f):
                if line.startswith("##FASTA"):
                    break
                if line.startswith("#"):
                    continue
                x = line.split("\t", 3)
                if len(x) > 2 and x[2] in c:
                    c[x[2]] += 1
        else:
            c = {"genome": g}
        counts.append(c)
    S2 = S2.merge(pd.DataFrame(counts).rename(columns={k: "prokka_" + k for k in ("CDS", "tRNA", "rRNA", "tmRNA")}),
                  on="genome", how="left")
    write(out, "S2", {"assembly_annotation": S2})

    # ---------------------------------------------------------------- S3
    rows = []
    for tool in ("panaroo", "roary"):
        p = root / "04_pangenome" / tool / "summary_statistics.txt"
        if not p.exists():
            print("  [i] missing", p)
            continue
        for line in open(p):
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3:
                rows.append({"tool": tool, "category": f[0], "definition": f[1], "genes": f[2]})
    write(out, "S3", {"summary": pd.DataFrame(rows)})

    # ---------------------------------------------------------------- S4
    sh4 = {}
    L = read(root / "07_functional/lineage_genes/lineage_specific_genes.tsv", sep="\t")
    if L is not None:
        sh4["genes"] = L
    BK = read(root / "07_functional/lineage_genes/blocks.tsv", sep="\t")
    if BK is not None:
        sh4["blocks"] = BK
    if sh4:
        write(out, "S4", sh4)

    # ---------------------------------------------------------------- S5
    G = read(root / "06_pangwas/convergence_v2/gene_branch_convergence.tsv", sep="\t")
    if G is not None:
        write(out, "S5", {"gene_convergence": G})

    # ---------------------------------------------------------------- S6
    M = read(root / "07_ml/C_nonB/metrics.tsv", sep="\t")
    P = read(root / "07_ml/C_nonB/per_study.tsv", sep="\t")
    sh = {k: v for k, v in (("metrics", M), ("per_study", P)) if v is not None}
    if sh:
        write(out, "S6", sh)
    print("written:", out)


if __name__ == "__main__":
    main()
