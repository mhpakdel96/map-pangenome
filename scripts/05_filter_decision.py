#!/usr/bin/env python3
"""
05_filter_decision.py

Merge the QUAST and CheckM2 results into one per-genome QC table, then show
how many genomes each candidate threshold combination would remove.

This script DECIDES NOTHING and DELETES NOTHING. It produces the evidence
needed to pick thresholds; 05_apply_filter.py then applies whatever is chosen.

Usage
-----
    python 05_filter_decision.py                       # explore thresholds
    python 05_filter_decision.py --contigs 800 --n50 10000 \
                                 --completeness 90 --contamination 5

Outputs
-------
    02_qc/genome_qc_master.tsv    one row per genome, all QC metrics
    02_qc/threshold_scan.tsv      genomes removed under each threshold combo
"""


import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/cfs/earth/scratch/tabasran/MAP")
QUAST = ROOT / "02_qc/quast/all_genomes/transposed_report.tsv"
CHECKM = ROOT / "02_qc/checkm2/all_genomes/quality_report.tsv"
META = ROOT / "01_metadata/Supplementary_Table_S1_clean.tsv"   # optional
OUT_MASTER = ROOT / "02_qc/genome_qc_master.tsv"
OUT_SCAN = ROOT / "02_qc/threshold_scan.tsv"


def load_quast():
    if not QUAST.exists():
        sys.exit(f"[!] QUAST output not found: {QUAST}\n    Run 03_quast.sbatch first.")
    df = pd.read_csv(QUAST, sep="\t")
    df = df.rename(columns={df.columns[0]: "genome"})
    df["genome"] = df["genome"].astype(str).str.replace(r"\.fasta$", "", regex=True)

    keep = {
        "# contigs": "contigs",
        "Total length": "length",
        "GC (%)": "gc",
        "N50": "n50",
        "L50": "l50",
        "Largest contig": "largest_contig",
    }
    missing = [c for c in keep if c not in df.columns]
    if missing:
        sys.exit(f"[!] QUAST report missing expected columns: {missing}")

    out = df[["genome"] + list(keep)].rename(columns=keep)
    for c in keep.values():
        out[c] = pd.to_numeric(out[c], errors="coerce")

    dup = out["genome"].duplicated().sum()
    if dup:
        print(f"[!] {dup} duplicate genome IDs in QUAST output — keeping first of each.")
        out = out[~out["genome"].duplicated(keep="first")]
    return out


def load_checkm():
    if not CHECKM.exists():
        print(f"[i] CheckM2 output not found ({CHECKM}). Continuing without it.")
        return None
    df = pd.read_csv(CHECKM, sep="\t")
    ren = {"Name": "genome", "Completeness": "completeness",
           "Contamination": "contamination"}
    df = df.rename(columns=ren)
    if "genome" not in df.columns:
        sys.exit("[!] CheckM2 report has no 'Name' column — unexpected format.")
    df["genome"] = df["genome"].astype(str).str.replace(r"\.fasta$", "", regex=True)
    cols = ["genome", "completeness", "contamination"]
    for extra in ("Contig_N50", "Genome_Size", "GC_Content", "Total_Coding_Sequences"):
        if extra in df.columns:
            cols.append(extra)
    return df[cols]


def describe(df, cols):
    pct = [1, 5, 25, 50, 75, 95, 99]
    print(f"\n{'metric':>16} " + " ".join(f"{'p'+str(p):>10}" for p in pct))
    for c in cols:
        if c not in df.columns:
            continue
        v = df[c].dropna()
        if v.empty:
            continue
        q = np.percentile(v, pct)
        print(f"{c:>16} " + " ".join(f"{x:>10.0f}" for x in q))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--contigs", type=int, default=None,
                    help="max contigs allowed")
    ap.add_argument("--n50", type=int, default=None,
                    help="min N50 allowed")
    ap.add_argument("--completeness", type=float, default=None,
                    help="min CheckM2 completeness %%")
    ap.add_argument("--contamination", type=float, default=None,
                    help="max CheckM2 contamination %%")
    args = ap.parse_args()

    q = load_quast()
    c = load_checkm()
    df = q if c is None else q.merge(c, on="genome", how="outer")

    df["source"] = np.where(df["genome"].str.startswith("GCA"), "reference", "de_novo")

    if META.exists():
        m = pd.read_csv(META, sep="\t")
        idcol = next((x for x in ("accession", "genome", "Isolate") if x in m.columns), None)
        if idcol:
            m = m.rename(columns={idcol: "genome"})
            df = df.merge(m, on="genome", how="left")

    OUT_MASTER.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_MASTER, sep="\t", index=False)

    print("=" * 63)
    print(f" Genomes in QC table : {len(df)}")
    print(f"   de novo           : {(df['source'] == 'de_novo').sum()}")
    print(f"   reference (GCA)   : {(df['source'] == 'reference').sum()}")
    print(f" Written to          : {OUT_MASTER}")
    print("=" * 63)

    describe(df, ["contigs", "n50", "length", "gc", "completeness", "contamination"])

    # Reference genomes reported separately: they were never QC'd before and
    # they carry the highest CDS counts in the cohort.
    if (df["source"] == "reference").any():
        print("\n--- reference (GCA) genomes only ---")
        describe(df[df["source"] == "reference"],
                 ["contigs", "n50", "completeness", "contamination"])

    # ---------------------------------------------------------------- scan
    grid_contigs = [None, 1200, 1000, 800, 600, 500, 400, 300]
    grid_n50 = [None, 5000, 7500, 10000, 15000, 20000, 25000]
    rows = []
    for gc_ in grid_contigs:
        for gn in grid_n50:
            if gc_ is None and gn is None:
                continue
            bad = pd.Series(False, index=df.index)
            if gc_ is not None:
                bad |= df["contigs"] > gc_
            if gn is not None:
                bad |= df["n50"] < gn
            rows.append({"max_contigs": gc_, "min_n50": gn,
                         "removed": int(bad.sum()), "retained": int((~bad).sum())})
    scan = pd.DataFrame(rows)
    scan.to_csv(OUT_SCAN, sep="\t", index=False)

    print("\n=== assembly threshold scan (contiguity only) ===")
    for gc_ in [1000, 800, 600, 500, 400]:
        for gn in [5000, 10000, 15000, 20000]:
            r = scan[(scan.max_contigs == gc_) & (scan.min_n50 == gn)]
            if not r.empty:
                print(f"  contigs<={gc_:<5} n50>={gn:<6} -> removed {r.removed.iloc[0]:>4}"
                      f"  retained {r.retained.iloc[0]:>5}")

    if c is not None:
        print("\n=== CheckM2 threshold scan ===")
        for comp in [85, 90, 93, 95]:
            for cont in [10, 5, 3, 2]:
                bad = (df["completeness"] < comp) | (df["contamination"] > cont)
                print(f"  completeness>={comp:<3} contamination<={cont:<3}"
                      f" -> removed {int(bad.sum()):>4}  retained {int((~bad).sum()):>5}")

    # ---------------------------------------------------- chosen thresholds
    chosen = [a is not None for a in (args.contigs, args.n50,
                                      args.completeness, args.contamination)]
    if any(chosen):
        fail = pd.Series(False, index=df.index)
        reasons = {}
        if args.contigs is not None:
            reasons["contigs"] = df["contigs"] > args.contigs
        if args.n50 is not None:
            reasons["n50"] = df["n50"] < args.n50
        if args.completeness is not None and "completeness" in df:
            reasons["completeness"] = df["completeness"] < args.completeness
        if args.contamination is not None and "contamination" in df:
            reasons["contamination"] = df["contamination"] > args.contamination

        for r in reasons.values():
            fail |= r.fillna(False)

        df["fail_reason"] = [
            ",".join(k for k, v in reasons.items() if bool(v.iloc[i]))
            for i in range(len(df))
        ]
        df["passes_filter"] = ~fail
        df.to_csv(OUT_MASTER, sep="\t", index=False)

        print("\n" + "=" * 63)
        print(" APPLIED THRESHOLDS")
        for k, v in vars(args).items():
            if v is not None:
                print(f"   {k}: {v}")
        print(f" Removed  : {int(fail.sum())}")
        print(f" Retained : {int((~fail).sum())}")
        print("\n Breakdown by reason (genomes may fail more than one):")
        for k, v in reasons.items():
            print(f"   {k:>14}: {int(v.fillna(False).sum())}")
        print("\n Retained by source:")
        print(df[~fail].groupby("source").size().to_string())
        print("=" * 63)

        removed = df[fail][["genome", "source", "contigs", "n50", "fail_reason"]]
        path = ROOT / "02_qc/genomes_failing_filter.tsv"
        removed.to_csv(path, sep="\t", index=False)
        print(f"\n Failing genomes written to: {path}")
        print(removed.head(25).to_string(index=False))


if __name__ == "__main__":
    main()
