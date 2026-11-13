#!/usr/bin/env python3
"""
36_read_ternary.py

Three-state gene calls from read breadth, and the same negative control.

Why. 35_read_presence.py showed that a single breadth cut-off (>= 80%) makes
read-based gene counts depend on sequencing depth (Spearman +0.74): at low
depth part of a gene that IS present goes uncovered and falls below 80%.
A genuinely absent gene has (near) zero breadth; a present gene at low depth
has intermediate breadth. So:
    present  breadth >= --present  (default 80)
    absent   breadth <  --absent   (default 20)
    NA       in between — not called, left out of every test for that gene
NA cells are not imputed. Reported: share of NA, category counts on called
cells, genes per genome vs N50 / depth, and the within-study negative
control (low vs high N50, low vs high depth; Fisher on called cells,
Bonferroni per study) next to the numbers of union80 and binary reads from 35.

Outputs (08_reads/summary/)
  presence_reads_ternary.tsv   1 / 0 / NA, kept clusters >= 250 bp x genomes
  ternary_report.txt
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, spearmanr

ROOT = Path("/cfs/earth/scratch/tabasran/MAP")


def fisher_count(M, grp, alpha=0.05):
    """M: genes x genomes float (1/0/NaN). NaN cells excluded per gene."""
    known = ~np.isnan(M)
    pres = (M == 1)
    a = (pres & grp).sum(axis=1)
    b = (known & ~pres & grp).sum(axis=1)
    c = (pres & ~grp).sum(axis=1)
    d = (known & ~pres & ~grp).sum(axis=1)
    var = (a + c > 0) & (b + d > 0) & (a + b > 0) & (c + d > 0)
    idx = np.where(var)[0]
    n = len(idx)
    sig = sum(1 for i in idx if fisher_exact([[a[i], b[i]], [c[i], d[i]]])[1] < alpha / max(n, 1))
    return sig, n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=str(ROOT / "08_reads/summary"))
    ap.add_argument("--ena", default=str(ROOT / "01_metadata/ena_raw.tsv"))
    ap.add_argument("--present", type=float, default=80.0)
    ap.add_argument("--absent", type=float, default=20.0)
    ap.add_argument("--min-study", type=int, default=60)
    args = ap.parse_args()
    d = Path(args.dir)
    R = []

    def say(s=""):
        print(s, flush=True)
        R.append(s)

    B = pd.read_csv(d / "breadth.tsv.gz", sep="\t", index_col=0)
    Q = pd.read_csv(d / "genome_qc_reads.tsv", sep="\t", index_col=0).reindex(B.columns)
    G = list(B.columns)
    b = B.values
    M = np.where(b >= args.present, 1.0, np.where(b < args.absent, 0.0, np.nan))
    pd.DataFrame(M, index=B.index, columns=G).to_csv(d / "presence_reads_ternary.tsv", sep="\t",
                                                     na_rep="NA", float_format="%.0f")
    say("ternary calls: present >= {}%, absent < {}%, NA in between".format(args.present, args.absent))
    na = np.isnan(M)
    say("cells: present {:.2f}%  absent {:.2f}%  NA {:.2f}%".format(
        100 * (M == 1).mean(), 100 * (M == 0).mean(), 100 * na.mean()))
    hist = pd.cut(pd.Series(b.ravel()), [-0.1, 0, 20, 50, 80, 99.9, 100]).value_counts().sort_index()
    say("breadth distribution of all cells:\n" + (100 * hist / hist.sum()).round(3).to_string())
    na_g = na.mean(axis=0)
    say("NA share per genome: median {:.4f}, p95 {:.4f}; Spearman with depth {:+.3f}".format(
        np.median(na_g), np.quantile(na_g, .95), spearmanr(na_g, Q["depth"])[0]))

    # genes per genome on clusters with no NA in that genome is not comparable across
    # genomes; use the share of called clusters that are present
    called = (~na).sum(axis=0)
    share = (M == 1).sum(axis=0) / called
    f = np.nanmean(M, axis=1)
    acc = (f < 0.99) & (f > 0)
    acc_share = np.nansum(M[acc] == 1, axis=0) / (~na[acc]).sum(axis=0)
    say("share of called clusters present, per genome — Spearman: all genes N50 {:+.3f} depth {:+.3f} | "
        "accessory N50 {:+.3f} depth {:+.3f}".format(
            spearmanr(share, Q["n50"])[0], spearmanr(share, Q["depth"])[0],
            spearmanr(acc_share, Q["n50"])[0], spearmanr(acc_share, Q["depth"])[0]))
    cat = pd.Series(np.select([f >= .99, f >= .95, f >= .15, f > 0],
                              ["core", "soft_core", "shell", "cloud"], "absent"))
    say("categories (frequency over called cells): " + ", ".join(
        "{} {}".format(k, int(v)) for k, v in cat.value_counts().reindex(
            ["core", "soft_core", "shell", "cloud", "absent"]).fillna(0).items()))

    ena = pd.read_csv(args.ena, sep="\t", dtype=str).drop_duplicates("run_accession").set_index("run_accession")
    study = pd.Series([ena.at[g, "study_accession"] if g in ena.index else "NA" for g in G], index=G)
    say("\nnegative control within studies (significant / tested); compare 35: "
        "union80 115 & 63, binary reads 212 & 260")
    say("{:<13} {:>4}  {:>16}  {:>16}".format("study", "n", "N50 split", "depth split"))
    tot = [0, 0]
    for s, n in study.value_counts().items():
        if n < args.min_study or s == "NA":
            continue
        m = (study == s).values
        cells = []
        for k, q in enumerate((Q["n50"].values[m], Q["depth"].values[m])):
            sg, nt = fisher_count(M[:, m], q < np.median(q))
            tot[k] += sg
            cells.append("{:>6} of {:<6}".format(sg, nt))
        say("{:<13} {:>4}  ".format(s, int(m.sum())) + "  ".join(cells))
    say("{:<13} {:>4}  {:>16}  {:>16}".format("TOTAL", "", tot[0], tot[1]))
    (d / "ternary_report.txt").write_text("\n".join(R) + "\n")
    say("\nwritten: {}/presence_reads_ternary.tsv, ternary_report.txt".format(d))


if __name__ == "__main__":
    main()
