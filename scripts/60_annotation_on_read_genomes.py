#!/usr/bin/env python3
"""
60_annotation_on_read_genomes.py

Recompute the ANNOTATION presence call on the 1,747 genomes that have reads,
so that Table 1 of the manuscript can rest on a single genome set instead of
mixing 1,792-genome and 1,747-genome rows.

Prints core / soft core / shell / cloud / accessory and Heaps' gamma for the
annotation call on three cluster sets, so the one that matches Figure 2B
(accessory 1,711) can be identified and used.

Nothing is written to the project tree; the output is text only.

Usage
  conda activate ml_env
  python 60_annotation_on_read_genomes.py
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/cfs/earth/scratch/tabasran/MAP")
SEED = 20261110


def find(*cands):
    for c in cands:
        p = ROOT / c
        if p.exists():
            return p
    return None


def categories(M):
    """M: boolean DataFrame, rows = clusters, columns = genomes."""
    n = M.shape[1]
    freq = M.sum(axis=1) / n
    core = int((freq >= 0.99).sum())
    soft = int(((freq >= 0.95) & (freq < 0.99)).sum())
    shell = int(((freq >= 0.15) & (freq < 0.95)).sum())
    cloud = int(((freq > 0) & (freq < 0.15)).sum())
    absent = int((freq == 0).sum())
    return core, soft, shell, cloud, absent


def heaps_gamma(M, n_perm=50, start=10, seed=SEED):
    """Mean over n_perm random genome orders; fit log P on log n from `start`."""
    rng = np.random.default_rng(seed)
    A = M.to_numpy(dtype=bool)
    n_gen = A.shape[1]
    curves = np.zeros((n_perm, n_gen), dtype=np.int64)
    for i in range(n_perm):
        order = rng.permutation(n_gen)
        seen = np.zeros(A.shape[0], dtype=bool)
        for j, g in enumerate(order):
            seen |= A[:, g]
            curves[i, j] = seen.sum()
    mean_curve = curves.mean(axis=0)
    x = np.log(np.arange(1, n_gen + 1)[start - 1:])
    y = np.log(mean_curve[start - 1:])
    slope = np.polyfit(x, y, 1)[0]
    return slope, mean_curve[-1]


def main():
    rtab = find("04_pangenome/panaroo/gene_presence_absence.Rtab",
                "04_pangenome/panaroo/results/gene_presence_absence.Rtab")
    if rtab is None:
        sys.exit("[!] Panaroo gene_presence_absence.Rtab not found under "
                 + str(ROOT / "04_pangenome/panaroo"))
    reads = find("08_reads/read_paths.tsv")
    if reads is None:
        sys.exit("[!] 08_reads/read_paths.tsv not found")

    print("Panaroo matrix :", rtab)
    print("read genomes   :", reads)

    P = pd.read_csv(rtab, sep="\t", index_col=0)
    P = P.astype(bool)
    print("\nPanaroo matrix: {} clusters x {} genomes".format(*P.shape))

    # genome ids that have reads
    R = pd.read_csv(reads, sep="\t", header=None, dtype=str)
    first = R.iloc[:, 0].tolist()
    if first and first[0].lower() in ("genome", "sample", "id", "accession"):
        first = first[1:]
    want = {Path(g).name.split(".")[0] for g in first}

    cols = {c: Path(c).name.split(".")[0] for c in P.columns}
    keep = [c for c, s in cols.items() if s in want]
    print("read-genome ids in read_paths.tsv : {}".format(len(want)))
    print("matched to Panaroo columns        : {}".format(len(keep)))
    if len(keep) < 1700:
        print("\n[!] too few matched - the naming differs between the two files.")
        print("    first 5 Panaroo columns :", list(P.columns[:5]))
        print("    first 5 read ids        :", first[:5])
        sys.exit(1)

    Q = P[keep]

    # cluster sets to try
    sets = {"all Panaroo clusters": P.index}
    final = find("04_pangenome/dna_presence/presence_final.Rtab")
    if final is not None:
        idx = pd.read_csv(final, sep="\t", index_col=0, usecols=[0]).index
        sets[">=250 bp ({} clusters)".format(len(idx))] = P.index.intersection(idx)
    short = find("04_pangenome/dna_presence/presence_final_short.Rtab")
    if short is not None:
        idx2 = pd.read_csv(short, sep="\t", index_col=0, usecols=[0]).index
        kept = P.index.intersection(idx.union(idx2)) if final is not None else None
        if kept is not None:
            sets["kept ({} clusters)".format(len(kept))] = kept

    hdr = "{:<30} {:>8} {:>7} {:>7} {:>7} {:>7} {:>10} {:>8}"
    row = "{:<30} {:>8} {:>7} {:>7} {:>7} {:>7} {:>10} {:>8.4f}"
    for label, which in sets.items():
        sub = Q.loc[which]
        c, s, sh, cl, ab = categories(sub)
        g, total = heaps_gamma(sub)
        print("\n=== annotation call, {} genomes, {} ===".format(sub.shape[1], label))
        print(hdr.format("", "clusters", "core", "soft", "shell", "cloud", "accessory", "gamma"))
        print(row.format("", sub.shape[0], c, s, sh, cl, s + sh + cl, g))
        if ab:
            print("    ({} clusters present in no genome of this subset)".format(ab))

    print("\nThe row whose accessory equals 1,711 is the one Figure 2B plots.")


if __name__ == "__main__":
    main()
