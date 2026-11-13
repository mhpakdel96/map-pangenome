#!/usr/bin/env python3
"""
25_presence_consensus.py

Compare gene-presence calls and pick the one least tied to assembly quality.
No new alignment: reads the outputs of 24_dna_pangenome.py.

What 24 showed
  * annotation calls: genes per genome FALL with N50 (Spearman -0.83) —
    fragmented assemblies carry split/extra gene calls
  * DNA calls (>= 80% covered): genes per genome RISE with N50 (+0.80) —
    in fragmented assemblies part of a gene can be missing from the contigs
    (assembly gap), and minimap2 misses very short genes (clusters that go
    core -> 'absent' have median length ~160 bp)
  So each single call carries its own assembly-quality bias, in opposite
  directions.

Calls compared (for kept clusters = not a redundant fragment, has reference)
  ann            Panaroo annotation
  dna80          DNA covered >= 0.80
  dna50          DNA covered >= 0.50 (tolerates assembly gaps inside a gene)
  union80        ann OR dna80   (absent only if both say absent)
  union50        ann OR dna50
For each: genes per genome vs N50 and vs contigs (Spearman), the same within
the genomes with N50 >= --hq-n50 (where gaps are rarer), category
counts, and Heaps' gamma. Clusters shorter than --short bp are reported
separately because the DNA call is unreliable for them.

Outputs (04_pangenome/dna_presence/)
  consensus_report.txt
  presence_<call>.Rtab   for the call chosen with --write (default union50)
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path("/cfs/earth/scratch/tabasran/MAP")


def category(f):
    return np.select([f >= 0.99, f >= 0.95, f >= 0.15, f > 0],
                     ["core", "soft_core", "shell", "cloud"], "absent")


def gamma(M, perms, rng):
    n = M.shape[1]
    grid = [x for x in (10, 20, 50, 100, 200, 500, 1000, n) if x <= n]
    sizes = np.zeros((perms, len(grid)))
    for p in range(perms):
        order = rng.permutation(n)
        cum = np.zeros(M.shape[0], dtype=bool)
        k = 0
        for j, c in enumerate(order, 1):
            cum |= M[:, c]
            if j == grid[k]:
                sizes[p, k] = cum.sum()
                k += 1
                if k == len(grid):
                    break
    mean = sizes.mean(axis=0)
    return float(np.polyfit(np.log(grid), np.log(mean), 1)[0]), mean[-1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=str(ROOT / "04_pangenome/dna_presence"))
    ap.add_argument("--genes", default=str(ROOT / "04_pangenome/panaroo/gene_presence_absence_roary.csv"))
    ap.add_argument("--qc", default=str(ROOT / "02_qc/genome_qc_master.tsv"))
    ap.add_argument("--hq-n50", type=int, default=30000)
    ap.add_argument("--short", type=int, default=250)
    ap.add_argument("--perm", type=int, default=30)
    ap.add_argument("--write", default="union50")
    args = ap.parse_args()
    d = Path(args.dir)
    R = []

    def say(s=""):
        print(s)
        R.append(s)

    T = pd.read_csv(d / "gene_table.tsv", sep="\t").set_index("gene")
    C = pd.read_csv(d / "coverage.tsv.gz", sep="\t", index_col=0)
    gm = pd.read_csv(args.genes, dtype=str, keep_default_na=False).set_index("Gene")
    iso = list(C.columns)
    kept = T.index[T["kept"]]
    A = (gm.loc[kept, iso] != "").values
    D = C.loc[kept, iso].values
    calls = {"ann": A, "dna80": D >= 0.8, "dna50": D >= 0.5,
             "union80": A | (D >= 0.8), "union50": A | (D >= 0.5)}
    short = (T.loc[kept, "len_bp"] < args.short).values
    qc = pd.read_csv(args.qc, sep="\t").set_index("genome").reindex(iso)
    hq = (qc["n50"] >= args.hq_n50).values
    say("kept clusters {} (of which < {} bp: {}) | genomes {} (N50 >= {}: {})".format(
        len(kept), args.short, int(short.sum()), len(iso), args.hq_n50, int(hq.sum())))

    say("\nSpearman of genes per genome with N50 (all genomes | N50 >= {} only); "
        "closer to 0 = less tied to assembly quality".format(args.hq_n50))
    say("{:<9} {:>16} {:>16} {:>18} {:>12}".format(
        "call", "all genes", "accessory*", "accessory, HQ", "median genes"))
    rng = np.random.RandomState(1)
    rows = {}
    for name, M in calls.items():
        f = M.mean(axis=1)
        acc = category(f) != "core"
        g_all = M.sum(axis=0)
        g_acc = M[acc].sum(axis=0)
        r1 = spearmanr(g_all, qc["n50"])[0]
        r2 = spearmanr(g_acc, qc["n50"])[0]
        r3 = spearmanr(g_acc[hq], qc["n50"].values[hq])[0]
        say("{:<9} {:>+16.3f} {:>+16.3f} {:>+18.3f} {:>12.0f}".format(name, r1, r2, r3, np.median(g_all)))
        rows[name] = M
    say("* accessory = not core (< 99%) under that call")

    say("\ncategory counts (kept clusters; 'absent' = in no genome under that call)")
    cats = ["core", "soft_core", "shell", "cloud", "absent"]
    say("{:<9} ".format("call") + " ".join("{:>9}".format(c) for c in cats) + "   gamma   P(1792)")
    for name, M in rows.items():
        c = pd.Series(category(M.mean(axis=1))).value_counts()
        g, p = gamma(M, args.perm, rng)
        say("{:<9} ".format(name) + " ".join("{:>9}".format(int(c.get(x, 0))) for x in cats) +
            "   {:.4f}  {:.0f}".format(g, p))

    say("\nsame, clusters >= {} bp only".format(args.short))
    for name, M in rows.items():
        Ms = M[~short]
        c = pd.Series(category(Ms.mean(axis=1))).value_counts()
        g, p = gamma(Ms, args.perm, rng)
        say("{:<9} ".format(name) + " ".join("{:>9}".format(int(c.get(x, 0))) for x in cats) +
            "   {:.4f}  {:.0f}".format(g, p))

    M = rows[args.write]
    pd.DataFrame(M.astype(int), index=kept, columns=iso).to_csv(
        d / "presence_{}.Rtab".format(args.write), sep="\t", index_label="Gene")
    (d / "consensus_report.txt").write_text("\n".join(R) + "\n")
    say("\nwritten: {}/presence_{}.Rtab, consensus_report.txt".format(d, args.write))


if __name__ == "__main__":
    main()
