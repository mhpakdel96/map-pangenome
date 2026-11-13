#!/usr/bin/env python3
"""
35_read_presence.py

Summarise the read mapping of 34_read_mapping.sbatch.

PART 1  Gene presence from reads (independent of assembly and gene calling)
  present = breadth of read coverage over the cluster >= --min-breadth %.
  Same gene set as the DNA pangenome (kept clusters >= 250 bp, from
  04_pangenome/dna_presence/gene_table.tsv). Genomes with low K-10 depth or
  breadth are flagged (breadth is unreliable at low depth).
  Reported: category counts, Heaps' gamma, agreement with union80, and genes
  per genome vs N50 AND vs read depth (Spearman) — read calls must not
  simply trade the assembly bias for a depth bias.
PART 2  Negative control (as script 26): within every study with >= 60
  genomes, Fisher test low- vs high-N50 half, and low- vs high-DEPTH half;
  Bonferroni per study. Compared for union80 and reads.
PART 3  Mixture index: at core SNP sites with >= --min-dp reads, the minor
  allele fraction (second-largest AD / total). A genome with many sites at
  >= --mixed-af is a candidate mixed infection or contamination.
PART 4  The 21 'homoplastic' SNP sites (31_homoplasy_check.py): allele
  fractions from reads in carriers and non-carriers — are the assembly calls
  clean, or do carriers/non-carriers show mixed reads?

Outputs (08_reads/summary/)
  presence_reads.Rtab   0/1, kept clusters >= 250 bp x genomes with reads
  breadth.tsv.gz        breadth % (same shape)
  genome_qc_reads.tsv   per genome: K-10 depth/breadth, flags, mixture index
  homoplasy_reads.tsv   per site: allele fractions in carriers / non-carriers
  report.txt
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, spearmanr

ROOT = Path("/cfs/earth/scratch/tabasran/MAP")


def category(f):
    return np.select([f >= 0.99, f >= 0.95, f >= 0.15, f > 0],
                     ["core", "soft_core", "shell", "cloud"], "absent")


def gamma(M, perms, rng):
    n = M.shape[1]
    grid = [x for x in (10, 20, 50, 100, 200, 500, 1000, n) if x <= n]
    S = np.zeros((perms, len(grid)))
    for p in range(perms):
        cum = np.zeros(M.shape[0], dtype=bool)
        k = 0
        for j, c in enumerate(rng.permutation(n), 1):
            cum |= M[:, c]
            if j == grid[k]:
                S[p, k] = cum.sum()
                k += 1
                if k == len(grid):
                    break
    return float(np.polyfit(np.log(grid), np.log(S.mean(axis=0)), 1)[0])


def fisher_count(M, grp, alpha=0.05):
    a = M[:, grp].sum(axis=1)
    b = grp.sum() - a
    c = M[:, ~grp].sum(axis=1)
    d = (~grp).sum() - c
    var = (a + c > 0) & (b + d > 0)
    idx = np.where(var)[0]
    n = len(idx)
    sig = 0
    for i in idx:
        if fisher_exact([[a[i], b[i]], [c[i], d[i]]])[1] < alpha / max(n, 1):
            sig += 1
    return sig, n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reads", default=str(ROOT / "08_reads"))
    ap.add_argument("--dna", default=str(ROOT / "04_pangenome/dna_presence"))
    ap.add_argument("--qc", default=str(ROOT / "02_qc/genome_qc_master.tsv"))
    ap.add_argument("--ena", default=str(ROOT / "01_metadata/ena_raw.tsv"))
    ap.add_argument("--snps", default=str(ROOT / "06_pangwas/pyseer/variants_snps.Rtab"))
    ap.add_argument("--homoplasy", default=str(ROOT / "06_pangwas/pyseer/homoplasy_check.tsv"))
    ap.add_argument("--min-breadth", type=float, default=80.0)
    ap.add_argument("--min-depth", type=float, default=15.0, help="flag genomes below this K-10 depth")
    ap.add_argument("--min-dp", type=int, default=10)
    ap.add_argument("--mixed-af", type=float, default=0.2)
    ap.add_argument("--short", type=int, default=250)
    ap.add_argument("--min-study", type=int, default=60)
    ap.add_argument("--perm", type=int, default=30)
    args = ap.parse_args()
    rd = Path(args.reads)
    out = rd / "summary"
    out.mkdir(parents=True, exist_ok=True)
    R = []

    def say(s=""):
        print(s, flush=True)
        R.append(s)

    paths = pd.read_csv(rd / "read_paths.tsv", sep="\t")
    G = [g for g in paths["genome"] if (rd / "pan_cov" / (g + ".tsv")).exists()
         and (rd / "k10" / (g + ".depth")).exists()]
    say("genomes with read results: {} of {} in read_paths.tsv".format(len(G), len(paths)))

    T = pd.read_csv(Path(args.dna) / "gene_table.tsv", sep="\t").set_index("gene")
    kept = T.index[T["kept"] & (T["len_bp"] >= args.short)]

    # ------------------------------------------------ load
    B = {}
    dep = {}
    for k, g in enumerate(G, 1):
        c = pd.read_csv(rd / "pan_cov" / (g + ".tsv"), sep="\t", usecols=[0, 5],
                        names=["rname", "coverage"], header=0)
        B[g] = c.set_index("rname")["coverage"]
        d = open(rd / "k10" / (g + ".depth")).read().split()
        dep[g] = (float(d[0]), float(d[1]))
        if k % 300 == 0:
            print("  loaded {} / {}".format(k, len(G)), flush=True)
    Bm = pd.DataFrame(B).reindex(kept).fillna(0.0)
    Bm.to_csv(out / "breadth.tsv.gz", sep="\t", float_format="%.1f")
    P = (Bm.values >= args.min_breadth)
    pd.DataFrame(P.astype(int), index=kept, columns=G).to_csv(out / "presence_reads.Rtab", sep="\t", index_label="Gene")

    qc = pd.read_csv(args.qc, sep="\t").set_index("genome").reindex(G)
    D = pd.DataFrame(dep, index=["depth", "breadth"]).T.reindex(G)
    low = (D["depth"] < args.min_depth) | (D["breadth"] < 95)
    say("K-10 read depth: median {:.1f}, p5 {:.1f}, p95 {:.1f}; flagged (depth < {} or breadth < 95%): {}".format(
        D["depth"].median(), D["depth"].quantile(.05), D["depth"].quantile(.95), args.min_depth, int(low.sum())))

    # ------------------------------------------------ PART 1
    U = pd.read_csv(Path(args.dna) / "presence_final.Rtab", sep="\t", index_col=0).reindex(index=kept, columns=G)
    Uv = U.values.astype(bool)
    say("\nPART 1 — presence calls on {} clusters (>= {} bp), {} genomes".format(len(kept), args.short, len(G)))
    cats = ["core", "soft_core", "shell", "cloud", "absent"]
    rng = np.random.RandomState(1)
    say("{:<26} ".format("call") + " ".join("{:>9}".format(c) for c in cats) + "    gamma")
    for lab, M in (("union80 (assembly)", Uv), ("reads", P), ("reads, unflagged genomes", P[:, ~low.values])):
        c = pd.Series(category(M.mean(axis=1))).value_counts()
        say("{:<26} ".format(lab) + " ".join("{:>9}".format(int(c.get(x, 0))) for x in cats) +
            "    {:.4f}".format(gamma(M, args.perm, rng)))
    agree = (Uv == P)
    say("agreement union80 vs reads: {:.2f}% of cells | union80 1 & reads 0: {:.2f}% | union80 0 & reads 1: {:.2f}%".format(
        100 * agree.mean(), 100 * (Uv & ~P).mean(), 100 * (~Uv & P).mean()))
    say("genes per genome vs assembly/read quality (Spearman):")
    for lab, M in (("union80", Uv), ("reads", P)):
        f = M.mean(axis=1)
        acc = (f < 0.99) & (f > 0)
        n_all, n_acc = M.sum(axis=0), M[acc].sum(axis=0)
        say("  {:<8} all genes: N50 {:+.3f}  depth {:+.3f} | accessory: N50 {:+.3f}  depth {:+.3f}".format(
            lab, spearmanr(n_all, qc["n50"])[0], spearmanr(n_all, D["depth"])[0],
            spearmanr(n_acc, qc["n50"])[0], spearmanr(n_acc, D["depth"])[0]))

    # ------------------------------------------------ PART 2
    ena = pd.read_csv(args.ena, sep="\t", dtype=str).drop_duplicates("run_accession").set_index("run_accession")
    study = pd.Series([ena.at[g, "study_accession"] if g in ena.index else "NA" for g in G], index=G)
    say("\nPART 2 — negative control within studies (Fisher, Bonferroni per study); significant / tested")
    say("{:<13} {:>4}   {:>18} {:>18}   {:>18} {:>18}".format(
        "study", "n", "N50 split: union80", "reads", "depth split: union80", "reads"))
    tot = np.zeros(4, int)
    for s, n in study.value_counts().items():
        if n < args.min_study or s == "NA":
            continue
        m = (study == s).values
        cells = []
        for q in (qc["n50"].values[m], D["depth"].values[m]):
            g_low = q < np.median(q)
            for M in (Uv, P):
                cells.append(fisher_count(M[:, m], g_low))
        for i, (sg, _) in enumerate(cells):
            tot[i] += sg
        say("{:<13} {:>4}   ".format(s, int(m.sum())) + "  ".join(
            "{:>8} of {:<7}".format(sg, nt) for sg, nt in cells))
    say("{:<13} {:>4}   ".format("TOTAL", "") + "  ".join("{:>18}".format(t) for t in tot))

    # ------------------------------------------------ PART 3 + 4
    say("\nPART 3 — mixture index at core SNP sites (depth >= {}, minor allele >= {})".format(
        args.min_dp, args.mixed_af))
    hom = pd.read_csv(args.homoplasy, sep="\t") if Path(args.homoplasy).exists() else None
    V = pd.read_csv(args.snps, sep="\t", index_col=0)
    hsites = {}
    if hom is not None:
        for v in hom["variant"].unique():
            p = v.split("_")
            hsites[int(p[-3])] = (v, p[-2], p[-1])
    mix_rows, hom_af = [], {v: {} for v, _, _ in hsites.values()}
    for k, g in enumerate(G, 1):
        n_ok = n_mix = 0
        for line in open(rd / "k10" / (g + ".ad.tsv")):
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            ad = [int(x) for x in f[3].split(",") if x not in (".", "")]
            tot_ = sum(ad)
            pos = int(f[0])
            if pos in hsites and tot_ >= args.min_dp:
                v, ref, alt = hsites[pos]
                alts = f[2].split(",")
                a = ad[1 + alts.index(alt)] if alt in alts and 1 + alts.index(alt) < len(ad) else 0
                hom_af[v][g] = a / tot_
            if tot_ < args.min_dp:
                continue
            n_ok += 1
            s = sorted(ad, reverse=True)
            if len(s) > 1 and s[1] / tot_ >= args.mixed_af:
                n_mix += 1
        mix_rows.append({"genome": g, "sites_dp_ok": n_ok, "mixed_sites": n_mix,
                         "mixed_frac": round(n_mix / n_ok, 5) if n_ok else np.nan})
        if k % 300 == 0:
            print("  AD parsed {} / {}".format(k, len(G)), flush=True)
    MX = pd.DataFrame(mix_rows).set_index("genome")
    Q = pd.concat([D, MX], axis=1)
    Q["flag_low_depth"] = low
    Q["n50"] = qc["n50"]
    Q["study"] = study
    Q.to_csv(out / "genome_qc_reads.tsv", sep="\t")
    say("mixed sites per genome: median {:.0f}, p95 {:.0f}, max {:.0f}".format(
        Q["mixed_sites"].median(), Q["mixed_sites"].quantile(.95), Q["mixed_sites"].max()))
    for thr in (5, 20, 50, 200):
        say("  genomes with >= {:>3} mixed sites: {}".format(thr, int((Q["mixed_sites"] >= thr).sum())))
    say("  top 10:\n" + Q.sort_values("mixed_sites", ascending=False).head(10)[
        ["mixed_sites", "sites_dp_ok", "depth", "n50", "study"]].to_string())

    if hom is not None:
        say("\nPART 4 — homoplastic SNP sites: ALT allele fraction from reads")
        rows = []
        for v, d_ in hom_af.items():
            x = V.loc[v]
            af = pd.Series(d_)
            car = [g for g in af.index if x.get(g, 0) == 1]
            non = [g for g in af.index if x.get(g, 0) == 0]
            a_c, a_n = af[car], af[non]
            rows.append({"variant": v, "alt_carriers_called": len(car),
                         "median_af_carriers": round(float(a_c.median()), 3) if len(car) else np.nan,
                         "carriers_mixed(0.2-0.8)": int(((a_c >= .2) & (a_c <= .8)).sum()),
                         "carriers_af<0.2": int((a_c < .2).sum()),
                         "noncarriers_af>=0.2": int((a_n >= .2).sum()),
                         "noncarriers_af>=0.8": int((a_n >= .8).sum()),
                         "locus": hom.loc[hom["variant"] == v, "locus"].iloc[0]})
        H = pd.DataFrame(rows)
        H.to_csv(out / "homoplasy_reads.tsv", sep="\t", index=False)
        pd.set_option("display.width", 250)
        pd.set_option("display.max_colwidth", 60)
        say(H.to_string(index=False))
    (out / "report.txt").write_text("\n".join(R) + "\n")
    say("\nwritten: {}".format(out))


if __name__ == "__main__":
    main()
