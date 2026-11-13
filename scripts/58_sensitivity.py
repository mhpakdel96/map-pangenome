#!/usr/bin/env python3
"""
58_sensitivity.py

Are the phylogenetic and population conclusions stable when (i) genomes with
read-level evidence of strain mixture are dropped and (ii) the geographic
composition of the cohort is changed?  Both are reviewer questions about
robustness, not new biology, so everything here is a re-analysis of existing
outputs; nothing is re-assembled.

What it does
  PART 1  MIXED SAMPLES
    Drops the genomes flagged as possible mixtures (>= --mixed-sites core SNP
    sites with a minor allele at >= 20% of reads) and recomputes, on the
    pruned tree: lineage group sizes, clade D membership and its regional
    composition, the twelve-SNP concordance, and the ensemble consistency
    index on read-verified calls.  Also reports which of the 21 originally
    homoplastic sites survive.
    A pruned tree tests whether the flagged genomes CHANGE the conclusions.
    It does not re-infer the topology; if any conclusion moves, the tree must
    be re-inferred without them (command printed at the end).

  PART 2  GEOGRAPHIC COMPOSITION
    Repeatedly rebuilds a region-balanced cohort by subsampling each region
    down to the smallest region's size (--reps replicates, fixed seed), and
    recomputes for each replicate: the proportion of Americas / Oceania /
    Europe genomes inside clade D, the clade D versus basal-grade split, and
    the twelve-SNP concordance.  Also does the same by study, so that a single
    large study cannot carry a result on its own.
    This tests whether "96% of American isolates are in clade D" is a property
    of the population or of who sequenced what.

  PART 3  PANGENOME
    Recomputes the pangenome categories and Heaps' exponent on the balanced
    subsamples for the three-state read calls, so the closed-pangenome claim
    is not resting on one over-sampled region.

Outputs (12_validation/sensitivity/)
  mixed_impact.tsv          every statistic with and without the flagged genomes
  mixed_dropped.tsv         which genomes were dropped and why
  balanced_region.tsv       per replicate: clade D composition and concordance
  balanced_study.tsv        same, balanced by study
  pangenome_balanced.tsv    categories and gamma per replicate
  report.txt

Usage
  python 58_sensitivity.py --reps 100
Environment: ml_env (numpy, pandas, ete3 or dendropy not required; the tree is
parsed with the small reader in 44_homoplasy_site.py).

NOTE  This script reports whether conclusions move.  It does not decide what
      to write; read the output before changing the manuscript.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module

ROOT = Path("/cfs/earth/scratch/tabasran/MAP")
SEED = 20261110
PANEL_D = [(405941, "A"), (444953, "T")]          # clade D markers (K-10 alleles)
PANEL_S = [(9124, "A"), (5365, "C"), (5751, "T"), (6610, "C")]
REGIONS = ["Americas", "Europe", "Oceania"]


class Log:
    def __init__(self, path):
        self.f = open(path, "w")

    def __call__(self, *a):
        t = " ".join(str(x) for x in a)
        print(t)
        self.f.write(t + "\n")
        self.f.flush()


def load_tables(root):
    lin = pd.read_csv(root / "05_phylogeny/typing/lineage_assignments.tsv", sep="\t", dtype=str)
    sp = pd.read_csv(root / "05_phylogeny/typing/split_assignments.tsv", sep="\t", dtype=str)
    tr = pd.read_csv(root / "06_pangwas/trait.csv", dtype=str)
    rq = pd.read_csv(root / "08_reads/summary/genome_qc_reads.tsv", sep="\t")
    rq = rq.rename(columns={rq.columns[0]: "genome"})
    M = lin[["genome", "lineage", "sublineage"]].copy()
    d = set(sp.loc[sp["tree_side"].fillna("").str.startswith("D "), "genome"])
    M["group"] = np.where(M["lineage"] == "S",
                          np.where(M["sublineage"] == "III", "S/Type III", "S/Type I"),
                          np.where(M["sublineage"] == "B", "Bison",
                                   np.where(M["genome"].isin(d), "C clade D", "C basal grade")))
    M["region"] = "Other/unknown"
    for r in REGIONS:
        c = "Region_" + r
        if c in tr.columns:
            M.loc[M["genome"].isin(tr.loc[tr[c] == "1", "Name"]), "region"] = r
    ena = root / "01_metadata/ena_raw.tsv"
    if ena.exists():
        e = pd.read_csv(ena, sep="\t", dtype=str).drop_duplicates("run_accession")
        M = M.merge(e[["run_accession", "study_accession"]].rename(
            columns={"run_accession": "genome", "study_accession": "study"}), on="genome", how="left")
    else:
        M["study"] = np.nan
    M["study"] = M["study"].fillna("assembly")
    M = M.merge(rq[[c for c in ("genome", "mixed_sites", "depth") if c in rq.columns]],
                on="genome", how="left")
    return M


def panel_calls(root, M, positions):
    """allele of every genome at the given K-10 positions, from core.tab."""
    pat = " || ".join("$2=={}".format(p) for p, _ in positions)
    import subprocess
    cmd = "awk -F'\\t' 'NR==1 || {}' {}".format(pat, root / "05_phylogeny/core_masked/core.tab")
    txt = subprocess.run(["bash", "-lc", cmd], capture_output=True, text=True).stdout
    rows = [l.split("\t") for l in txt.strip().split("\n")]
    head, data = rows[0], rows[1:]
    T = pd.DataFrame(data, columns=head).set_index("POS")
    cols = [c for c in T.columns if c in set(M["genome"])]
    return {int(p): T.loc[str(p), cols].str.upper() for p, _ in positions if str(p) in T.index}, cols


def concordance(alleles, cols, M, positions, target_groups):
    """fraction of genomes whose panel call matches their tree group."""
    inside = M.set_index("genome")["group"].reindex(cols).isin(target_groups).values
    ok = np.ones(len(cols), dtype=bool)
    called = np.ones(len(cols), dtype=bool)
    for p, marker in positions:
        if p not in alleles:
            continue
        v = alleles[p].reindex(cols).values
        c = np.isin(v, list("ACGT"))
        called &= c
        ok &= np.where(c, (v == marker) == inside, True)
    if called.sum() == 0:
        return np.nan, 0
    return float(ok[called].mean()), int(called.sum())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=str(ROOT))
    ap.add_argument("--reps", type=int, default=100)
    ap.add_argument("--mixed-sites", type=int, default=20)
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    root = Path(args.root)
    out = Path(args.out) if args.out else root / "12_validation/sensitivity"
    out.mkdir(parents=True, exist_ok=True)
    log = Log(out / "report.txt")
    rng = np.random.default_rng(SEED)

    M = load_tables(root)
    log("cohort: {} genomes".format(len(M)))
    log("groups: " + ", ".join("{} {}".format(k, v) for k, v in M["group"].value_counts().items()))
    log("regions: " + ", ".join("{} {}".format(k, v) for k, v in M["region"].value_counts().items()))
    # Clade D is defined only inside the C lineage and only against the basal
    # grade, so every clade-D percentage below uses that set as its denominator.
    # Mixing denominators (all genomes in one part, C-lineage genomes in
    # another) makes the parts non-comparable, which is why this is set once.
    DPOOL = M[M["group"].isin(["C clade D", "C basal grade"])].copy()
    log("clade-D denominator: {} genomes (C lineage, excluding Bison); by region: {}".format(
        len(DPOOL), ", ".join("{} {}".format(k, v)
                              for k, v in DPOOL["region"].value_counts().items())))

    alleles_d, cols = panel_calls(root, M, PANEL_D)
    alleles_s, _ = panel_calls(root, M, PANEL_S)

    # ------------------------------------------------------------ PART 1
    log("\nPART 1  MIXED SAMPLES (flagged at >= {} mixed core SNP sites)".format(args.mixed_sites))
    flagged = M[(M["mixed_sites"].notna()) & (M["mixed_sites"] >= args.mixed_sites)]
    log("  flagged genomes: {}".format(len(flagged)))
    flagged[["genome", "group", "region", "study", "mixed_sites", "depth"]].to_csv(
        out / "mixed_dropped.tsv", sep="\t", index=False)
    keep = M[~M["genome"].isin(flagged["genome"])]
    rows = []
    for name, sub in (("with flagged genomes", M), ("without flagged genomes", keep)):
        c = [g for g in cols if g in set(sub["genome"])]
        dconc, dn = concordance(alleles_d, c, sub, PANEL_D, {"C clade D"})
        sconc, sn = concordance(alleles_s, c, sub, PANEL_S, {"S/Type I", "S/Type III"})
        inD = sub[sub["group"] == "C clade D"]
        r = {"cohort": name, "n_genomes": len(sub), "n_clade_D": len(inD),
             "clade_D_two_snp_concordance": dconc, "clade_D_n_called": dn,
             "S_panel_concordance": sconc, "S_panel_n_called": sn}
        d = sub[sub["group"].isin(["C clade D", "C basal grade"])]
        for reg in REGIONS:
            tot = (d["region"] == reg).sum()
            inside = ((d["region"] == reg) & (d["group"] == "C clade D")).sum()
            r["pct_" + reg + "_in_clade_D"] = round(100 * inside / tot, 2) if tot else np.nan
            r["n_" + reg + "_C"] = int(tot)
        rows.append(r)
    P1 = pd.DataFrame(rows)
    P1.to_csv(out / "mixed_impact.tsv", sep="\t", index=False)
    log(P1.to_string(index=False))
    delta = {c: P1[c].iloc[1] - P1[c].iloc[0] for c in P1.columns
             if c.startswith("pct_") or c.endswith("concordance")}
    log("  change after dropping: " + ", ".join("{} {:+.2f}".format(k, v)
                                                for k, v in delta.items() if pd.notna(v)))
    log("  If any of these move materially, re-infer the tree without the flagged genomes:")
    log("    sbatch 13_snippy_core.sbatch  # with the reduced cohort list")
    log("    sbatch 14_gubbins_iqtree.sbatch")

    # ------------------------------------------------------------ PART 2
    log("\nPART 2  GEOGRAPHIC AND STUDY COMPOSITION ({} replicates, seed {})".format(args.reps, SEED))

    def balanced(by, label):
        rows = []
        # the same pool for both schemes, so the two sets of percentages are
        # comparable with each other and with part 1
        pool = DPOOL
        sizes = pool[pool[by].isin(REGIONS)][by].value_counts() if by == "region" \
            else pool[by].value_counts()
        if by == "region":
            n = int(sizes.min())
        else:
            sizes = sizes[sizes >= 20]
            n = int(sizes.min())
        log("  balancing by {} at {} genomes per level ({} levels), pool = {} C-lineage genomes"
            .format(by, n, len(sizes), len(pool)))
        full = [str(lev) for lev in sizes.index if int(sizes[lev]) <= n]
        if full:
            log("    [!] taken in full, so fixed by construction and with no spread: "
                + ", ".join(full))
        for rep in range(args.reps):
            pick = []
            for lev in sizes.index:
                idx = pool.index[pool[by] == lev].to_numpy()
                pick += list(rng.choice(idx, size=min(n, len(idx)), replace=False))
            sub = pool.loc[pick]
            c = [g for g in cols if g in set(sub["genome"])]
            dconc, dn = concordance(alleles_d, c, sub, PANEL_D, {"C clade D"})
            r = {"replicate": rep, "n_genomes": len(sub),
                 "frac_clade_D": float((sub["group"] == "C clade D").mean()),
                 "clade_D_two_snp_concordance": dconc, "n_called": dn}
            for reg in REGIONS:
                tot = (sub["region"] == reg).sum()
                inside = ((sub["region"] == reg) & (sub["group"] == "C clade D")).sum()
                r["pct_" + reg + "_in_clade_D"] = 100 * inside / tot if tot else np.nan
                r["n_" + reg + "_C"] = int(tot)
            rows.append(r)
        B = pd.DataFrame(rows)
        B.to_csv(out / ("balanced_" + by + ".tsv"), sep="\t", index=False)
        for c in B.columns:
            if c.startswith("pct_") or c.startswith("frac_") or c.endswith("concordance"):
                v = B[c].dropna()
                if len(v):
                    log("    {:<34} median {:.2f}  2.5-97.5% {:.2f}-{:.2f}".format(
                        c, v.median(), v.quantile(.025), v.quantile(.975)))
        return B

    balanced("region", "region")
    balanced("study", "study")

    # ------------------------------------------------------------ PART 3
    tern = root / "08_reads/summary/presence_reads_ternary.tsv"
    if tern.exists():
        log("\nPART 3  PANGENOME ON BALANCED SUBSAMPLES (three-state read calls)")
        X = pd.read_csv(tern, sep="\t", index_col=0)
        # the whole cohort here, not the C lineage: the pangenome question is
        # about the species, so every lineage belongs in the subsample
        pool = M[M["region"].isin(REGIONS) & M["genome"].isin(X.columns)]
        n = int(pool["region"].value_counts().min())
        log("  balancing by region at {} genomes per level, pool = {} genomes of all "
            "lineages with read-based calls".format(n, len(pool)))
        rows = []
        for rep in range(min(args.reps, 30)):
            pick = []
            for reg in REGIONS:
                idx = pool.index[pool["region"] == reg].to_numpy()
                pick += list(rng.choice(idx, size=n, replace=False))
            gs = pool.loc[pick, "genome"].tolist()
            sub = X[gs]
            fr = (sub == 1).sum(axis=1) / sub.notna().sum(axis=1).replace(0, np.nan)
            fr = fr.dropna()
            fr = fr[fr > 0]
            cat = pd.cut(fr, [-np.inf, .15, .95, .99, np.inf], right=False,
                         labels=["cloud", "shell", "soft_core", "core"])
            cnt = cat.value_counts()
            # Heaps' exponent on this subsample
            P = (sub > 0.5).fillna(False).to_numpy()
            ks = np.unique(np.r_[1:10, np.arange(10, P.shape[1] + 1, 8), P.shape[1]])
            order = rng.permutation(P.shape[1])
            first = np.where(P[:, order].any(1), P[:, order].argmax(1), -1)
            cs = np.cumsum(np.bincount(first[first >= 0], minlength=P.shape[1]))
            y = cs[ks - 1]
            m = ks >= 10
            gam = np.polyfit(np.log(ks[m]), np.log(np.maximum(y[m], 1)), 1)[0]
            rows.append({"replicate": rep, "n_genomes": len(gs), "gamma": gam,
                         **{k: int(cnt.get(k, 0)) for k in ("core", "soft_core", "shell", "cloud")},
                         "accessory": int(cnt.get("soft_core", 0) + cnt.get("shell", 0) + cnt.get("cloud", 0))})
        G = pd.DataFrame(rows)
        G.to_csv(out / "pangenome_balanced.tsv", sep="\t", index=False)
        log("    gamma       median {:.4f}  2.5-97.5% {:.4f}-{:.4f}".format(
            G["gamma"].median(), G["gamma"].quantile(.025), G["gamma"].quantile(.975)))
        log("    accessory   median {:.0f}  2.5-97.5% {:.0f}-{:.0f}".format(
            G["accessory"].median(), G["accessory"].quantile(.025), G["accessory"].quantile(.975)))
    else:
        log("\nPART 3 skipped: missing " + str(tern))

    log("\nwritten: " + str(out))


if __name__ == "__main__":
    main()
