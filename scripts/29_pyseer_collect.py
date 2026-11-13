#!/usr/bin/env python3
"""
29_pyseer_collect.py

Collect the pyseer runs (28_pyseer_run.sbatch) and apply a pattern-based
Bonferroni threshold per run: alpha / (number of unique variant patterns),
as recommended by the pyseer authors (count_patterns.py computes the same).

Per run: variants tested, unique patterns, threshold, significant variants
(split by direction of beta). The negative-control runs (trait neg_lowN50,
no covariates) show how much residual assembly-quality signal the variant
matrix still carries; a clean matrix gives ~0.

Outputs (06_pangwas/pyseer/)
  run_summary.tsv
  significant.tsv     all significant variants, all runs, with gene
                      annotation for gene variants
"""

import argparse
from pathlib import Path

import pandas as pd

ROOT = Path("/cfs/earth/scratch/tabasran/MAP")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=str(ROOT / "06_pangwas/pyseer"))
    ap.add_argument("--gene-table", default=str(ROOT / "04_pangenome/dna_presence/gene_table.tsv"))
    ap.add_argument("--alpha", type=float, default=0.05)
    args = ap.parse_args()
    d = Path(args.dir)
    man = pd.read_csv(d / "manifest.tsv", sep="\t")
    ann = pd.read_csv(args.gene_table, sep="\t").set_index("gene")["annotation"]

    summ, sig, missing = [], [], []
    for _, m in man.iterrows():
        res, pat = d / "runs" / (m["run"] + ".tsv"), d / "runs" / (m["run"] + ".patterns")
        if not res.exists() or not pat.exists():
            missing.append(m["run"])
            continue
        r = pd.read_csv(res, sep="\t")
        npat = len(set(l.strip() for l in open(pat) if l.strip()))
        thr = args.alpha / max(npat, 1)
        r["significant"] = r["lrt-pvalue"] < thr
        s = r[r["significant"]].copy()
        summ.append({"run": m["run"], "subset": m["subset"], "trait": m["trait"],
                     "variants": m["variants"], "pos": m["pos"], "neg": m["neg"],
                     "covariates": m["covariates"], "tested": len(r), "patterns": npat,
                     "threshold": "{:.2e}".format(thr), "significant": len(s),
                     "beta>0": int((s["beta"] > 0).sum()), "beta<0": int((s["beta"] < 0).sum()),
                     "min_p": "{:.1e}".format(r["lrt-pvalue"].min()) if len(r) else "NA"})
        if len(s):
            s.insert(0, "run", m["run"])
            s.insert(1, "trait", m["trait"])
            s.insert(2, "subset", m["subset"])
            s.insert(3, "variants_type", m["variants"])
            if m["variants"] == "genes":
                s["annotation"] = s["variant"].map(ann)
            sig.append(s)
    S = pd.DataFrame(summ)
    S.to_csv(d / "run_summary.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    print(S.drop(columns=["run"]).to_string(index=False))
    if missing:
        print("\n[!] runs without results ({}): {}".format(len(missing), missing))
    if sig:
        A = pd.concat(sig, ignore_index=True)
        A.to_csv(d / "significant.tsv", sep="\t", index=False)
        print("\nsignificant variants: {} -> {}/significant.tsv".format(len(A), d))
        cols = [c for c in ("trait", "subset", "variants_type", "variant", "af", "lrt-pvalue",
                            "beta", "annotation") if c in A.columns]
        print(A.sort_values("lrt-pvalue")[cols].head(40).to_string(index=False))
    else:
        print("\nno significant variants in any run")


if __name__ == "__main__":
    main()
