#!/usr/bin/env python3
"""
20_scoary_collect.py

Collect the Scoary runs of 19_scoary_run.sbatch, apply the reported filter,
and test every primary hit for robustness.

Filter (stated in code, default as in the Methods):
    Bonferroni_p < 0.05  AND  Worst_pairwise_comp_p < 0.05
Bonferroni removes chance hits across the genes tested in that run; the worst
pairwise p (Scoary's pairwise comparisons on the phylogeny) removes hits that
are only clonal structure. Example from C_nonB__Region_Americas: group_603 has
naive p 1e-182 but worst pairwise p 0.87 — one lineage, not an association.

Robustness of each PRIMARY hit (subset C_nonB):
  in_<subset>   for every sensitivity run of the same trait:
                  "sig"  passes the same filter, same direction
                  "ns"   tested but not significant (or opposite direction)
                  "-"    not tested (gene invariant there, or run skipped)
  study_consistent  presence frequency of the gene per study (studies with
                >= --min-study isolates in the run and known trait), computed
                separately among trait-positive and trait-negative isolates.
                "yes" if every study's positive-group frequency is on the same
                side of every study's negative-group frequency as the overall
                association; for region traits, where each study is from one
                region, this compares positive-region studies against
                negative-region studies. "insufficient" if < 2 studies on
                either side.
Direction: Odds_ratio > 1 = gene enriched in trait-positive isolates.

Outputs (06_pangwas/scoary_summary/)
  run_summary.tsv     per run: genes tested, naive / Bonferroni / filter hits
  all_hits.tsv        every row Scoary wrote (naive p < 0.05), all runs
  significant.tsv     rows passing the filter, all runs
  primary_hits.tsv    C_nonB hits with robustness columns and per-study
                      frequencies
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

ROOT = Path("/cfs/earth/scratch/tabasran/MAP")
NEED = ["Gene", "Annotation", "Number_pos_present_in", "Number_neg_present_in",
        "Number_pos_not_present_in", "Number_neg_not_present_in", "Odds_ratio",
        "Naive_p", "Bonferroni_p", "Benjamini_H_p", "Max_Pairwise_comparisons",
        "Best_pairwise_comp_p", "Worst_pairwise_comp_p"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default=str(ROOT / "06_pangwas/scoary_runs"))
    ap.add_argument("--ena", default=str(ROOT / "01_metadata/ena_raw.tsv"))
    ap.add_argument("--out", default=str(ROOT / "06_pangwas/scoary_summary"))
    ap.add_argument("--bonf", type=float, default=0.05)
    ap.add_argument("--worst-pw", type=float, default=0.05)
    ap.add_argument("--primary", default="C_nonB")
    ap.add_argument("--min-study", type=int, default=10)
    args = ap.parse_args()

    runs = Path(args.runs)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    man = pd.read_csv(runs / "manifest.tsv", sep="\t")

    frames, summ, missing = [], [], []
    for _, m in man.iterrows():
        f = list((runs / m["run"]).glob("*.results.csv"))
        if len(f) != 1:
            missing.append(m["run"])
            continue
        r = pd.read_csv(f[0])
        lack = [c for c in NEED if c not in r.columns]
        if lack:
            sys.exit("[!] {} lacks columns {}".format(f[0], lack))
        r.insert(0, "run", m["run"])
        r.insert(1, "subset", m["subset"])
        r.insert(2, "trait", m["trait"])
        r["direction"] = (r["Odds_ratio"] > 1).map({True: "enriched", False: "depleted"})
        r["pass"] = (r["Bonferroni_p"] < args.bonf) & (r["Worst_pairwise_comp_p"] < args.worst_pw)
        frames.append(r)
        summ.append({"run": m["run"], "subset": m["subset"], "trait": m["trait"],
                     "pos": m["pos"], "neg": m["neg"], "genes_tested": m["genes_tested"],
                     "naive_p<0.05": len(r),
                     "bonferroni<{}".format(args.bonf): int((r["Bonferroni_p"] < args.bonf).sum()),
                     "PASS(bonf&worstPW)": int(r["pass"].sum()),
                     "pass_enriched": int((r["pass"] & (r["direction"] == "enriched")).sum()),
                     "pass_depleted": int((r["pass"] & (r["direction"] == "depleted")).sum())})
    if missing:
        print("[!] runs without exactly one results.csv ({}): {}".format(len(missing), missing))
    allr = pd.concat(frames, ignore_index=True)
    allr.to_csv(out / "all_hits.tsv", sep="\t", index=False)
    sig = allr[allr["pass"]]
    sig.to_csv(out / "significant.tsv", sep="\t", index=False)
    sm = pd.DataFrame(summ)
    sm.to_csv(out / "run_summary.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    print("=== per run (filter: Bonferroni < {} AND worst pairwise < {}) ===".format(
        args.bonf, args.worst_pw))
    print(sm.drop(columns=["run"]).to_string(index=False))

    # ------------------------------------------------ robustness of primary
    prim = sig[sig["subset"] == args.primary].copy()
    print("\nprimary ({}) hits passing the filter: {}".format(args.primary, len(prim)))
    if prim.empty:
        return
    others = [s for s in man["subset"].unique() if s != args.primary]
    key = allr.set_index(["subset", "trait", "Gene"])
    tested = {}
    for _, m in man.iterrows():
        g = pd.read_csv(runs / m["run"] / "genes.csv", usecols=[0], dtype=str).iloc[:, 0]
        tested[(m["subset"], m["trait"])] = set(g)
    for s in others:
        col = []
        for _, h in prim.iterrows():
            k = (s, h["trait"], h["Gene"])
            if (s, h["trait"]) not in tested or h["Gene"] not in tested[(s, h["trait"])]:
                col.append("-")
            elif k in key.index:
                row = key.loc[k]
                ok = bool(row["pass"]) and row["direction"] == h["direction"]
                col.append("sig" if ok else "ns")
            else:
                col.append("ns")            # tested, naive p >= 0.05
        prim["in_" + s] = col

    # per-study frequencies inside the primary run
    ena = pd.read_csv(args.ena, sep="\t", dtype=str).drop_duplicates("run_accession").set_index("run_accession")
    cons, detail = [], []
    cache = {}
    for _, h in prim.iterrows():
        rd = runs / "{}__{}".format(args.primary, h["trait"])
        if rd not in cache:
            gm = pd.read_csv(rd / "genes.csv", dtype=str, keep_default_na=False).set_index("Gene")
            tr = pd.read_csv(rd / "trait.csv", dtype=str).set_index("Name").iloc[:, 0]
            iso = list(tr.index)
            st = pd.Series([ena.at[i, "study_accession"] if i in ena.index else "assembly"
                            for i in iso], index=iso)
            cache[rd] = (gm, tr, st)
        gm, tr, st = cache[rd]
        p = (gm.loc[h["Gene"], list(tr.index)].str.strip() != "")
        fp, fn = [], []
        for s_, idx in st.groupby(st).groups.items():
            if s_ == "assembly":
                continue
            t = tr.loc[idx]
            for lab, store in (("1", fp), ("0", fn)):
                sel = t.index[t == lab]
                if len(sel) >= args.min_study:
                    store.append((s_, float(p.loc[sel].mean()), len(sel)))
        if len(fp) < 2 or len(fn) < 2:
            c = "insufficient"
        elif h["direction"] == "enriched":
            c = "yes" if min(x[1] for x in fp) > max(x[1] for x in fn) else "no"
        else:
            c = "yes" if max(x[1] for x in fp) < min(x[1] for x in fn) else "no"
        cons.append(c)
        detail.append("pos[" + "; ".join("{}:{:.2f}(n={})".format(*x) for x in sorted(fp, key=lambda x: -x[1])) +
                      "] neg[" + "; ".join("{}:{:.2f}(n={})".format(*x) for x in sorted(fn, key=lambda x: -x[1])) + "]")
    prim["study_consistent"] = cons
    prim["study_frequencies"] = detail

    keep = ["trait", "Gene", "Annotation", "direction", "Number_pos_present_in",
            "Number_pos_not_present_in", "Number_neg_present_in", "Number_neg_not_present_in",
            "Odds_ratio", "Bonferroni_p", "Worst_pairwise_comp_p", "Max_Pairwise_comparisons"] + \
           ["in_" + s for s in others] + ["study_consistent", "study_frequencies"]
    prim = prim[keep].sort_values(["trait", "Bonferroni_p"])
    prim.to_csv(out / "primary_hits.tsv", sep="\t", index=False)

    print("\n=== primary hits: robustness ===")
    for t, g in prim.groupby("trait"):
        line = "  {:<22} n={:<4}".format(t, len(g))
        for s in others:
            c = g["in_" + s].value_counts()
            line += "  {}: sig {} / ns {} / - {}".format(
                s, c.get("sig", 0), c.get("ns", 0), c.get("-", 0))
        c = g["study_consistent"].value_counts()
        line += "  | study-consistent yes {} / no {} / insuff {}".format(
            c.get("yes", 0), c.get("no", 0), c.get("insufficient", 0))
        print(line)
    print("\nwritten: {}/run_summary.tsv, all_hits.tsv, significant.tsv, primary_hits.tsv".format(out))


if __name__ == "__main__":
    main()
