#!/usr/bin/env python3
"""
21_batch_check.py

Are the primary Scoary hits biology, or sequencing/study batch?

Why this is needed. In the C_nonB runs, 118 Region_Americas and 80
Region_Europe hits pass Bonferroni + worst-pairwise, but almost none survive
in the MiSeq-only subset and none are study-consistent. Scoary's pairwise
test controls for clonal structure, NOT for technical batch: a gene that is
missed by assembly/annotation in one sequencing batch is absent on many
unrelated branches, which looks like independent gains/losses and passes the
pairwise test. Platform and region are almost fully confounded (NovaSeq 6000
and HiSeq 2500 = Americas only; NextSeq 500 = Europe only), so this has to be
tested directly.

For every primary hit (06_pangwas/scoary_summary/primary_hits.tsv), within the
isolates of its C_nonB run:

  delta_all       freq(trait+) - freq(trait-)
  delta_<platform>  the same difference computed INSIDE one platform, for
                  every platform with >= --min-n positives AND negatives.
                  A biological association should keep its sign and roughly
                  its size inside a platform.
  pos_platform_range / neg_platform_range
                  among trait-positive (or -negative) isolates only: max-min
                  gene frequency across platforms with >= --min-n isolates,
                  and the chi-square p of presence x platform. If the gene
                  varies between platforms WITHIN one trait group as much as
                  between the trait groups, platform explains the hit.
  n50_carrier / n50_noncarrier   median N50 of carriers and non-carriers

Verdict per hit
  within_platform_supported  every eligible platform has the same sign and
                  |delta| >= 0.5 |delta_all| (needs >= 1 eligible platform)
  batch_suspect   max(pos_platform_range, neg_platform_range) >= |delta_all|
  untestable      no platform with enough positives and negatives AND no
                  within-group platform contrast

Outputs (06_pangwas/scoary_summary/)
  batch_check.tsv           one row per primary hit
  batch_check_summary.txt   counts per trait
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

try:
    from scipy.stats import chi2_contingency
except ImportError:
    sys.exit("[!] scipy is required (conda install scipy in this env)")

ROOT = Path("/cfs/earth/scratch/tabasran/MAP")


def chi_p(presence, groups):
    tab = pd.crosstab(groups, presence)
    if tab.shape[0] < 2 or tab.shape[1] < 2:
        return float("nan")
    return float(chi2_contingency(tab.values)[1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default=str(ROOT / "06_pangwas/scoary_runs"))
    ap.add_argument("--summary", default=str(ROOT / "06_pangwas/scoary_summary"))
    ap.add_argument("--ena", default=str(ROOT / "01_metadata/ena_raw.tsv"))
    ap.add_argument("--qc", default=str(ROOT / "02_qc/genome_qc_master.tsv"))
    ap.add_argument("--primary", default="C_nonB")
    ap.add_argument("--min-n", type=int, default=20)
    args = ap.parse_args()

    summ = Path(args.summary)
    hits = pd.read_csv(summ / "primary_hits.tsv", sep="\t")
    ena = pd.read_csv(args.ena, sep="\t", dtype=str).drop_duplicates("run_accession").set_index("run_accession")
    qc = pd.read_csv(args.qc, sep="\t").set_index("genome")

    rows = []
    for trait, hs in hits.groupby("trait"):
        rd = Path(args.runs) / "{}__{}".format(args.primary, trait)
        gm = pd.read_csv(rd / "genes.csv", dtype=str, keep_default_na=False).set_index("Gene")
        tr = pd.read_csv(rd / "trait.csv", dtype=str).set_index("Name").iloc[:, 0].astype(int)
        iso = list(tr.index)
        plat = pd.Series([ena.at[i, "instrument_model"] if i in ena.index else "assembly"
                          for i in iso], index=iso)
        n50 = qc.reindex(iso)["n50"]
        # platforms usable for each test
        elig = [p for p in plat.unique() if p != "assembly"
                and (tr[plat == p] == 1).sum() >= args.min_n
                and (tr[plat == p] == 0).sum() >= args.min_n]
        for _, h in hs.iterrows():
            pres = (gm.loc[h["Gene"], iso].str.strip() != "").astype(int)
            d_all = pres[tr == 1].mean() - pres[tr == 0].mean()
            r = {"trait": trait, "Gene": h["Gene"], "Annotation": h["Annotation"],
                 "direction": h["direction"], "delta_all": round(d_all, 3)}
            sup = []
            for p in sorted(elig):
                m = plat == p
                d = pres[m & (tr == 1)].mean() - pres[m & (tr == 0)].mean()
                r["delta_" + p.replace(" ", "_")] = round(d, 3)
                sup.append((d * d_all > 0) and abs(d) >= 0.5 * abs(d_all))
            for lab, val in (("pos", 1), ("neg", 0)):
                sel = tr == val
                ps = [p for p in plat[sel].unique() if p != "assembly"
                      and (plat[sel] == p).sum() >= args.min_n]
                if len(ps) >= 2:
                    s2 = sel & plat.isin(ps)
                    f = pres[s2].groupby(plat[s2]).mean()
                    r[lab + "_platform_range"] = round(float(f.max() - f.min()), 3)
                    r[lab + "_platform_p"] = chi_p(pres[s2], plat[s2])
                    r[lab + "_platform_freqs"] = "; ".join(
                        "{}:{:.2f}".format(k, v) for k, v in f.sort_values().items())
                else:
                    r[lab + "_platform_range"] = float("nan")
                    r[lab + "_platform_p"] = float("nan")
                    r[lab + "_platform_freqs"] = ""
            r["n50_carrier"] = float(n50[pres == 1].median()) if pres.sum() else float("nan")
            r["n50_noncarrier"] = float(n50[pres == 0].median()) if (pres == 0).sum() else float("nan")
            rng = max(x for x in (r["pos_platform_range"], r["neg_platform_range"], -1)
                      if x == x)
            r["within_platform_supported"] = ("yes" if all(sup) else "no") if sup else "untestable"
            r["batch_suspect"] = ("yes" if rng >= abs(d_all) else "no") if rng >= 0 else "untestable"
            rows.append(r)

    out = pd.DataFrame(rows)
    out.to_csv(summ / "batch_check.tsv", sep="\t", index=False)

    lines = ["platforms counted only with >= {} isolates per group".format(args.min_n), ""]
    lines.append("{:<22} {:>5}  {:>22}  {:>22}  {:>8}".format(
        "trait", "hits", "within-platform yes/no/-", "batch-suspect yes/no/-", "robust*"))
    for t, g in out.groupby("trait"):
        a = g["within_platform_supported"].value_counts()
        b = g["batch_suspect"].value_counts()
        robust = int(((g["within_platform_supported"] == "yes") & (g["batch_suspect"] != "yes")).sum())
        lines.append("{:<22} {:>5}  {:>22}  {:>22}  {:>8}".format(
            t, len(g), "{}/{}/{}".format(a.get("yes", 0), a.get("no", 0), a.get("untestable", 0)),
            "{}/{}/{}".format(b.get("yes", 0), b.get("no", 0), b.get("untestable", 0)), robust))
    lines.append("\n* robust = supported inside every eligible platform AND not batch-suspect")
    dcols = [c for c in out.columns if c.startswith("delta_")]
    lines.append("\nmedian |delta| by trait (all isolates vs inside each platform):")
    for t, g in out.groupby("trait"):
        lines.append("  {:<22} ".format(t) + "  ".join(
            "{}={:.2f}".format(c[6:], g[c].abs().median()) for c in dcols if g[c].notna().any()))
    lines.append("\nN50 (median) carriers vs non-carriers, by direction:")
    for (t, d), g in out.groupby(["trait", "direction"]):
        lines.append("  {:<22} {:<9} n={:<4} carriers {:>9.0f}   non-carriers {:>9.0f}".format(
            t, d, len(g), g["n50_carrier"].median(), g["n50_noncarrier"].median()))
    txt = "\n".join(lines)
    (summ / "batch_check_summary.txt").write_text(txt + "\n")
    print(txt)
    print("\nwritten: {}/batch_check.tsv, batch_check_summary.txt".format(summ))


if __name__ == "__main__":
    main()
