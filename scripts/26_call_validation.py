#!/usr/bin/env python3
"""
26_call_validation.py

Validate the chosen gene-presence call with a NEGATIVE and a POSITIVE
control, then write the final presence matrix used downstream.

Chosen call (from 25_presence_consensus.py): union80 = present if Panaroo
annotates the gene OR >= 80% of its DNA is in the assembly. It was the only
call whose genes-per-genome count is ~uncorrelated with N50 (Spearman +0.09
all genes, +0.01 accessory, -0.08 within N50 >= 30 kb), while annotation
(-0.83) and DNA-only (+0.80) carry opposite assembly-quality biases.
Clusters < --short bp are set aside: minimap2 cannot detect them, so for
them union80 = annotation and keeps its bias.

NEGATIVE control — assembly quality inside one study.
  For every study with >= --min-study genomes in the cohort, genomes are
  split at that study's median N50 (same population, same lab, same
  platform; only assembly quality differs). Per gene, Fisher's exact test
  low- vs high-N50 half; Bonferroni over genes variable in that study.
  A good call gives ~0 significant genes. Reported for annotation, DNA-only
  and union80.
POSITIVE control — real biology.
  S vs C lineage (whole cohort). Genes with Bonferroni-significant Fisher p
  AND |freq_S - freq_C| >= 0.5. A good call must KEEP these (the call must
  not erase true lineage-specific gene content).

Outputs (04_pangenome/dna_presence/)
  validation_report.txt
  presence_final.Rtab        union80, clusters >= --short bp (the matrix for
                             association and ML)
  presence_final_short.Rtab  union80 for the short clusters (reported, not
                             used for inference)
  lineage_specific_genes.tsv S vs C genes under each call
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

ROOT = Path("/cfs/earth/scratch/tabasran/MAP")


def fisher_scan(M, grp, alpha):
    """M genes x genomes (bool), grp bool per genome. Returns (n tested, n sig, table)."""
    a = M[:, grp].sum(axis=1)
    b = grp.sum() - a
    c = M[:, ~grp].sum(axis=1)
    d = (~grp).sum() - c
    var = (a + c > 0) & (b + d > 0)
    idx = np.where(var)[0]
    p = np.ones(M.shape[0])
    for i in idx:
        p[i] = fisher_exact([[a[i], b[i]], [c[i], d[i]]])[1]
    n = len(idx)
    sig = var & (p < alpha / max(n, 1))
    f1 = a / max(grp.sum(), 1)
    f0 = c / max((~grp).sum(), 1)
    return n, sig, p, f1, f0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=str(ROOT / "04_pangenome/dna_presence"))
    ap.add_argument("--genes", default=str(ROOT / "04_pangenome/panaroo/gene_presence_absence_roary.csv"))
    ap.add_argument("--qc", default=str(ROOT / "02_qc/genome_qc_master.tsv"))
    ap.add_argument("--ena", default=str(ROOT / "01_metadata/ena_raw.tsv"))
    ap.add_argument("--lineages", default=str(ROOT / "05_phylogeny/typing/lineage_assignments.tsv"))
    ap.add_argument("--short", type=int, default=250)
    ap.add_argument("--min-study", type=int, default=60)
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--min-diff", type=float, default=0.5)
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
    longg = kept[(T.loc[kept, "len_bp"] >= args.short).values]
    shortg = kept[(T.loc[kept, "len_bp"] < args.short).values]

    def calls(genes):
        A = (gm.loc[genes, iso] != "").values
        D = C.loc[genes, iso].values >= 0.8
        return {"annotation": A, "dna80": D, "union80": A | D}
    K = calls(longg)
    say("clusters used: {} (>= {} bp); set aside: {} short".format(len(longg), args.short, len(shortg)))

    qc = pd.read_csv(args.qc, sep="\t").set_index("genome").reindex(iso)
    ena = pd.read_csv(args.ena, sep="\t", dtype=str).drop_duplicates("run_accession").set_index("run_accession")
    study = pd.Series([ena.at[i, "study_accession"] if i in ena.index else "assembly" for i in iso], index=iso)

    # ------------------------------------------------ negative control
    say("\nNEGATIVE control: low- vs high-N50 half WITHIN each study "
        "(Fisher, Bonferroni per study); fewer significant genes = better")
    studies = [s for s, n in study.value_counts().items() if n >= args.min_study and s != "assembly"]
    say("{:<14} {:>5} {:>10} {:>10}   ".format("study", "n", "N50 low", "N50 high") +
        "  ".join("{:>14}".format(k) for k in K))
    tot = {k: 0 for k in K}
    for s in studies:
        m = (study == s).values
        n50 = qc["n50"].values[m]
        med = np.median(n50)
        low = n50 < med
        cells = []
        for k, M in K.items():
            n, sig, _, _, _ = fisher_scan(M[:, m], low, args.alpha)
            tot[k] += int(sig.sum())
            cells.append("{:>5} of {:<6}".format(int(sig.sum()), n))
        say("{:<14} {:>5} {:>10.0f} {:>10.0f}   ".format(
            s, int(m.sum()), np.median(n50[low]), np.median(n50[~low])) + "  ".join(cells))
    say("{:<14} {:>5} {:>10} {:>10}   ".format("TOTAL", "", "", "") +
        "  ".join("{:>14}".format(tot[k]) for k in K))

    # ------------------------------------------------ positive control
    lin = pd.read_csv(args.lineages, sep="\t", dtype=str).set_index("genome").reindex(iso)
    S = (lin["lineage"] == "S").values
    say("\nPOSITIVE control: S ({}) vs C ({}) — genes with Bonferroni p < {} and |freq diff| >= {}".format(
        int(S.sum()), int((~S).sum()), args.alpha, args.min_diff))
    rows = []
    sets = {}
    for k, M in K.items():
        n, sig, p, fS, fC = fisher_scan(M, S, args.alpha)
        hit = sig & (np.abs(fS - fC) >= args.min_diff)
        sets[k] = set(longg[hit])
        say("  {:<11} tested {:>5}   significant {:>4}   of which |diff| >= {}: {:>4}  "
            "(S-only/enriched {:>3}, C-only/enriched {:>3})".format(
                k, n, int(sig.sum()), args.min_diff, int(hit.sum()),
                int((hit & (fS > fC)).sum()), int((hit & (fC > fS)).sum())))
        for i in np.where(hit)[0]:
            rows.append({"call": k, "gene": longg[i], "annotation": T.at[longg[i], "annotation"],
                         "len_bp": T.at[longg[i], "len_bp"], "freq_S": round(fS[i], 3),
                         "freq_C": round(fC[i], 3), "p": p[i]})
    pd.DataFrame(rows).to_csv(d / "lineage_specific_genes.tsv", sep="\t", index=False)
    a, u = sets["annotation"], sets["union80"]
    say("  overlap annotation vs union80: both {} | annotation only {} | union80 only {}".format(
        len(a & u), len(a - u), len(u - a)))

    # ------------------------------------------------ final matrices
    pd.DataFrame(K["union80"].astype(int), index=longg, columns=iso).to_csv(
        d / "presence_final.Rtab", sep="\t", index_label="Gene")
    Ks = calls(shortg)
    pd.DataFrame(Ks["union80"].astype(int), index=shortg, columns=iso).to_csv(
        d / "presence_final_short.Rtab", sep="\t", index_label="Gene")
    f = K["union80"].mean(axis=1)
    say("\nfinal matrix presence_final.Rtab: {} clusters x {} genomes | core(>=99%) {} | "
        "soft {} | shell {} | cloud {}".format(
            len(longg), len(iso), int((f >= .99).sum()), int(((f >= .95) & (f < .99)).sum()),
            int(((f >= .15) & (f < .95)).sum()), int(((f > 0) & (f < .15)).sum())))
    (d / "validation_report.txt").write_text("\n".join(R) + "\n")
    say("written: {}/validation_report.txt, presence_final.Rtab, presence_final_short.Rtab, "
        "lineage_specific_genes.tsv".format(d))


if __name__ == "__main__":
    main()
