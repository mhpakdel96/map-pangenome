#!/usr/bin/env python3
"""
30_snp_hits.py

What are the significant core-SNP associations from pyseer?

Context. The negative control (neg_lowN50, 29_pyseer_collect.py) gave 1
significant core SNP but ~90 significant genes: the SNP matrix is clean, the
gene matrix still carries assembly-quality signal. So the interpretable
results are the SNPs. The traits have h^2 = 1.00 (region, host) or 0.5-0.9
(source): phenotype largely predicted by relatedness. A significant SNP can
then be
  (a) a CLADE MARKER — arose once, carried by one clade that happens to be
      concentrated in one region/study (phylogeography, not adaptation), or
  (b) HOMOPLASTIC — arose several times independently in trait-positive
      isolates (convergence; the signature of a possibly adaptive variant).
This script tells them apart and annotates each SNP.

Per significant SNP (per run):
  origins      number of maximal monophyletic groups of carriers on the ML
               tree rooted on the S|C branch (1 = clade marker, >1 = arose
               repeatedly, assuming no reversions)
  largest_clade_frac   fraction of carriers in the largest of those groups
  carriers / pos_carriers  carriers, and how many are trait-positive
  top_study, top_platform  where the carriers come from (share)
  locus        K-10 CDS hit (locus_tag, gene, product, codon position) or
               the flanking genes if intergenic

Outputs: 06_pangwas/pyseer/snp_hits_annotated.tsv and a printed summary.
"""

import argparse
import importlib.util
import re
from bisect import bisect_left
from pathlib import Path

import pandas as pd

ROOT = Path("/cfs/earth/scratch/tabasran/MAP")


def load27(path):
    sp = importlib.util.spec_from_file_location("p27", path)
    m = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(m)
    return m


def rooted_masks(adj, name, S, tips):
    """Root at the S|C edge; return (leafmask, parentmask) for every node.
    Masks are Python ints with one bit per tip (order of `tips`)."""
    idx = {t: k for k, t in enumerate(tips)}
    allm = (1 << len(tips)) - 1
    S, edge = set(S), None
    for a in adj:
        for b, _ in adj[a]:
            if a < b:
                st, seen, side = [b], {a, b}, set()
                while st:
                    x = st.pop()
                    if x in name:
                        side.add(name[x])
                    for y, _l in adj[x]:
                        if y not in seen:
                            seen.add(y)
                            st.append(y)
                if side == S or set(tips) - side == S:
                    edge = (a, b)
                    break
        if edge:
            break
    if edge is None:
        raise SystemExit("[!] no branch separates S from C")
    out = []
    for start, other in (edge, edge[::-1]):
        order, parent = [], {start: other}
        st = [start]
        while st:
            x = st.pop()
            order.append(x)
            for y, _l in adj[x]:
                if y != parent[x]:
                    parent[y] = x
                    st.append(y)
        mask = {}
        for x in reversed(order):
            kids = [y for y, _l in adj[x] if y != parent[x]]
            mask[x] = 0 if kids else (1 << idx[name[x]])
            for y in kids:
                mask[x] |= mask[y]
        for x in order:
            out.append((mask[x], allm if x == start else mask[parent[x]]))
    return out


def supported_splits(text, tips, thr):
    """Leaf masks below every internal branch whose UFBoot >= thr (original
    rooting of the Newick; a split is the same whichever side is taken)."""
    idx = {t: k for k, t in enumerate(tips)}
    text = text.strip().rstrip(";")
    stack, masks, i, n = [[]], [], 0, len(text)
    while i < n:
        c = text[i]
        if c == "(":
            stack.append([])
            i += 1
        elif c == ",":
            i += 1
        elif c == ")":
            m = 0
            for x in stack.pop():
                m |= x
            i += 1
            j = i
            while j < n and text[j] not in ",():;":
                j += 1
            lab = text[i:j].strip()
            try:
                sup = float(lab.split("/")[-1])
            except ValueError:
                sup = None
            if sup is not None and sup >= thr:
                masks.append(m)
            stack[-1].append(m)
            while j < n and text[j] not in ",()":
                j += 1
            i = j
        else:
            j = i
            while j < n and text[j] not in ",():":
                j += 1
            stack[-1].append(1 << idx[text[i:j].strip().strip("'\"")])
            while j < n and text[j] not in ",()":
                j += 1
            i = j
    return masks


def conflicts(X, U, splits):
    """Number of supported splits incompatible with carrier set X within U."""
    k = 0
    for m in splits:
        A = m & U
        B = U & ~m
        if not A or not B:
            continue
        if (X & A) and (X & B) and (A & ~X) and (B & ~X):
            k += 1
    return k


def popcount(x):
    return bin(x).count("1")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=str(ROOT / "06_pangwas/pyseer"))
    ap.add_argument("--tree", default=str(ROOT / "05_phylogeny/iqtree/tree_final.nwk"))
    ap.add_argument("--lineages", default=str(ROOT / "05_phylogeny/typing/lineage_assignments.tsv"))
    ap.add_argument("--ena", default=str(ROOT / "01_metadata/ena_raw.tsv"))
    ap.add_argument("--gff", default=str(ROOT / "03_genomes/gff/GCA_000007865.1.gff"))
    ap.add_argument("--script27", default=str(ROOT / "scripts/27_pyseer_prep.py"))
    ap.add_argument("--support", type=float, default=95.0)
    args = ap.parse_args()
    d = Path(args.dir)

    sig = pd.read_csv(d / "significant.tsv", sep="\t")
    sig = sig[(sig["variants_type"] == "snps") & (sig["trait"] != "neg_lowN50")].copy()
    V = pd.read_csv(d / "variants_snps.Rtab", sep="\t", index_col=0)
    tips = list(V.columns)
    lin = pd.read_csv(args.lineages, sep="\t", dtype=str).set_index("genome").reindex(tips)
    S = list(lin.index[lin["lineage"] == "S"])
    m27 = load27(args.script27)
    adj, name = m27.parse(open(args.tree).read())
    nodes = rooted_masks(adj, name, S, tips)
    idx = {t: k for k, t in enumerate(tips)}
    splits = supported_splits(open(args.tree).read(), tips, args.support)
    print("internal branches with UFBoot >= {:g}: {}".format(args.support, len(splits)))

    ena = pd.read_csv(args.ena, sep="\t", dtype=str).drop_duplicates("run_accession").set_index("run_accession")
    study = pd.Series([ena.at[g, "study_accession"] if g in ena.index else "assembly" for g in tips], index=tips)
    plat = pd.Series([ena.at[g, "instrument_model"] if g in ena.index else "assembly" for g in tips], index=tips)

    # K-10 CDS for annotation
    cds = []
    for line in open(args.gff):
        if line.startswith("##FASTA"):
            break
        f = line.rstrip("\n").split("\t")
        if line.startswith("#") or len(f) < 9 or f[2] != "CDS":
            continue
        at = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        cds.append((int(f[3]), int(f[4]), f[6], at.get("locus_tag", ""), at.get("gene", ""),
                    at.get("product", "")))
    cds.sort()
    starts = [c[0] for c in cds]

    def annotate(pos):
        i = bisect_left(starts, pos + 1) - 1
        for j in (i, i - 1):
            if 0 <= j < len(cds) and cds[j][0] <= pos <= cds[j][1]:
                s, e, strand, lt, g, p = cds[j]
                off = (pos - s) if strand == "+" else (e - pos)
                return "CDS", "{} {} | {} | codon pos {}".format(lt, g, p, off % 3 + 1)
        left = cds[i] if i >= 0 else None
        right = cds[i + 1] if i + 1 < len(cds) else None
        return "intergenic", "between {} {} and {} {}".format(
            left[3] if left else "-", left[4] if left else "", right[3] if right else "-",
            right[4] if right else "")

    rows = []
    for _, h in sig.iterrows():
        pheno = pd.read_csv(d / "pheno" / "{}__{}.tsv".format(h["subset"], h["trait"]), sep="\t")
        samp = list(pheno["sample"])
        y = pheno.set_index("sample").iloc[:, 0]
        x = V.loc[h["variant"], samp]
        if x.mean() > 0.5:                   # use the minor allele as 'carrier'
            x = 1 - x
            allele = "ref"
        else:
            allele = "alt"
        car = list(x.index[x == 1])
        X = 0
        for g in car:
            X |= 1 << idx[g]
        sub = 0
        for g in samp:
            sub |= 1 << idx[g]
        # maximal clades: node mask (restricted to subset) inside X, parent's not
        groups = []
        for nm, pm in nodes:
            nm_s, pm_s = nm & sub, pm & sub
            if nm_s and (nm_s & ~X) == 0 and (pm_s & ~X) != 0:
                groups.append(popcount(nm_s))
        nconf = conflicts(X, sub, splits)
        pos = int(h["variant"].split("_")[-3])
        kind, locus = annotate(pos)
        st = study[car].value_counts()
        pl = plat[car].value_counts()
        rows.append({"trait": h["trait"], "subset": h["subset"], "variant": h["variant"],
                     "p": h["lrt-pvalue"], "beta": round(h["beta"], 3), "minor": allele,
                     "carriers": len(car), "pos_carriers": int(y[car].sum()),
                     "pos_rate_carriers": round(float(y[car].mean()), 3),
                     "pos_rate_all": round(float(y.mean()), 3),
                     "origins": len(groups),
                     "conflicts_supported": nconf,
                     "verdict": "clade-compatible" if nconf == 0 else "homoplastic",
                     "largest_clade_frac": round(max(groups) / len(car), 3) if groups else 0,
                     "top_study": "{} ({:.0%})".format(st.index[0], st.iloc[0] / len(car)),
                     "top_platform": "{} ({:.0%})".format(pl.index[0], pl.iloc[0] / len(car)),
                     "kind": kind, "locus": locus})
    R = pd.DataFrame(rows).sort_values(["subset", "trait", "p"])
    R.to_csv(d / "snp_hits_annotated.tsv", sep="\t", index=False)

    print("significant SNP hits (excluding the negative control): {}".format(len(R)))
    print("\nper run: hits | 1 origin on ML tree | compatible with all UFBoot>={:g} clades | "
          "truly homoplastic | median carriers".format(args.support))
    for (s, t), g in R.groupby(["subset", "trait"]):
        print("  {:<7} {:<20} {:>4} | {:>4} | {:>4} | {:>4} | {:>6.0f}".format(
            s, t, len(g), int((g["origins"] == 1).sum()),
            int((g["verdict"] == "clade-compatible").sum()),
            int((g["verdict"] == "homoplastic").sum()), g["carriers"].median()))
    pd.set_option("display.width", 260)
    pd.set_option("display.max_colwidth", 70)
    print("\n--- C_nonB hits incompatible with well-supported clades (candidate convergent) ---")
    hom = R[(R["subset"] == "C_nonB") & (R["verdict"] == "homoplastic")]
    print(hom[["trait", "variant", "p", "carriers", "pos_rate_carriers", "pos_rate_all",
               "origins", "conflicts_supported", "largest_clade_frac", "top_study", "kind",
               "locus"]].to_string(index=False)
          if len(hom) else "  none")
    print("\n--- C_nonB clade-compatible hits (one row per carrier group) ---")
    cm = R[(R["subset"] == "C_nonB") & (R["verdict"] == "clade-compatible")].drop_duplicates(["trait", "carriers", "top_study"])
    print(cm[["trait", "variant", "carriers", "pos_rate_carriers", "pos_rate_all", "top_study",
              "top_platform", "kind", "locus"]].head(30).to_string(index=False))
    print("\nwritten: {}/snp_hits_annotated.tsv".format(d))


if __name__ == "__main__":
    main()
