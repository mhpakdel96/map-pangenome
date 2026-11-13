#!/usr/bin/env python3
"""
27_pyseer_prep.py

Inputs for the association analysis with pyseer (linear mixed model).
Replaces the Scoary analysis, whose hits were shown to be assembly/annotation
artefacts (scripts 21-23). Differences from Scoary:
  * variants = (a) accessory genes from the DNA-validated union80 matrix
    (04_pangenome/dna_presence/presence_final.Rtab, clusters >= 250 bp) and
    (b) core-genome SNPs from snippy-core (core.tab, mapping-based, so not
    dependent on gene calling)
  * population structure: LMM random effect with a phylogenetic similarity
    matrix (shared branch length from the root to the MRCA of each pair, as
    pyseer's phylogeny_distance.py --lmm), computed here from the IQ-TREE ML
    tree rooted in the middle of the branch that separates S from C
  * covariates: log10 N50 (quantitative) and sequencing platform
    (categorical; platforms with < --min-platform genomes pooled as 'other')
    so that residual assembly-quality/platform effects are absorbed
  * built-in NEGATIVE control phenotype 'neg_lowN50': within every study with
    >= --min-study genomes, 1 = below that study's median N50, 0 = above.
    Run WITHOUT covariates; with a clean variant matrix it should give ~no
    hits. (The same contrast gave 115 Fisher hits on union80 in script 26.)

Subsets: C_nonB (primary, as before) and all.

Outputs (06_pangwas/pyseer/)
  kinship_<subset>.tsv    similarity matrix for --similarity
  covariates.tsv          sample, log10_n50, platform   (use: 2q 3)
  variants_genes.Rtab     accessory genes (1-99% in the cohort)
  variants_snps.Rtab      one row per core SNP allele (POS_REF_ALT)
  pheno/<subset>__<trait>.tsv  one phenotype, samples with NA dropped
  manifest.tsv            one line per pyseer run (array input)
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/cfs/earth/scratch/tabasran/MAP")


# ---------------------------------------------------------------- tree
def parse(text):
    """Return adjacency {node: [(nbr, length)]}, tip names {node: name}."""
    text = text.strip().rstrip(";")
    adj, name = {}, {}
    nid = [0]

    def new():
        nid[0] += 1
        adj[nid[0]] = []
        return nid[0]
    root = new()
    stack, cur, last = [], root, None
    i, n = 0, len(text)

    def label(i):
        j = i
        while j < n and text[j] not in ",():;":
            j += 1
        return text[i:j].strip(), j
    while i < n:
        c = text[i]
        if c == "(":
            ch = new()
            adj[cur].append([ch, 0.0])
            adj[ch].append([cur, 0.0])
            stack.append(cur)
            cur = ch
            i += 1
        elif c == ",":
            par = stack[-1]
            ch = new()
            adj[par].append([ch, 0.0])
            adj[ch].append([par, 0.0])
            cur = ch
            i += 1
        elif c == ")":
            cur = stack.pop()
            i += 1
            _, i = label(i)                      # support label, ignored
        elif c == ":":
            v, i = label(i + 1)
            L = float(v) if v else 0.0
            par = stack[-1] if stack else None
            if par is not None:
                for e in adj[cur]:
                    if e[0] == par:
                        e[1] = L
                for e in adj[par]:
                    if e[0] == cur:
                        e[1] = L
        else:
            lab, i = label(i)
            name[cur] = lab.strip("'\"")
    return adj, name


def leaves_below(adj, name, a, b):
    """Tip names on b's side of edge a-b."""
    out, st, seen = set(), [b], {a, b}
    while st:
        x = st.pop()
        if x in name:
            out.add(name[x])
        for y, _ in adj[x]:
            if y not in seen:
                seen.add(y)
                st.append(y)
    return out


def similarity(adj, name, S):
    """Root at the midpoint of the edge separating set S from the rest;
    C[i,j] = distance from root to MRCA(i,j); C[i,i] = root-to-tip."""
    S = set(S)
    edge = None
    for a in adj:
        for b, L in adj[a]:
            if a < b:
                side = leaves_below(adj, name, a, b)
                if side == S or (set(name.values()) - side) == S:
                    edge = (a, b, L)
                    break
        if edge:
            break
    if edge is None:
        sys.exit("[!] no branch separates the S lineage from the rest")
    a, b, L = edge
    tips = sorted(name.values())
    idx = {t: k for k, t in enumerate(tips)}
    C = np.zeros((len(tips), len(tips)))
    # iterative DFS from each side of the root edge, depth starts at L/2
    for start, other in ((a, b), (b, a)):
        order, parent, depth = [], {start: other}, {start: L / 2.0}
        st = [start]
        while st:
            x = st.pop()
            order.append(x)
            for y, l in adj[x]:
                if y != parent[x]:
                    parent[y] = x
                    depth[y] = depth[x] + l
                    st.append(y)
        below = {}
        for x in reversed(order):                 # post-order
            kids = [y for y, _ in adj[x] if y != parent[x]]
            if not kids:
                below[x] = [idx[name[x]]]
                C[below[x][0], below[x][0]] = depth[x]
                continue
            groups = [below[y] for y in kids]
            for gi in range(len(groups)):
                for gj in range(gi + 1, len(groups)):
                    C[np.ix_(groups[gi], groups[gj])] = depth[x]
                    C[np.ix_(groups[gj], groups[gi])] = depth[x]
            below[x] = [k for g in groups for k in g]
            for y in kids:
                del below[y]
    # pairs across the root share depth 0 (already 0)
    return pd.DataFrame(C, index=tips, columns=tips), L


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree", default=str(ROOT / "05_phylogeny/iqtree/tree_final.nwk"))
    ap.add_argument("--genes", default=str(ROOT / "04_pangenome/dna_presence/presence_final.Rtab"))
    ap.add_argument("--coretab", default=str(ROOT / "05_phylogeny/core_masked/core.tab"))
    ap.add_argument("--traits", default=str(ROOT / "06_pangwas/trait.csv"))
    ap.add_argument("--lineages", default=str(ROOT / "05_phylogeny/typing/lineage_assignments.tsv"))
    ap.add_argument("--qc", default=str(ROOT / "02_qc/genome_qc_master.tsv"))
    ap.add_argument("--ena", default=str(ROOT / "01_metadata/ena_raw.tsv"))
    ap.add_argument("--out", default=str(ROOT / "06_pangwas/pyseer"))
    ap.add_argument("--min-n", type=int, default=20)
    ap.add_argument("--min-platform", type=int, default=20)
    ap.add_argument("--min-study", type=int, default=60)
    args = ap.parse_args()
    out = Path(args.out)
    if (out / "manifest.tsv").exists():
        sys.exit("[!] {} already prepared — move it away first".format(out))
    (out / "pheno").mkdir(parents=True, exist_ok=True)

    tr = pd.read_csv(args.traits, dtype=str, keep_default_na=False).set_index("Name")
    iso = list(tr.index)
    lin = pd.read_csv(args.lineages, sep="\t", dtype=str).set_index("genome").reindex(iso)
    lab = lin["lineage"] + "/" + lin["sublineage"]
    S = set(lab.index[lin["lineage"] == "S"])

    # -------------------------------------------------- kinship
    adj, name = parse(open(args.tree).read())
    if set(name.values()) != set(iso):
        sys.exit("[!] tree tips != cohort")
    K, L = similarity(adj, name, S)
    print("kinship: rooted on the S|C branch (length {:.6g}); diag median {:.4g}".format(
        L, float(np.median(np.diag(K.values)))))
    subsets = {"all": iso, "C_nonB": [g for g in iso if lab[g] == "C/non-B"]}
    for s, members in subsets.items():
        K.loc[members, members].to_csv(out / "kinship_{}.tsv".format(s), sep="\t",
                                       float_format="%.6g")
        print("  kinship_{}.tsv  {} x {}".format(s, len(members), len(members)))

    # -------------------------------------------------- covariates
    qc = pd.read_csv(args.qc, sep="\t").set_index("genome").reindex(iso)
    ena = pd.read_csv(args.ena, sep="\t", dtype=str).drop_duplicates("run_accession").set_index("run_accession")
    plat = pd.Series([ena.at[g, "instrument_model"] if g in ena.index else "assembly" for g in iso], index=iso)
    vc = plat.value_counts()
    plat = plat.where(plat.map(vc) >= args.min_platform, "other").str.replace(" ", "_")
    cov = pd.DataFrame({"sample": iso, "log10_n50": np.log10(qc["n50"].values).round(4),
                        "platform": plat.values})
    cov.to_csv(out / "covariates.tsv", sep="\t", index=False)
    print("covariates: platform levels {}".format(dict(plat.value_counts())))

    # -------------------------------------------------- variants
    G = pd.read_csv(args.genes, sep="\t", index_col=0)
    G = G[iso]
    f = G.mean(axis=1)
    Ga = G[(f > 0.01) & (f < 0.99)]
    Ga.to_csv(out / "variants_genes.Rtab", sep="\t", index_label="Gene")
    print("genes: {} accessory (1-99%) of {}".format(len(Ga), len(G)))

    ct = pd.read_csv(args.coretab, sep="\t", dtype=str)
    # snippy-core versions name the reference-allele column differently
    refcol = next((c for c in ("Reference", "REF", "Ref", "reference", "ref") if c in ct.columns), None)
    if refcol is None or not {"CHR", "POS"} <= set(ct.columns):
        sys.exit("[!] core.tab header not recognised; first columns: {}".format(list(ct.columns[:6])))
    print("core.tab: reference column '{}'".format(refcol))
    miss = [g for g in iso if g not in ct.columns]
    if miss:
        sys.exit("[!] core.tab lacks samples: {}".format(miss[:5]))
    rows, names = [], []
    A = ct[iso].apply(lambda c: c.str.upper())
    for k in range(len(ct)):
        ref = ct.at[k, refcol].upper()
        al = A.iloc[k].values
        for alt in sorted(set(al) - {ref, "N", "-"}):
            names.append("{}_{}_{}_{}".format(ct.at[k, "CHR"], ct.at[k, "POS"], ref, alt))
            rows.append((al == alt).astype(np.int8))
    V = pd.DataFrame(np.vstack(rows), index=names, columns=iso)
    V.to_csv(out / "variants_snps.Rtab", sep="\t", index_label="Variant")
    print("snps: {} sites -> {} allele rows".format(len(ct), len(V)))

    # -------------------------------------------------- phenotypes
    man = []

    def write(subset, trait, series, covs):
        known = series.dropna()
        pos, neg = int((known == 1).sum()), int((known == 0).sum())
        if min(pos, neg) < args.min_n:
            print("  skip {}__{}: pos {} neg {}".format(subset, trait, pos, neg))
            return
        p = out / "pheno" / "{}__{}.tsv".format(subset, trait)
        pd.DataFrame({"sample": known.index, trait: known.astype(int).values}).to_csv(p, sep="\t", index=False)
        for v in ("genes", "snps"):
            man.append((subset, trait, v, pos, neg, covs))

    for s, members in subsets.items():
        for t in tr.columns:
            x = pd.to_numeric(tr.loc[members, t].replace("NA", np.nan), errors="coerce")
            write(s, t, x, "2q 3")
    # negative control
    study = pd.Series([ena.at[g, "study_accession"] if g in ena.index else "assembly" for g in iso], index=iso)
    for s, members in subsets.items():
        neg = pd.Series(np.nan, index=members)
        for st, n in study[members].value_counts().items():
            if st == "assembly" or n < args.min_study:
                continue
            m = [g for g in members if study[g] == st]
            med = qc.loc[m, "n50"].median()
            neg[m] = (qc.loc[m, "n50"] < med).astype(int).values
        write(s, "neg_lowN50", neg, "none")

    with open(out / "manifest.tsv", "w") as fh:
        fh.write("run\tsubset\ttrait\tvariants\tpos\tneg\tcovariates\n")
        for subset, trait, v, pos, neg, covs in man:
            fh.write("{0}__{1}__{2}\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(subset, trait, v, pos, neg, covs))
    print("runs: {} -> {}/manifest.tsv".format(len(man), out))


if __name__ == "__main__":
    main()
