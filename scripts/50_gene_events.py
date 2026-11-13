#!/usr/bin/env python3
"""
50_gene_events.py

Details of the convergent genes found by 49 (default: BH q < 0.05 in
06_pangwas/convergence_v2/gene_branch_convergence.tsv), above all icl
(isocitrate lyase: 15 independent branch events, all nonsynonymous; site
4421235 C>T on 10 branches).

For every read-verified SNP allele in the selected genes:
  amino-acid change in K-10 coordinates (e.g. A123V), direction of each
  change (gain of ALT or reversion to REF), and for EVERY branch with a
  change: number of tips below, and lineage / region / host / study of the
  genomes that carry the changed state below that branch.
Per gene: distinct amino-acid positions hit, amino-acid changes that arose
more than once (same residue, same change = strongest signal), and other
K-10 CDS with the same product (possible paralogs -> mapping risk).

  python 50_gene_events.py [--genes icl,pafA,BDOFFIBO_04561]
Output: 06_pangwas/convergence_v2/gene_events_detail.tsv and printed report
Environment: ml_env. Imports 44_homoplasy_site.py.
"""

import argparse
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/cfs/earth/scratch/tabasran/MAP")


def load(path, nm):
    sp = importlib.util.spec_from_file_location(nm, path)
    m = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(m)
    return m


def reconstruct(children, post, root, tipidx, M):
    """Fitch down-pass + top-down assignment; returns states {node: int8 array}."""
    s0, s1 = {}, {}
    for x in post:
        kids = children.get(x, [])
        if not kids:
            r = M[tipidx[x]]
            s0[x] = (r == 0) | (r == -1)
            s1[x] = (r == 1) | (r == -1)
            continue
        a0, a1 = s0[kids[0]], s1[kids[0]]
        for k in kids[1:]:
            b0, b1 = s0[k], s1[k]
            i0, i1 = a0 & b0, a1 & b1
            e = ~(i0 | i1)
            a0 = np.where(e, a0 | b0, i0)
            a1 = np.where(e, a1 | b1, i1)
        s0[x], s1[x] = a0, a1
    state = {root: np.where(s0[root], 0, 1).astype(np.int8)}
    stack = [root]
    while stack:
        p = stack.pop()
        for c in children.get(p, []):
            st = state[p]
            ok = np.where(st == 0, s0[c], s1[c])
            state[c] = np.where(ok, st, 1 - st).astype(np.int8)
            stack.append(c)
    return state


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes", default="", help="locus tags or gene names, comma-separated")
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--min-dp", type=int, default=10)
    ap.add_argument("--min-af", type=float, default=0.8)
    ap.add_argument("--conv", default=str(ROOT / "06_pangwas/convergence_v2/gene_branch_convergence.tsv"))
    ap.add_argument("--tree", default=str(ROOT / "05_phylogeny/iqtree/tree_final.nwk"))
    ap.add_argument("--coretab", default=str(ROOT / "05_phylogeny/core_masked/core.tab"))
    ap.add_argument("--reads", default=str(ROOT / "08_reads/k10"))
    ap.add_argument("--ref", default=str(ROOT / "03_genomes/fasta/GCA_000007865.1.fasta"))
    ap.add_argument("--gff", default=str(ROOT / "03_genomes/gff/GCA_000007865.1.gff"))
    ap.add_argument("--lineages", default=str(ROOT / "05_phylogeny/typing/lineage_assignments.tsv"))
    ap.add_argument("--traits", default=str(ROOT / "06_pangwas/trait.csv"))
    ap.add_argument("--ena", default=str(ROOT / "01_metadata/ena_raw.tsv"))
    ap.add_argument("--s44", default=str(Path(__file__).resolve().parent / "44_homoplasy_site.py"))
    ap.add_argument("--out", default=str(ROOT / "06_pangwas/convergence_v2"))
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    s44 = load(args.s44, "s44")

    # ------------------------------------------------ genes
    cds = []
    for line in open(args.gff):
        if line.startswith("##FASTA"):
            break
        f = line.rstrip("\n").split("\t")
        if line.startswith("#") or len(f) < 9 or f[2] != "CDS":
            continue
        at = dict(x.split("=", 1) for x in f[8].split(";") if "=" in x)
        cds.append((int(f[3]), int(f[4]), f[6], at.get("locus_tag", "?"), at.get("gene", ""),
                    at.get("product", "").replace("%2C", ",")))
    C = pd.DataFrame(sorted(cds), columns=["start", "end", "strand", "locus_tag", "gene", "product"])
    if args.genes:
        want = [x.strip() for x in args.genes.split(",") if x.strip()]
        sel = C[C["locus_tag"].isin(want) | C["gene"].isin(want)]
    else:
        G = pd.read_csv(args.conv, sep="\t")
        sel = C[C["locus_tag"].isin(G.loc[G["q"] < args.fdr, "locus_tag"])]
    if sel.empty:
        raise SystemExit("[!] no gene selected")
    print("genes: {}".format(", ".join("{} {}".format(a, b or "-") for a, b in zip(sel["locus_tag"], sel["gene"]))))

    seq, cur = {}, None
    for line in open(args.ref):
        if line.startswith(">"):
            cur = line[1:].split()[0]
            seq[cur] = []
        elif cur:
            seq[cur].append(line.strip().upper())
    g = "".join(next(iter(seq.values())))

    def aa_change(r, pos, alt):
        s, e, strand = r["start"], r["end"], r["strand"]
        if strand == "+":
            off = pos - s
            cs = s + (off // 3) * 3
            codon = g[cs - 1:cs + 2]
            cp = off % 3
            ac = codon[:cp] + alt + codon[cp + 1:]
        else:
            off = e - pos
            ce = e - (off // 3) * 3
            codon = s44.rc(g[ce - 3:ce])
            cp = off % 3
            ac = codon[:cp] + alt.translate(s44.COMP) + codon[cp + 1:]
        a1, a2 = s44.CODE.get(codon, "?"), s44.CODE.get(ac, "?")
        return "{}{}{}".format(a1, off // 3 + 1, a2), ("S" if a1 == a2 else "N"), "{}>{}".format(codon, ac)

    # ------------------------------------------------ tree + metadata
    adj, name, sup = s44.parse_newick(open(args.tree).read())
    tips = sorted(name.values())
    ti = {t: k for k, t in enumerate(tips)}
    lin = pd.read_csv(args.lineages, sep="\t", dtype=str).set_index("genome")
    lab = (lin["lineage"] + "/" + lin["sublineage"]).reindex(tips).fillna("?")
    children, parent, post, root, edge = s44.root_on_split(adj, name, [t for t in tips if lab[t].startswith("S/")])
    tipidx = {nd: ti[nm] for nd, nm in name.items()}
    below = {}
    for nd in post:
        kids = children.get(nd, [])
        below[nd] = [tipidx[nd]] if not kids else [i for k in kids for i in below[k]]
    tr = pd.read_csv(args.traits, dtype=str).set_index("Name")
    reg = pd.Series("unknown", index=tips)
    for c in ("Region_Americas", "Region_Europe", "Region_Oceania"):
        if c in tr.columns:
            reg[[t for t in tips if tr[c].get(t) == "1"]] = c.replace("Region_", "")
    host = pd.Series("NA", index=tips)
    if "hostCattle" in tr.columns:
        host[[t for t in tips if tr["hostCattle"].get(t) == "1"]] = "cattle"
        host[[t for t in tips if tr["hostCattle"].get(t) == "0"]] = "non-cattle"
    ena = pd.read_csv(args.ena, sep="\t", dtype=str).drop_duplicates("run_accession").set_index("run_accession")
    study = pd.Series([ena.at[t, "study_accession"] if t in ena.index else "assembly" for t in tips], index=tips)

    # ------------------------------------------------ genotypes at selected sites, read-verified
    ct = pd.read_csv(args.coretab, sep="\t", dtype=str)
    refcol = "REF" if "REF" in ct.columns else "Reference"
    ct["POS"] = ct["POS"].astype(int)
    inside = np.zeros(len(ct), dtype=bool)
    for _, r in sel.iterrows():
        inside |= (ct["POS"] >= r["start"]).values & (ct["POS"] <= r["end"]).values
    sub = ct[inside].reset_index(drop=True)
    P = sub["POS"].values
    BA = sub[tips].apply(lambda c: c.str.upper()).fillna("N").values
    want = set(P.tolist())
    for t in tips:
        f = Path(args.reads) / (t + ".ad.tsv")
        if not f.exists():
            continue
        rc = {}
        for line in open(f):
            x = line.rstrip("\n").split("\t")
            if len(x) < 4 or len(x[1]) != 1:
                continue
            p = int(x[0])
            if p not in want or p in rc:
                continue
            al = [x[1].upper()] + [a.upper() for a in x[2].split(",")]
            ad = [int(v) if v not in (".", "") else 0 for v in x[3].split(",")]
            tot = sum(ad)
            b = int(np.argmax(ad)) if ad else 0
            rc[p] = al[b] if tot >= args.min_dp and b < len(al) and len(al[b]) == 1 and ad[b] / tot >= args.min_af else None
        j = ti[t]
        for k, p in enumerate(P):
            if BA[k, j] not in ("N", "-", "?") and rc.get(p) != BA[k, j]:
                BA[k, j] = "N"
    meta, cols = [], []
    for k in range(len(P)):
        row = BA[k]
        miss = np.isin(row, ["N", "-", "?"])
        for alt in sorted(set(row[~miss]) - {sub.at[k, refcol].upper()}):
            m = (row == alt).astype(np.int8)
            m[miss] = -1
            meta.append((int(P[k]), sub.at[k, refcol].upper(), alt))
            cols.append(m)
    if not cols:
        raise SystemExit("[!] no verified SNP in the selected genes")
    M = np.vstack(cols).T
    state = reconstruct(children, post, root, tipidx, M)

    def describe(idx):
        gl = [tips[i] for i in idx]
        f = lambda s: "; ".join("{}:{}".format(a, b) for a, b in s[gl].value_counts().head(3).items())
        return f(lab), f(reg), f(host), f(study)

    rows = []
    for _, r in sel.iterrows():
        print("\n" + "=" * 110)
        print("{} {} | {} | {}-{} ({}) {} bp".format(r["locus_tag"], r["gene"] or "-", r["product"],
                                                     r["start"], r["end"], r["strand"], r["end"] - r["start"] + 1))
        para = C[(C["product"] == r["product"]) & (C["locus_tag"] != r["locus_tag"])]
        if len(para) and r["product"] not in ("hypothetical protein", ""):
            print("  other K-10 CDS with the same product: " + "; ".join(
                "{} {} {}-{}".format(a, b or "-", s, e) for a, b, s, e in
                zip(para["locus_tag"], para["gene"], para["start"], para["end"])))
        aa_seen = {}
        for j, (pos, ref, alt) in enumerate(meta):
            if not (r["start"] <= pos <= r["end"]):
                continue
            aa, eff, cod = aa_change(r, pos, alt)
            car = int((M[:, j] == 1).sum())
            evs = []
            for nd in post:
                if nd == root:
                    continue
                p = parent[nd]
                if state[nd][j] != state[p][j]:
                    newst = int(state[nd][j])
                    idx = [i for i in below[nd] if M[i, j] == newst]
                    lg, rg, hs, sd = describe(idx)
                    evs.append({"tips_below": len(below[nd]), "carriers_below": len(idx),
                                "direction": "gain " + alt if newst == 1 else "reversion to " + ref,
                                "terminal": nd in name,
                                "genome": name.get(nd, ""), "lineages": lg, "regions": rg, "hosts": hs, "studies": sd})
            if not evs:
                continue
            aa_seen.setdefault(aa, 0)
            aa_seen[aa] += len([e for e in evs if e["direction"].startswith("gain")])
            print("  {} {}>{}  {} {} ({})  carriers {}  events {}".format(pos, ref, alt, aa, eff, cod, car, len(evs)))
            for e in sorted(evs, key=lambda z: -z["carriers_below"]):
                print("      - {:<22} {:>4} genomes{} | {} | {} | {} | {}".format(
                    e["direction"], e["carriers_below"],
                    " ({})".format(e["genome"]) if e["terminal"] else "",
                    e["lineages"], e["regions"], e["hosts"], e["studies"]))
                rows.append(dict(locus_tag=r["locus_tag"], gene=r["gene"], product=r["product"], pos=pos,
                                 ref=ref, alt=alt, aa=aa, effect=eff, codon=cod, carriers=car, **e))
        rec = {k: v for k, v in aa_seen.items() if v >= 2}
        print("  distinct amino-acid changes: {} | arose >= 2 times: {}".format(
            len(aa_seen), ", ".join("{} x{}".format(k, v) for k, v in rec.items()) or "none"))
    pd.DataFrame(rows).to_csv(out / "gene_events_detail.tsv", sep="\t", index=False)
    print("\nwritten: {}/gene_events_detail.tsv".format(out))


if __name__ == "__main__":
    main()
