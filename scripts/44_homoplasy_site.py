#!/usr/bin/env python3
"""
44_homoplasy_site.py

Why is one core SNP homoplastic? Default: AE016958.1_374739_T_C, the only
pattern whose carriers were REJECTED as one clade by the AU test (37:
dlogL 615, p-AU 0.0033, also rejected by KH and SH). Alleles are clean in
reads (35, PART 4) and Gubbins found no block at the site (31). Remaining
explanations and what each predicts:

  R  undetected RECOMBINATION (import of a fragment carrying the allele):
     neighbouring SNPs share the same incongruent carrier set (linkage),
     carriers have more SNPs in the window, a Gubbins block nearby
  M  recurrent MUTATION / convergence: the site is isolated (neighbours follow
     the tree), carriers fall into several independent groups
  P  PARALOG / mapping artefact: site near a masked repeat (IS, PE/PPE), read
     depth at the site deviates in carriers (copy-number difference)
  T  the ML tree is locally wrong: many UNLINKED sites across the genome share
     the same carrier set

Analyses
  1. locus: K-10 CDS, codon position, amino-acid change; nearest masked
     repeat (05_phylogeny/core_masked/mask_genes.tsv)
  2. allele counts per lineage (which allele is ancestral, given S|C root)
  3. Fitch parsimony on the ML tree for EVERY core SNP allele (binary,
     N/- = missing) -> genome-wide distribution and the rank of this site;
     saved once as 05_phylogeny/homoplasy/fitch_all_sites.tsv
  4. independent carrier groups: maximal clades whose tips all carry the
     minor allele, on the tree rooted on the S|C branch, with UFBoot support,
     region, study and lineage of each group
  5. neighbours within +-window: carrier overlap (Jaccard) with this site and
     their own Fitch score; SNP density in the window, carriers vs others
  6. genome-wide sites with a similar carrier set (Jaccard >= --jaccard)
  7. Gubbins blocks overlapping the window
  8. reads (08_reads/k10): allele fraction and depth ratio (site depth /
     genome mean depth) at this site and at the window sites, carriers vs
     non-carriers

  python 44_homoplasy_site.py [--variant AE016958.1_374739_T_C] [--window 10000]
Outputs: 05_phylogeny/homoplasy/<variant>/{report.txt, groups.tsv,
         neighbours.tsv, similar_patterns.tsv, reads.tsv}
Environment: ml_env (numpy, pandas, scipy).
"""

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/cfs/earth/scratch/tabasran/MAP")
CODE = {}
_b = "TCAG"
_aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for _i, _x in enumerate(_b):
    for _j, _y in enumerate(_b):
        for _k, _z in enumerate(_b):
            CODE[_x + _y + _z] = _aa[16 * _i + 4 * _j + _k]
COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(COMP)[::-1]


# ---------------------------------------------------------------- tree
def parse_newick(text):
    """Unrooted adjacency {node: {nbr: length}}, tip names {node: name},
    support {frozenset(edge): label} (IQ-TREE writes UFBoot as the label
    of the node below an edge)."""
    text = text.strip().rstrip(";")
    adj, name, sup = {}, {}, {}
    nid = [0]

    def new():
        nid[0] += 1
        adj[nid[0]] = {}
        return nid[0]

    root = new()
    stack, cur = [], root
    i, n = 0, len(text)

    def token(i):
        if i < n and text[i] in "'\"":
            q = text[i]
            j = text.index(q, i + 1)
            return text[i + 1:j], j + 1
        j = i
        while j < n and text[j] not in ",():;":
            j += 1
        return text[i:j].strip(), j

    def link(a, b):
        adj[a][b] = 0.0
        adj[b][a] = 0.0

    while i < n:
        c = text[i]
        if c == "(":
            ch = new()
            link(cur, ch)
            stack.append(cur)
            cur = ch
            i += 1
        elif c == ",":
            ch = new()
            link(stack[-1], ch)
            cur = ch
            i += 1
        elif c == ")":
            # after the pop, `cur` is the internal node just closed; its
            # parent (if any) is on top of the stack
            cur = stack.pop()
            i += 1
            lab, i = token(i)
            if lab and stack:
                sup[frozenset((cur, stack[-1]))] = lab
        elif c == ":":
            v, i = token(i + 1)
            L = float(v) if v else 0.0
            par = stack[-1] if stack else None
            if par is not None:
                adj[cur][par] = L
                adj[par][cur] = L
        elif c.isspace():
            i += 1
        else:
            lab, i = token(i)
            name[cur] = lab
    return adj, name, sup


def root_on_split(adj, name, S):
    """Root in the middle of the edge separating tip set S from the rest.
    Returns children {node: [kids]}, parent {node: parent}, postorder list,
    root id, and the edge (a, b)."""
    S = set(S)
    alltips = set(name.values())

    def side(a, b):
        out, st, seen = set(), [b], {a, b}
        while st:
            x = st.pop()
            if x in name:
                out.add(name[x])
            for y in adj[x]:
                if y not in seen:
                    seen.add(y)
                    st.append(y)
        return out

    edge = None
    for a in adj:
        for b in adj[a]:
            if a < b:
                s = side(a, b)
                if s == S or alltips - s == S:
                    edge = (a, b)
                    break
        if edge:
            break
    if edge is None:
        sys.exit("[!] no branch separates lineage S from the rest")
    R = max(adj) + 1
    children, parent = {R: [edge[0], edge[1]]}, {edge[0]: R, edge[1]: R}
    order, st = [], [edge[0], edge[1]]
    while st:
        x = st.pop()
        order.append(x)
        kids = [y for y in adj[x] if y != parent[x] and not (x in edge and y in edge)]
        children[x] = kids
        for y in kids:
            parent[y] = x
            st.append(y)
    post = list(reversed(order)) + [R]
    return children, parent, post, R, edge


# ---------------------------------------------------------------- Fitch
def fitch(children, post, tipidx, M):
    """M: tips x sites int8 (1 carrier, 0 not, -1 missing). Returns minimum
    number of state changes per site. Nodes with > 2 children are folded
    pairwise (exact for binary trees; IQ-TREE trees are binary)."""
    nsite = M.shape[1]
    s0, s1 = {}, {}
    cost = np.zeros(nsite, dtype=np.int32)
    for x in post:
        kids = children.get(x, [])
        if not kids:
            row = M[tipidx[x]]
            s0[x] = (row == 0) | (row == -1)
            s1[x] = (row == 1) | (row == -1)
            continue
        a0, a1 = s0.pop(kids[0]), s1.pop(kids[0])
        for k in kids[1:]:
            b0, b1 = s0.pop(k), s1.pop(k)
            i0, i1 = a0 & b0, a1 & b1
            empty = ~(i0 | i1)
            cost += empty
            a0 = np.where(empty, a0 | b0, i0)
            a1 = np.where(empty, a1 | b1, i1)
        s0[x], s1[x] = a0, a1
    return cost


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variant", default="AE016958.1_374739_T_C")
    ap.add_argument("--window", type=int, default=10000)
    ap.add_argument("--jaccard", type=float, default=0.8)
    ap.add_argument("--tree", default=str(ROOT / "05_phylogeny/iqtree/tree_final.nwk"))
    ap.add_argument("--coretab", default=str(ROOT / "05_phylogeny/core_masked/core.tab"))
    ap.add_argument("--mask", default=str(ROOT / "05_phylogeny/core_masked/mask_genes.tsv"))
    ap.add_argument("--ref", default=str(ROOT / "03_genomes/fasta/GCA_000007865.1.fasta"))
    ap.add_argument("--gff", default=str(ROOT / "03_genomes/gff/GCA_000007865.1.gff"))
    ap.add_argument("--gubbins", default=str(ROOT / "05_phylogeny/gubbins/gubbins.recombination_predictions.gff"))
    ap.add_argument("--lineages", default=str(ROOT / "05_phylogeny/typing/lineage_assignments.tsv"))
    ap.add_argument("--traits", default=str(ROOT / "06_pangwas/trait.csv"))
    ap.add_argument("--ena", default=str(ROOT / "01_metadata/ena_raw.tsv"))
    ap.add_argument("--reads", default=str(ROOT / "08_reads/k10"))
    ap.add_argument("--no-reads", action="store_true")
    ap.add_argument("--out", default=str(ROOT / "05_phylogeny/homoplasy"))
    args = ap.parse_args()
    base = Path(args.out)
    out = base / args.variant
    out.mkdir(parents=True, exist_ok=True)
    R = []

    def say(s=""):
        print(s, flush=True)
        R.append(s)

    try:
        chrom, pos, vref, valt = args.variant.rsplit("_", 3)
        pos = int(pos)
    except ValueError:
        sys.exit("[!] variant must look like CHR_POS_REF_ALT")
    say("variant {}  (chrom {}, position {}, {} -> {}), window +-{} bp".format(
        args.variant, chrom, pos, vref, valt, args.window))

    # ------------------------------------------------ metadata
    lin = pd.read_csv(args.lineages, sep="\t", dtype=str).set_index("genome")
    lab = (lin["lineage"] + "/" + lin["sublineage"])
    tr = pd.read_csv(args.traits, dtype=str).set_index("Name")
    ena = pd.read_csv(args.ena, sep="\t", dtype=str).drop_duplicates("run_accession").set_index("run_accession")

    # ------------------------------------------------ core.tab -> allele rows
    ct = pd.read_csv(args.coretab, sep="\t", dtype=str)
    refcol = "REF" if "REF" in ct.columns else ("Reference" if "Reference" in ct.columns else None)
    if refcol is None or not {"CHR", "POS"} <= set(ct.columns):
        sys.exit("[!] core.tab needs CHR, POS and REF/Reference; has {}".format(list(ct.columns)[:6]))
    adj, name, sup = parse_newick(open(args.tree).read())
    tips = sorted(name.values())
    miss = [t for t in tips if t not in ct.columns]
    if miss:
        sys.exit("[!] {} tree tips not in core.tab, e.g. {}".format(len(miss), miss[:3]))
    if set(tips) != set(lab.index):
        say("[i] tree tips {} | lineage table {} (differences are ignored)".format(len(tips), len(lab)))
    B = ct[tips].apply(lambda c: c.str.upper()).values          # sites x tips
    refs = ct[refcol].str.upper().values
    P = ct["POS"].astype(int).values
    CH = ct["CHR"].values
    rows_meta, rows_M = [], []
    for k in range(len(ct)):
        al = B[k]
        missing = np.isin(al, ["N", "-", "?", ""]) | pd.isna(al)
        for alt in sorted(set(al[~missing]) - {refs[k]}):
            m = (al == alt).astype(np.int8)
            m[missing] = -1
            rows_meta.append(("{}_{}_{}_{}".format(CH[k], P[k], refs[k], alt), P[k], refs[k], alt))
            rows_M.append(m)
    M = np.vstack(rows_M).T                                     # tips x alleles
    meta = pd.DataFrame(rows_meta, columns=["variant", "pos", "ref", "alt"])
    if args.variant not in set(meta["variant"]):
        sys.exit("[!] {} not among the {} core SNP alleles".format(args.variant, len(meta)))
    vi = int(np.where(meta["variant"].values == args.variant)[0][0])
    say("core SNP sites {} -> allele rows {}; tips {}".format(len(ct), len(meta), len(tips)))

    # ------------------------------------------------ 1. locus
    say("\n1. LOCUS")
    seq = {}
    cur = None
    for line in open(args.ref):
        if line.startswith(">"):
            cur = line[1:].split()[0]
            seq[cur] = []
        elif cur:
            seq[cur].append(line.strip().upper())
    seq = {k: "".join(v) for k, v in seq.items()}
    g = seq.get(chrom) or next(iter(seq.values()))
    if g[pos - 1] != vref:
        say("[!] K-10 base at {} is {}, variant says {}".format(pos, g[pos - 1], vref))
    cds = []
    for line in open(args.gff):
        if line.startswith("##FASTA"):
            break
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS":
            continue
        at = dict(x.split("=", 1) for x in f[8].split(";") if "=" in x)
        cds.append((int(f[3]), int(f[4]), f[6], at.get("locus_tag", at.get("ID", "?")),
                    at.get("gene", ""), at.get("product", "")))
    cds.sort()
    hit = [c for c in cds if c[0] <= pos <= c[1]]
    locus = ""
    for s, e, strand, lt, gene, prod in hit:
        if strand == "+":
            off = pos - s
            cs = s + (off // 3) * 3
            codon = g[cs - 1:cs + 2]
            cp = off % 3
            alt_codon = codon[:cp] + valt + codon[cp + 1:]
        else:
            off = e - pos
            ce = e - (off // 3) * 3
            codon = rc(g[ce - 3:ce])
            cp = off % 3
            alt_codon = codon[:cp] + valt.translate(COMP) + codon[cp + 1:]
        aa1, aa2 = CODE.get(codon, "?"), CODE.get(alt_codon, "?")
        kind = "synonymous" if aa1 == aa2 else ("nonsense" if aa2 == "*" else "missense")
        aapos = off // 3 + 1
        locus = "{} {} ({}), strand {}, codon {} -> {}, {}{}{} ({}), codon position {}, CDS {}-{}".format(
            lt, gene or "-", prod, strand, codon, alt_codon, aa1, aapos, aa2, kind, cp + 1, s, e)
        say("  in CDS: " + locus)
    if not hit:
        left = [c for c in cds if c[1] < pos]
        right = [c for c in cds if c[0] > pos]
        lt = left[-1] if left else None
        rt = right[0] if right else None
        locus = "intergenic between {} and {}".format(
            "{} {} ({}, {} bp away)".format(lt[3], lt[4], lt[5], pos - lt[1]) if lt else "-",
            "{} {} ({}, {} bp away)".format(rt[3], rt[4], rt[5], rt[0] - pos) if rt else "-")
        say("  " + locus)
    if Path(args.mask).exists():
        mk = pd.read_csv(args.mask, sep="\t", header=None, names=["chr", "s0", "e", "product"])
        d = np.where(pos < mk["s0"] + 1, mk["s0"] + 1 - pos, np.where(pos > mk["e"], pos - mk["e"], 0))
        j = int(np.argmin(d))
        say("  nearest masked repeat gene: {} at {}-{}, {} bp away".format(
            mk.at[j, "product"], mk.at[j, "s0"] + 1, mk.at[j, "e"], int(d[j])))
        near = mk[d <= args.window]
        if len(near):
            say("  masked genes within +-{} bp: {}".format(args.window, "; ".join(
                "{} ({}-{})".format(r["product"], r["s0"] + 1, r["e"]) for _, r in near.iterrows())))
    else:
        say("  [i] {} not found: repeat proximity skipped".format(args.mask))

    # ------------------------------------------------ carriers
    x = M[:, vi]
    tip = np.array(tips)
    n1, n0, nm = int((x == 1).sum()), int((x == 0).sum()), int((x == -1).sum())
    minor = 1 if n1 <= n0 else 0
    car = set(tip[x == minor])
    non = set(tip[x == 1 - minor])
    minor_allele = valt if minor == 1 else "not-{} (mostly {})".format(valt, vref)
    say("\n2. ALLELES")
    say("  {} carriers {} | {} carriers {} | missing {} -> minor allele: {} ({} genomes)".format(
        valt, n1, "non-" + valt, n0, nm, minor_allele, len(car)))
    tab = pd.DataFrame({"lineage": lab.reindex(tip).fillna("?").values,
                        "state": np.where(x == 1, valt, np.where(x == 0, "other", "missing"))})
    say(pd.crosstab(tab["lineage"], tab["state"]).to_string())
    say("  (the tree is rooted between S and C: the allele fixed in S and C/Bison is"
        " the likely ancestral state)")

    # ------------------------------------------------ tree + Fitch
    lsets = lab.reindex(tip)
    S = [t for t in tips if str(lsets.get(t, "")).startswith("S/")]
    children, parent, post, root, edge = root_on_split(adj, name, S)
    tipidx = {nd: tips.index(nm_) for nd, nm_ in name.items()}
    fpath = base / "fitch_all_sites.tsv"
    if fpath.exists():
        F = pd.read_csv(fpath, sep="\t")
        if len(F) != len(meta) or (F["variant"].values != meta["variant"].values).any():
            F = None
    else:
        F = None
    if F is None:
        say("\n[i] Fitch parsimony on {} allele rows ...".format(len(meta)))
        cost = fitch(children, post, tipidx, M)
        cnt1 = (M == 1).sum(axis=0)
        cnt0 = (M == 0).sum(axis=0)
        F = meta.copy()
        F["minor_count"] = np.minimum(cnt1, cnt0)
        F["fitch_changes"] = cost
        F.to_csv(fpath, sep="\t", index=False)
    say("\n3. FITCH PARSIMONY (minimum number of changes on the ML tree)")
    inf = F[F["minor_count"] >= 2]
    fc = int(F.at[vi, "fitch_changes"])
    say("  this site: {} changes (minor allele count {})".format(fc, int(F.at[vi, "minor_count"])))
    say("  genome-wide, sites with minor count >= 2: {} | 1 change {:.1%} | 2 {:.1%} | 3-5 {:.1%} | 6-20 {:.1%} | >20 {:.1%}".format(
        len(inf), (inf.fitch_changes == 1).mean(), (inf.fitch_changes == 2).mean(),
        inf.fitch_changes.between(3, 5).mean(), inf.fitch_changes.between(6, 20).mean(),
        (inf.fitch_changes > 20).mean()))
    say("  ensemble consistency index (all allele rows): {:.3f}".format(
        len(F) / F["fitch_changes"].clip(lower=1).sum()))
    say("  sites with >= {} changes: {} ({:.2%}); this site ranks {} of {}".format(
        fc, int((inf.fitch_changes >= fc).sum()), (inf.fitch_changes >= fc).mean(),
        int((inf.fitch_changes > fc).sum()) + 1, len(inf)))
    top = inf.sort_values("fitch_changes", ascending=False).head(15)
    say("  most homoplastic sites:\n" + top[["variant", "minor_count", "fitch_changes"]].to_string(index=False))

    # ------------------------------------------------ 4. carrier groups
    say("\n4. INDEPENDENT CARRIER GROUPS (tree rooted on the S|C branch)")
    bit = {t: 1 << k for k, t in enumerate(tips)}
    carm = sum(bit[t] for t in car)
    mism = sum(bit[t] for t in tip[x == -1])
    mask = {}
    for nd in post:
        kids = children.get(nd, [])
        mask[nd] = bit[name[nd]] if not kids else sum(mask[k] for k in kids)

    def allcar(nd):
        m = mask[nd] & ~mism
        return m != 0 and (m & ~carm) == 0

    groups = [nd for nd in post if nd != root and allcar(nd) and not allcar(parent[nd])]
    reg = pd.Series("unknown", index=tips)
    for c in ("Region_Americas", "Region_Europe", "Region_Oceania"):
        if c in tr.columns:
            reg[[t for t in tips if tr[c].get(t) == "1"]] = c.replace("Region_", "")
    study = pd.Series([ena.at[t, "study_accession"] if t in ena.index else "assembly" for t in tips], index=tips)

    def support(nd):
        p = parent.get(nd)
        if p is None or p == root:
            return sup.get(frozenset(edge), "")
        return sup.get(frozenset((nd, p)), "")

    G = []
    for nd in groups:
        mem = [t for t in tips if mask[nd] & bit[t] and t in car]
        p = parent[nd]
        G.append({"group_size": len(mem), "support": support(nd),
                  "parent_support": support(p) if p != root else "root",
                  "regions": "; ".join("{}:{}".format(k, v) for k, v in reg[mem].value_counts().items()),
                  "studies": "; ".join("{}:{}".format(k, v) for k, v in study[mem].value_counts().head(3).items()),
                  "lineages": "; ".join("{}:{}".format(k, v) for k, v in lsets[mem].value_counts().items()),
                  "members": ",".join(mem[:20]) + (",..." if len(mem) > 20 else "")})
    G = pd.DataFrame(G).sort_values("group_size", ascending=False)
    G.to_csv(out / "groups.tsv", sep="\t", index=False)
    sz = G["group_size"]
    say("  groups {} | singletons {} | sizes 2-9 {} | 10-49 {} | >=50 {} | largest {} ({:.0%} of carriers)".format(
        len(G), int((sz == 1).sum()), int(sz.between(2, 9).sum()), int(sz.between(10, 49).sum()),
        int((sz >= 50).sum()), int(sz.max()), sz.max() / len(car)))
    pd.set_option("display.width", 250)
    pd.set_option("display.max_colwidth", 60)
    say(G.drop(columns="members").head(20).to_string(index=False))

    # ------------------------------------------------ 5/6. neighbours and similar patterns
    say("\n5. NEIGHBOURING SNPs (+-{} bp)".format(args.window))
    cvec = np.array([t in car for t in tips])

    def jacc(col):
        m = M[:, col]
        mn = 1 if (m == 1).sum() <= (m == 0).sum() else 0
        s = (m == mn)
        inter = (s & cvec).sum()
        uni = (s | cvec).sum()
        return inter / uni if uni else 0.0, inter / cvec.sum(), (inter / s.sum()) if s.sum() else 0.0

    W = meta.index[(meta["pos"] - pos).abs() <= args.window]
    N = []
    for k in W:
        j, fa, fb = jacc(k)
        N.append({"variant": meta.at[k, "variant"], "distance": int(meta.at[k, "pos"] - pos),
                  "minor_count": int(F.at[k, "minor_count"]), "fitch_changes": int(F.at[k, "fitch_changes"]),
                  "jaccard": round(j, 3), "frac_focal_carriers": round(fa, 3), "frac_own_in_focal": round(fb, 3)})
    N = pd.DataFrame(N)
    N.to_csv(out / "neighbours.tsv", sep="\t", index=False)
    say(N.to_string(index=False) if len(N) else "  no other core SNP in the window")
    linked = N[(N["variant"] != args.variant) & (N["jaccard"] >= args.jaccard)]
    say("  neighbours sharing the carrier set (Jaccard >= {}): {}".format(args.jaccard, len(linked)))

    # SNP density in the window: non-reference calls per genome (focal site excluded)
    wsites = [k for k in range(len(ct)) if abs(P[k] - pos) <= args.window and P[k] != pos]
    if wsites:
        nonref = np.array([[(B[k, t] != refs[k]) and B[k, t] not in ("N", "-") for t in range(len(tips))]
                           for k in wsites]).sum(axis=0)
        C_non = np.array([str(lsets.get(t, "")) == "C/non-B" for t in tips])
        a_ = nonref[cvec & C_non]
        b_ = nonref[(~cvec) & C_non & (M[:, vi] != -1)]
        msg = "  window SNP calls per genome (C/non-B only; {} core sites): carriers median {:.0f} (mean {:.2f}) vs others median {:.0f} (mean {:.2f})".format(
            len(wsites), np.median(a_) if len(a_) else np.nan, a_.mean() if len(a_) else np.nan,
            np.median(b_) if len(b_) else np.nan, b_.mean() if len(b_) else np.nan)
        try:
            from scipy.stats import mannwhitneyu
            if len(a_) and len(b_):
                msg += ", Mann-Whitney p {:.3g}".format(mannwhitneyu(a_, b_).pvalue)
        except ImportError:
            pass
        say(msg)

    say("\n6. GENOME-WIDE SITES WITH A SIMILAR CARRIER SET (Jaccard >= {})".format(args.jaccard))
    sim = []
    for k in range(len(meta)):
        if F.at[k, "minor_count"] < 2:
            continue
        j, fa, fb = jacc(k)
        if j >= args.jaccard:
            sim.append({"variant": meta.at[k, "variant"], "pos": int(meta.at[k, "pos"]),
                        "distance": int(meta.at[k, "pos"] - pos), "jaccard": round(j, 3),
                        "fitch_changes": int(F.at[k, "fitch_changes"])})
    sim = pd.DataFrame(sim)
    sim.to_csv(out / "similar_patterns.tsv", sep="\t", index=False)
    if len(sim):
        say(sim.to_string(index=False))
        say("  similar sites within the window: {} | elsewhere: {}".format(
            int((sim["distance"].abs() <= args.window).sum()) - 1,
            int((sim["distance"].abs() > args.window).sum())))

    # ------------------------------------------------ 7. Gubbins
    say("\n7. GUBBINS BLOCKS overlapping {}-{}".format(pos - args.window, pos + args.window))
    nb = 0
    if Path(args.gubbins).exists():
        for line in open(args.gubbins):
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            s, e = int(f[3]), int(f[4])
            if e < pos - args.window or s > pos + args.window:
                continue
            m = re.search(r'taxa="([^"]*)"', f[8])
            taxa = set(m.group(1).split()) if m else set()
            nb += 1
            where = "contains the site" if s <= pos <= e else "{} bp from the site".format(
                min(abs(s - pos), abs(e - pos)))
            say("  block {}-{} ({} bp, {}): {} taxa, {} of them focal carriers".format(
                s, e, e - s + 1, where, len(taxa), len(taxa & car)))
    say("  blocks in window: {}".format(nb))

    # ------------------------------------------------ 8. reads
    if not args.no_reads and Path(args.reads).is_dir():
        say("\n8. READS: allele fraction and depth ratio (site depth / genome mean depth)")
        wpos = set(P[k] for k in wsites) if wsites else set()
        rows = []
        for t in tips:
            fa, fd = Path(args.reads) / (t + ".ad.tsv"), Path(args.reads) / (t + ".depth")
            if not fa.exists() or not fd.exists():
                continue
            try:
                md = float(open(fd).read().split()[0])
            except (ValueError, IndexError):
                continue
            if md <= 0:
                continue
            af = dep = np.nan
            wr, allr = [], []
            for line in open(fa):
                f = line.rstrip("\n").split("\t")
                if len(f) < 4:
                    continue
                ad = [int(v) for v in f[3].split(",") if v not in (".", "")]
                tot = sum(ad)
                p_ = int(f[0])
                allr.append(tot / md)
                if p_ == pos:
                    alts = f[2].split(",")
                    a = ad[1 + alts.index(valt)] if valt in alts and 1 + alts.index(valt) < len(ad) else 0
                    af = a / tot if tot else np.nan
                    dep = tot / md
                elif p_ in wpos:
                    wr.append(tot / md)
            rows.append({"genome": t, "carrier": t in car, "alt_fraction": af, "depth_ratio_site": dep,
                         "depth_ratio_window_median": np.median(wr) if wr else np.nan,
                         "depth_ratio_genome_median": np.median(allr) if allr else np.nan})
        Rd = pd.DataFrame(rows)
        Rd.to_csv(out / "reads.tsv", sep="\t", index=False)
        for lab_, sub in (("carriers", Rd[Rd["carrier"]]), ("non-carriers", Rd[~Rd["carrier"]])):
            say("  {:<13} n {:>4} | ALT fraction median {:.3f} | depth ratio at site median {:.2f} (p5 {:.2f}, p95 {:.2f}) | window {:.2f} | genome {:.2f}".format(
                lab_, len(sub), sub["alt_fraction"].median(), sub["depth_ratio_site"].median(),
                sub["depth_ratio_site"].quantile(.05), sub["depth_ratio_site"].quantile(.95),
                sub["depth_ratio_window_median"].median(), sub["depth_ratio_genome_median"].median()))
    else:
        say("\n8. READS: skipped")

    say("\nHOW TO READ")
    say("  R recombination: linked neighbours (section 5/6 within the window), more window SNPs in carriers, nearby Gubbins block")
    say("  M recurrent mutation: isolated site, several independent groups, neighbours follow the tree")
    say("  P paralog/mapping: masked repeat close by, depth ratio at the site differs between carriers and non-carriers")
    say("  T tree error: many UNLINKED sites with the same carrier set (section 6, outside the window)")
    (out / "report.txt").write_text("\n".join(R) + "\n")
    say("\nwritten: {}/".format(out))


if __name__ == "__main__":
    main()
