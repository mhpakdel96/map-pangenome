#!/usr/bin/env python3
"""
49_convergence_verified.py

Gene-level convergence, version 2: read-verified genotypes and events
counted per BRANCH.

Why a second version. 48 found 83 genes with excess SNPs, but many of them
are clusters of several SNPs within a few bp carried by the same one genome
(e.g. araQ_1, deaD_1, carD, yhdG_2) or the same small clade (menH_2). Such
clusters are one event (or an assembly/alignment artefact from snippy
--ctgs on fragmented assemblies), not repeated independent mutation. And
48 PART B's 'recurrent losses' were short genes absent mostly in genomes
with low read depth (median 18-30x vs 53x genome-wide): coverage dropout.

Changes here
  1. READ VERIFICATION of every assembly-based call: for each genome with
     reads, the allele with >= --min-af of >= --min-dp reads at that site
     (08_reads/k10/<ID>.ad.tsv) must equal the assembly call; otherwise the
     call is set to missing (discordant or mixed/low depth). Genomes without
     reads (NCBI assemblies) keep their assembly calls.
  2. Fitch DOWN-pass + top-down reconstruction places every change on a
     branch; per gene the number of DISTINCT branches with >= 1 change is
     the convergence count (several SNPs on one branch = one event).
     Poisson test vs gene length, BH FDR.
  3. PART B: absence calls from genomes below --min-depth x K-10 depth are
     set to NA before counting gene gain/loss changes.

  python 49_convergence_verified.py
Outputs: 06_pangwas/convergence_v2/{call_verification.tsv, site_events.tsv,
         gene_branch_convergence.tsv, gene_loss_depthfiltered.tsv, report.txt}
Environment: ml_env. Imports 44_homoplasy_site.py. Parsing the read tables
takes several minutes (srun, 2 cpu, 24 GB).
"""

import argparse
import bisect
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import poisson

ROOT = Path("/cfs/earth/scratch/tabasran/MAP")


def load44(path):
    sp = importlib.util.spec_from_file_location("s44", path)
    m = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(m)
    return m


def bh(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    o = np.argsort(p)
    q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1.0)


def reconstruct(children, post, root, tipidx, M):
    """Fitch down-pass + top-down assignment (parent state kept when allowed).
    M: tips x alleles int8 (1/0/-1). Returns (cost per allele,
    {node: bool array of alleles changing on the branch above node})."""
    s0, s1 = {}, {}
    cost = np.zeros(M.shape[1], dtype=np.int32)
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
            cost += e
            a0 = np.where(e, a0 | b0, i0)
            a1 = np.where(e, a1 | b1, i1)
        s0[x], s1[x] = a0, a1
    state = {root: np.where(s0[root], 0, 1).astype(np.int8)}
    change = {}
    stack = [root]
    while stack:
        p = stack.pop()
        for c in children.get(p, []):
            st = state[p]
            ok = np.where(st == 0, s0[c], s1[c])
            state[c] = np.where(ok, st, 1 - st).astype(np.int8)
            change[c] = state[c] != st
            stack.append(c)
    return cost, change


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-dp", type=int, default=10)
    ap.add_argument("--min-af", type=float, default=0.8)
    ap.add_argument("--min-depth", type=float, default=30.0, help="PART B: absences below this K-10 depth -> NA")
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--tree", default=str(ROOT / "05_phylogeny/iqtree/tree_final.nwk"))
    ap.add_argument("--coretab", default=str(ROOT / "05_phylogeny/core_masked/core.tab"))
    ap.add_argument("--reads", default=str(ROOT / "08_reads/k10"))
    ap.add_argument("--mask", default=str(ROOT / "05_phylogeny/core_masked/mask_genes.tsv"))
    ap.add_argument("--ref", default=str(ROOT / "03_genomes/fasta/GCA_000007865.1.fasta"))
    ap.add_argument("--gff", default=str(ROOT / "03_genomes/gff/GCA_000007865.1.gff"))
    ap.add_argument("--readqc", default=str(ROOT / "08_reads/summary/genome_qc_reads.tsv"))
    ap.add_argument("--ternary", default=str(ROOT / "08_reads/summary/presence_reads_ternary.tsv"))
    ap.add_argument("--genetable", default=str(ROOT / "04_pangenome/dna_presence/gene_table.tsv"))
    ap.add_argument("--lineages", default=str(ROOT / "05_phylogeny/typing/lineage_assignments.tsv"))
    ap.add_argument("--traits", default=str(ROOT / "06_pangwas/trait.csv"))
    ap.add_argument("--ena", default=str(ROOT / "01_metadata/ena_raw.tsv"))
    ap.add_argument("--s44", default=str(Path(__file__).resolve().parent / "44_homoplasy_site.py"))
    ap.add_argument("--out", default=str(ROOT / "06_pangwas/convergence_v2"))
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    s44 = load44(args.s44)
    R = []

    def say(s=""):
        print(s, flush=True)
        R.append(s)

    pd.set_option("display.width", 260)
    pd.set_option("display.max_colwidth", 45)

    # ------------------------------------------------ tree
    adj, name, sup = s44.parse_newick(open(args.tree).read())
    tips = sorted(name.values())
    ti = {t: k for k, t in enumerate(tips)}
    lin = pd.read_csv(args.lineages, sep="\t", dtype=str).set_index("genome")
    lab = (lin["lineage"] + "/" + lin["sublineage"]).reindex(tips).fillna("?")
    S = [t for t in tips if lab[t].startswith("S/")]
    children, parent, post, root, edge = s44.root_on_split(adj, name, S)
    tipidx = {nd: ti[nm] for nd, nm in name.items()}
    nleaf = {}
    for nd in post:
        kids = children.get(nd, [])
        nleaf[nd] = 1 if not kids else sum(nleaf[k] for k in kids)

    # ------------------------------------------------ 1. read verification
    ct = pd.read_csv(args.coretab, sep="\t", dtype=str)
    refcol = "REF" if "REF" in ct.columns else "Reference"
    P = ct["POS"].astype(int).values
    refs = ct[refcol].str.upper().values
    BA = ct[tips].apply(lambda c: c.str.upper()).fillna("N").values      # sites x tips
    pos_idx = {p: k for k, p in enumerate(P)}
    ver = BA.copy()
    stat = []
    rd = Path(args.reads)
    for t in tips:
        f = rd / (t + ".ad.tsv")
        if not f.exists():
            stat.append({"genome": t, "reads": False})
            continue
        rc = {}
        for line in open(f):
            x = line.rstrip("\n").split("\t")
            if len(x) < 4 or len(x[1]) != 1:
                continue
            p = int(x[0])
            if p not in pos_idx or p in rc:
                continue
            al = [x[1].upper()] + [a.upper() for a in x[2].split(",")]
            ad = [int(v) if v not in (".", "") else 0 for v in x[3].split(",")]
            tot = sum(ad)
            if tot < args.min_dp:
                rc[p] = None
                continue
            b = int(np.argmax(ad))
            rc[p] = al[b] if b < len(al) and len(al[b]) == 1 and ad[b] / tot >= args.min_af else None
        j = ti[t]
        n_ok = n_dis = n_amb = n_nocov = 0
        for k, p in enumerate(P):
            a = BA[k, j]
            if a in ("N", "-", "?"):
                continue
            r = rc.get(p, "NOLINE")
            if r == "NOLINE":
                n_nocov += 1
                ver[k, j] = "N"
            elif r is None:
                n_amb += 1
                ver[k, j] = "N"
            elif r != a:
                n_dis += 1
                ver[k, j] = "N"
            else:
                n_ok += 1
        stat.append({"genome": t, "reads": True, "concordant": n_ok, "discordant": n_dis,
                     "ambiguous_or_lowdp": n_amb, "no_read_line": n_nocov})
    ST = pd.DataFrame(stat).set_index("genome")
    ST.to_csv(out / "call_verification.tsv", sep="\t")
    rr = ST[ST["reads"]]
    say("1. READ VERIFICATION of assembly calls at {} core sites (reads >= {}, allele fraction >= {})".format(
        len(P), args.min_dp, args.min_af))
    say("  genomes with reads {} | without {} (assembly calls kept)".format(len(rr), int((~ST['reads']).sum())))
    tot = rr[["concordant", "discordant", "ambiguous_or_lowdp", "no_read_line"]].sum()
    say("  calls: concordant {:,} | discordant {:,} | ambiguous/low depth {:,} | no read line {:,}".format(*map(int, tot)))
    say("  genomes with most discordant calls:\n" + rr.sort_values("discordant", ascending=False).head(10).to_string())

    # ------------------------------------------------ 2. alleles, Fitch, branch placement
    meta, cols, cols_raw = [], [], []
    for k in range(len(P)):
        row, raw = ver[k], BA[k]
        miss = np.isin(row, ["N", "-", "?"])
        alts = sorted(set(raw[~np.isin(raw, ["N", "-", "?"])]) - {refs[k]})
        for alt in alts:
            m = (row == alt).astype(np.int8)
            m[miss] = -1
            mr = (raw == alt).astype(np.int8)
            mr[np.isin(raw, ["N", "-", "?"])] = -1
            meta.append((P[k], refs[k], alt))
            cols.append(m)
            cols_raw.append(mr)
    M = np.vstack(cols).T
    Mraw = np.vstack(cols_raw).T
    cost_raw = s44.fitch(children, post, tipidx, Mraw)
    cost, change = reconstruct(children, post, root, tipidx, M)
    A = pd.DataFrame(meta, columns=["pos", "ref", "alt"])
    A["carriers_raw"] = (Mraw == 1).sum(axis=0)
    A["carriers"] = (M == 1).sum(axis=0)
    A["fitch_raw"] = cost_raw
    A["fitch"] = cost
    A["lost"] = (A["carriers"] == 0) & (A["carriers_raw"] > 0)
    say("\n2. FITCH on read-verified calls")
    say("  alleles {} | no carrier left after verification {} | events raw {} -> verified {}".format(
        len(A), int(A["lost"].sum()), int(cost_raw.sum()), int(cost.sum())))
    inf = A[A["carriers"] >= 2]
    say("  verified alleles with >= 2 carriers: {} | 1 change {:.1%} | >= 2 {} | >= 3 {} | ensemble CI {:.3f}".format(
        len(inf), (inf.fitch == 1).mean(), int((inf.fitch >= 2).sum()), int((inf.fitch >= 3).sum()),
        int((A.carriers > 0).sum()) / max(int(A.loc[A.carriers > 0, "fitch"].clip(lower=1).sum()), 1)))
    say("  most homoplastic verified sites:\n" + A[A.carriers > 0].sort_values("fitch", ascending=False).head(15)
        [["pos", "ref", "alt", "carriers_raw", "carriers", "fitch_raw", "fitch"]].to_string(index=False))
    # branch of each change
    nodes = [nd for nd in post if nd != root]
    CH = np.vstack([change[nd] for nd in nodes])                 # nodes x alleles
    term = np.array([nd in name for nd in nodes])

    # ------------------------------------------------ 3. per gene
    seq, cur = {}, None
    for line in open(args.ref):
        if line.startswith(">"):
            cur = line[1:].split()[0]
            seq[cur] = []
        elif cur:
            seq[cur].append(line.strip().upper())
    g = "".join(next(iter(seq.values())))
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
    cds.sort()
    C = pd.DataFrame(cds, columns=["start", "end", "strand", "locus_tag", "gene", "product"])
    C["length"] = C["end"] - C["start"] + 1
    masked = []
    if Path(args.mask).exists():
        mk = pd.read_csv(args.mask, sep="\t", header=None, names=["chr", "s0", "e", "product"])
        masked = list(zip(mk["s0"] + 1, mk["e"]))
    C["masked"] = [any(not (e < a or s > b) for a, b in masked) for s, e in zip(C["start"], C["end"])]
    starts = C["start"].tolist()

    def cds_of(pos):
        out_, i = [], bisect.bisect_right(starts, pos) - 1
        while i >= 0 and C.at[i, "start"] > pos - 20000:
            if C.at[i, "start"] <= pos <= C.at[i, "end"] and not C.at[i, "masked"]:
                out_.append(i)
            i -= 1
        return out_

    def effect(i, pos, alt):
        s, e, strand = C.at[i, "start"], C.at[i, "end"], C.at[i, "strand"]
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
        return "?" if "?" in (a1, a2) else ("S" if a1 == a2 else "N")

    A["cds"] = [cds_of(p) for p in A["pos"]]
    inc = C.index[~C["masked"]]
    gene_alleles = {}
    for k, cl in enumerate(A["cds"]):
        for i in cl:
            gene_alleles.setdefault(i, []).append(k)
    Ltot = int(C.loc[inc, "length"].sum())
    recs = []
    for i in inc:
        ks = gene_alleles.get(i, [])
        if ks:
            br = CH[:, ks].any(axis=1)
            nb = int(br.sum())
            nterm = int((br & term).sum())
            per_branch = CH[:, ks].sum(axis=1)
            maxsnp = int(per_branch.max())
            effs = [effect(i, int(A.at[k, "pos"]), A.at[k, "alt"]) for k in ks]
            nN = sum(int(A.at[k, "fitch"]) for k, e in zip(ks, effs) if e == "N")
            nS = sum(int(A.at[k, "fitch"]) for k, e in zip(ks, effs) if e == "S")
            sites = len(set(int(A.at[k, "pos"]) for k in ks if A.at[k, "carriers"] > 0))
            # branch sizes (tips below) of the event branches
            sizes = sorted((nleaf[nodes[b]] for b in np.where(br)[0]), reverse=True)
        else:
            nb = nterm = maxsnp = nN = nS = sites = 0
            sizes = []
        recs.append({"cds": i, "branch_events": nb, "terminal_events": nterm, "sites": sites,
                     "max_snps_one_branch": maxsnp, "N_events": nN, "S_events": nS,
                     "event_clade_sizes": ",".join(map(str, sizes[:12]))})
    G = pd.DataFrame(recs).set_index("cds")
    G = C.loc[inc, ["locus_tag", "gene", "product", "start", "end", "length"]].join(G)
    Etot = int(G["branch_events"].sum())
    G["expected"] = Etot * G["length"] / Ltot
    G["ratio"] = (G["branch_events"] / G["expected"]).round(2)
    G["p"] = poisson.sf(G["branch_events"] - 1, G["expected"])
    G["q"] = bh(G["p"].values)
    G["expected"] = G["expected"].round(2)
    G = G.sort_values("p")
    G.to_csv(out / "gene_branch_convergence.tsv", sep="\t", index=False)
    A.assign(cds=A["cds"].map(lambda c: ",".join(C.loc[c, "locus_tag"]) if c else "")).to_csv(
        out / "site_events.tsv", sep="\t", index=False)
    say("\n3. GENE-LEVEL CONVERGENCE: distinct branches with >= 1 change per gene (verified calls)")
    say("  CDS {} (masked excluded) | total length {:,} | branch-events {} | genes with q < {}: {}".format(
        len(inc), Ltot, Etot, args.fdr, int((G["q"] < args.fdr).sum())))
    say(G.head(30)[["locus_tag", "gene", "product", "length", "sites", "branch_events", "terminal_events",
                    "expected", "ratio", "N_events", "S_events", "max_snps_one_branch", "event_clade_sizes",
                    "p", "q"]].to_string(index=False))

    # ------------------------------------------------ 4. PART B with depth filter
    if Path(args.ternary).exists():
        q = pd.read_csv(args.readqc, sep="\t", index_col=0)
        T3 = pd.read_csv(args.ternary, sep="\t", index_col=0)
        GT = pd.read_csv(args.genetable, sep="\t").set_index("gene") if Path(args.genetable).exists() else None
        tr = pd.read_csv(args.traits, dtype=str).set_index("Name")
        reg = pd.Series("unknown", index=tips)
        for c in ("Region_Americas", "Region_Europe", "Region_Oceania"):
            if c in tr.columns:
                reg[[t for t in tips if tr[c].get(t) == "1"]] = c.replace("Region_", "")
        ena = pd.read_csv(args.ena, sep="\t", dtype=str).drop_duplicates("run_accession").set_index("run_accession")
        study = pd.Series([ena.at[t, "study_accession"] if t in ena.index else "assembly" for t in tips], index=tips)
        common = [t for t in tips if t in T3.columns]
        X = pd.DataFrame(-1, index=tips, columns=T3.index, dtype=np.int8)
        X.loc[common] = T3[common].T.fillna(-1).astype(np.int8).values
        low = [t for t in common if q["depth"].get(t, 0) < args.min_depth]
        Xl = X.loc[low]
        X.loc[low] = Xl.where(Xl != 0, -1).values
        say("\n4. GENE GAIN/LOSS on read calls; absences from {} genomes with depth < {}x set to NA".format(
            len(low), args.min_depth))
        n1, n0 = (X == 1).sum(), (X == 0).sum()
        keep = X.columns[(n1 >= 2) & (n0 >= 2)]
        X = X[keep]
        fc = s44.fitch(children, post, tipidx, X.values)
        bit = {t: 1 << k for k, t in enumerate(tips)}
        msk = {}
        for nd in post:
            kids = children.get(nd, [])
            msk[nd] = bit[name[nd]] if not kids else sum(msk[k] for k in kids)
        rows = []
        for j, gn in enumerate(keep):
            x = X[gn]
            minor = 0 if n0[gn] <= n1[gn] else 1
            mt = [t for t in tips if x[t] == minor]
            mm = sum(bit[t] for t in mt)
            km = sum(bit[t] for t in tips if x[t] != -1)
            best = min(((((m ^ mm) & km).bit_count(), m) for nd, m in msk.items() if nd != root), key=lambda z: z[0])
            rows.append({"gene": gn,
                         "annotation": GT.at[gn, "annotation"] if GT is not None and gn in GT.index and "annotation" in GT.columns else "",
                         "len_bp": GT.at[gn, "len_bp"] if GT is not None and gn in GT.index else np.nan,
                         "minor_state": "absent" if minor == 0 else "present", "n_minor": len(mt),
                         "n_called": int((x != -1).sum()), "fitch_changes": int(fc[j]),
                         "closest_clade_size": bin(best[1]).count("1"), "differ_from_clade": best[0],
                         "minor_depth_median": round(float(q["depth"].reindex(mt).median()), 1),
                         "minor_lineages": "; ".join("{}:{}".format(k, v) for k, v in lab[mt].value_counts().head(3).items()),
                         "minor_regions": "; ".join("{}:{}".format(k, v) for k, v in reg[mt].value_counts().head(3).items()),
                         "minor_studies": "; ".join("{}:{}".format(k, v) for k, v in study[mt].value_counts().head(3).items())})
        L = pd.DataFrame(rows).sort_values(["fitch_changes", "n_minor"], ascending=[False, False])
        L.to_csv(out / "gene_loss_depthfiltered.tsv", sep="\t", index=False)
        say("  genes tested {} | 1 change {} | 2 {} | 3-5 {} | >= 6 {}".format(
            len(L), int((L.fitch_changes == 1).sum()), int((L.fitch_changes == 2).sum()),
            int(L.fitch_changes.between(3, 5).sum()), int((L.fitch_changes >= 6).sum())))
        say("  genes with >= 3 changes:\n" + L[L.fitch_changes >= 3].head(30).to_string(index=False))
        say("\n  lineage/clade-specific genes (1 change, or 2 changes with <= 1 genome off the clade):")
        ls = L[(L.fitch_changes == 1) | ((L.fitch_changes == 2) & (L.differ_from_clade <= 1))]
        say(ls[["gene", "annotation", "len_bp", "minor_state", "n_minor", "closest_clade_size",
                "minor_lineages", "minor_regions", "minor_studies"]].to_string(index=False))
    (out / "report.txt").write_text("\n".join(R) + "\n")
    say("\nwritten: {}/".format(out))


if __name__ == "__main__":
    main()
