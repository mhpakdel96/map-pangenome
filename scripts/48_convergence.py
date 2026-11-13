#!/usr/bin/env python3
"""
48_convergence.py

Signatures of selection that do not depend on trait labels.

Region and host are almost fully explained by phylogeny (h^2 ~ 1, 27-29), so
association tests cannot reveal adaptation in MAP. In clonal bacteria the
standard alternative is CONVERGENCE: a gene that mutates (or is lost)
independently on many branches of the tree more often than expected is a
candidate target of selection, whatever the trait.

PART A — gene-level SNP convergence (core SNPs, snippy-core core.tab)
  * minimum number of independent mutation events per SNP allele = Fitch
    changes on the ML tree; computed twice: all genomes, and with the
    possible mixed samples (>= --mixed mixed sites, 35) set to missing,
    because mixed samples create spurious extra changes (45). The second
    count is used for the tests.
  * events summed per K-10 CDS (Prokka GFF); CDS overlapping masked
    repeats (PE/PPE, IS, phage; 13) are excluded
  * expected events per gene = total events x gene length / total length of
    the included CDS; one-sided Poisson test, Benjamini-Hochberg FDR
  * per gene: nonsynonymous / synonymous events (relative to K-10 codons),
    number of distinct sites, largest number of events at one site
PART B — recurrent gene loss/gain (read-based ternary calls, 36)
  * for every accessory gene (called in >= 2 genomes in each state): Fitch
    changes of the presence pattern (NA and genomes without reads =
    missing), the closest clade and how many genomes differ from it, where
    the minority-state genomes come from, and their sequencing depth
    (absences at low depth are less trustworthy)
  1 change = lineage-specific gene; >= 3 changes = recurrent loss/gain

Caveats (state in Methods): Fitch gives a MINIMUM count; topological error
(UFBoot median 86) can add changes at clade-level sites; the per-gene
expectation assumes a uniform rate along the genome (no correction for
GC/context); callable length is approximated by CDS length.

  python 48_convergence.py
Outputs: 06_pangwas/convergence/{gene_snp_convergence.tsv, allele_events.tsv,
         gene_loss_convergence.tsv, report.txt}
Environment: ml_env. Imports 44_homoplasy_site.py.
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mixed", type=int, default=20)
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--tree", default=str(ROOT / "05_phylogeny/iqtree/tree_final.nwk"))
    ap.add_argument("--coretab", default=str(ROOT / "05_phylogeny/core_masked/core.tab"))
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
    ap.add_argument("--out", default=str(ROOT / "06_pangwas/convergence"))
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
    lin = pd.read_csv(args.lineages, sep="\t", dtype=str).set_index("genome")
    lab = (lin["lineage"] + "/" + lin["sublineage"]).reindex(tips).fillna("?")
    S = [t for t in tips if lab[t].startswith("S/")]
    children, parent, post, root, edge = s44.root_on_split(adj, name, S)
    tipidx = {nd: tips.index(nm) for nd, nm in name.items()}
    q = pd.read_csv(args.readqc, sep="\t", index_col=0) if Path(args.readqc).exists() else pd.DataFrame()
    mixed = [t for t in tips if "mixed_sites" in q.columns and q["mixed_sites"].get(t, 0) >= args.mixed]
    say("tree: {} tips | possible mixed samples set to missing in PART A: {}".format(len(tips), len(mixed)))

    # ------------------------------------------------ K-10 annotation
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
    masked = []
    if Path(args.mask).exists():
        mk = pd.read_csv(args.mask, sep="\t", header=None, names=["chr", "s0", "e", "product"])
        masked = list(zip(mk["s0"] + 1, mk["e"]))

    def overlaps_mask(s, e):
        return any(not (e < a or s > b) for a, b in masked)

    C = pd.DataFrame(cds, columns=["start", "end", "strand", "locus_tag", "gene", "product"])
    C["length"] = C["end"] - C["start"] + 1
    C["masked"] = [overlaps_mask(s, e) for s, e in zip(C["start"], C["end"])]
    starts = C["start"].tolist()

    def cds_of(pos):
        hits = []
        i = bisect.bisect_right(starts, pos) - 1
        while i >= 0 and C.at[i, "start"] > pos - 20000:
            if C.at[i, "start"] <= pos <= C.at[i, "end"]:
                hits.append(i)
            i -= 1
        return hits

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
        if "?" in (a1, a2):
            return "unknown"
        return "S" if a1 == a2 else "N"

    # ------------------------------------------------ PART A
    ct = pd.read_csv(args.coretab, sep="\t", dtype=str)
    refcol = "REF" if "REF" in ct.columns else "Reference"
    BA = ct[tips].apply(lambda c: c.str.upper()).fillna("N").values
    P = ct["POS"].astype(int).values
    refs = ct[refcol].str.upper().values
    meta, cols = [], []
    for k in range(len(ct)):
        row = BA[k]
        miss = np.isin(row, ["N", "-", "?"])
        for alt in sorted(set(row[~miss]) - {refs[k]}):
            m = (row == alt).astype(np.int8)
            m[miss] = -1
            meta.append((P[k], refs[k], alt))
            cols.append(m)
    M = np.vstack(cols).T
    say("core SNP sites {} -> allele rows {}".format(len(ct), M.shape[1]))
    f_all = s44.fitch(children, post, tipidx, M)
    M2 = M.copy()
    for t in mixed:
        M2[tips.index(t), :] = -1
    f_nm = s44.fitch(children, post, tipidx, M2)
    A = pd.DataFrame(meta, columns=["pos", "ref", "alt"])
    A["events_all"] = f_all
    A["events"] = f_nm
    A["carriers"] = (M == 1).sum(axis=0)
    A["cds"] = [cds_of(p) for p in A["pos"]]
    A["effect"] = [",".join(effect(i, p, a) for i in c) if c else "intergenic"
                   for c, p, a in zip(A["cds"], A["pos"], A["alt"])]
    say("events (Fitch, minimum): all genomes {} | mixed samples as missing {} | alleles with >= 2 events: {} -> {}".format(
        int(f_all.sum()), int(f_nm.sum()), int((f_all >= 2).sum()), int((f_nm >= 2).sum())))
    A2 = A.explode("cds")
    A2 = A2[A2["cds"].notna()]
    A2["cds"] = A2["cds"].astype(int)
    inc = C.index[~C["masked"]]
    A2 = A2[A2["cds"].isin(inc)]
    Ltot = int(C.loc[inc, "length"].sum())
    Etot = int(A2["events"].sum())
    say("CDS included {} of {} (masked excluded), total length {:,} bp, events in them {}".format(
        len(inc), len(C), Ltot, Etot))
    A2["eff1"] = [effect(i, p, a) for i, p, a in zip(A2["cds"], A2["pos"], A2["alt"])]
    nN = int(A2.loc[A2["eff1"] == "N", "events"].sum())
    nS = int(A2.loc[A2["eff1"] == "S", "events"].sum())
    say("genome-wide events: nonsynonymous {} | synonymous {} | N/S {:.2f}".format(nN, nS, nN / max(nS, 1)))
    grp = A2.groupby("cds")
    G = C.loc[inc, ["locus_tag", "gene", "product", "start", "end", "strand", "length"]].copy()
    G["sites"] = grp["pos"].nunique().reindex(G.index).fillna(0).astype(int)
    G["events"] = grp["events"].sum().reindex(G.index).fillna(0).astype(int)
    G["events_all"] = grp["events_all"].sum().reindex(G.index).fillna(0).astype(int)
    G["N_events"] = A2[A2["eff1"] == "N"].groupby("cds")["events"].sum().reindex(G.index).fillna(0).astype(int)
    G["S_events"] = A2[A2["eff1"] == "S"].groupby("cds")["events"].sum().reindex(G.index).fillna(0).astype(int)
    G["max_events_one_site"] = grp["events"].max().reindex(G.index).fillna(0).astype(int)
    G["sites_multi_event"] = A2[A2["events"] >= 2].groupby("cds")["pos"].nunique().reindex(G.index).fillna(0).astype(int)
    G["expected"] = Etot * G["length"] / Ltot
    G["ratio"] = (G["events"] / G["expected"]).round(2)
    G["p"] = poisson.sf(G["events"] - 1, G["expected"])
    G["q"] = bh(G["p"].values)
    G["expected"] = G["expected"].round(2)
    G = G.sort_values("p")
    G.to_csv(out / "gene_snp_convergence.tsv", sep="\t", index=False)
    A.drop(columns="cds").to_csv(out / "allele_events.tsv", sep="\t", index=False)
    say("\nPART A — genes with excess independent mutation events (Poisson vs length, BH q < {}): {}".format(
        args.fdr, int((G["q"] < args.fdr).sum())))
    say(G.head(25)[["locus_tag", "gene", "product", "length", "sites", "events", "expected", "ratio",
                    "N_events", "S_events", "max_events_one_site", "sites_multi_event", "p", "q"]].to_string(index=False))
    sig = G[G["q"] < args.fdr]
    for _, r in sig.head(15).iterrows():
        a = A2[A2["cds"] == r.name].sort_values("pos")
        say("  {} {}: ".format(r["locus_tag"], r["gene"] or "-") + "; ".join(
            "{}{}>{} {} ev{} car{}".format(p, rf, al, ef, ev, ca) for p, rf, al, ef, ev, ca in
            zip(a["pos"], a["ref"], a["alt"], a["eff1"], a["events"], a["carriers"])))

    # ------------------------------------------------ PART B
    if Path(args.ternary).exists():
        T3 = pd.read_csv(args.ternary, sep="\t", index_col=0)
        GT = pd.read_csv(args.genetable, sep="\t").set_index("gene") if Path(args.genetable).exists() else None
        tr = pd.read_csv(args.traits, dtype=str).set_index("Name")
        reg = pd.Series("unknown", index=tips)
        for c in ("Region_Americas", "Region_Europe", "Region_Oceania"):
            if c in tr.columns:
                reg[[t for t in tips if tr[c].get(t) == "1"]] = c.replace("Region_", "")
        ena = pd.read_csv(args.ena, sep="\t", dtype=str).drop_duplicates("run_accession").set_index("run_accession")
        study = pd.Series([ena.at[t, "study_accession"] if t in ena.index else "assembly" for t in tips], index=tips)
        depth = q["depth"] if "depth" in q.columns else pd.Series(dtype=float)
        bit = {t: 1 << k for k, t in enumerate(tips)}
        mask = {}
        for nd in post:
            kids = children.get(nd, [])
            mask[nd] = bit[name[nd]] if not kids else sum(mask[k] for k in kids)
        X = pd.DataFrame(-1, index=tips, columns=T3.index, dtype=np.int8)
        common = [t for t in tips if t in T3.columns]
        sub = T3[common].T
        X.loc[common] = sub.fillna(-1).astype(np.int8).values
        n1 = (X == 1).sum()
        n0 = (X == 0).sum()
        keep = X.columns[(n1 >= 2) & (n0 >= 2)]
        X = X[keep]
        fc = s44.fitch(children, post, tipidx, X.values)
        rows = []
        for j, gname in enumerate(keep):
            x = X[gname]
            minor = 0 if n0[gname] <= n1[gname] else 1
            mt = [t for t in tips if x[t] == minor]
            mm = sum(bit[t] for t in mt)
            km = sum(bit[t] for t in tips if x[t] != -1)
            best = min(((((m ^ mm) & km).bit_count(), m) for nd, m in mask.items() if nd != root),
                       key=lambda z: z[0])
            rows.append({"gene": gname, "annotation": GT.at[gname, "annotation"] if GT is not None and "annotation" in GT.columns and gname in GT.index else "",
                         "len_bp": GT.at[gname, "len_bp"] if GT is not None and gname in GT.index else np.nan,
                         "minor_state": "absent" if minor == 0 else "present",
                         "n_minor": len(mt), "n_called": int((x != -1).sum()),
                         "fitch_changes": int(fc[j]), "closest_clade_size": bin(best[1]).count("1"),
                         "differ_from_clade": best[0],
                         "minor_depth_median": round(float(depth.reindex(mt).median()), 1) if len(depth) else np.nan,
                         "minor_lineages": "; ".join("{}:{}".format(k, v) for k, v in lab[mt].value_counts().head(3).items()),
                         "minor_regions": "; ".join("{}:{}".format(k, v) for k, v in reg[mt].value_counts().head(3).items()),
                         "minor_studies": "; ".join("{}:{}".format(k, v) for k, v in study[mt].value_counts().head(3).items())})
        L = pd.DataFrame(rows).sort_values(["fitch_changes", "n_minor"], ascending=[False, False])
        L.to_csv(out / "gene_loss_convergence.tsv", sep="\t", index=False)
        say("\nPART B — accessory genes (read ternary calls), changes on the tree:")
        say("  genes tested {} | 1 change (lineage-specific) {} | 2 {} | 3-5 {} | >= 6 {}".format(
            len(L), int((L.fitch_changes == 1).sum()), int((L.fitch_changes == 2).sum()),
            int(L.fitch_changes.between(3, 5).sum()), int((L.fitch_changes >= 6).sum())))
        say("  genome median K-10 depth: {:.1f}".format(float(depth.median())) if len(depth) else "")
        say(L.head(30).to_string(index=False))
        say("\n  lineage-specific genes (1 change), by lineage of the minority state:")
        say(L[L.fitch_changes == 1][["gene", "annotation", "minor_state", "n_minor", "minor_lineages"]]
            .to_string(index=False))
    else:
        say("\nPART B skipped: {} not found".format(args.ternary))
    (out / "report.txt").write_text("\n".join(R) + "\n")
    say("\nwritten: {}/".format(out))


if __name__ == "__main__":
    main()
