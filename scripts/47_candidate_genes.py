#!/usr/bin/env python3
"""
47_candidate_genes.py

Manual check of the two gene candidates left from the Scoary analysis
(20-23): desA1_3 (acyl-ACP desaturase; more frequent in cattle isolates) and
glgX (glycogen debranching enzyme; depleted in the Americas). Both were
contig-edge genes, so their annotation-based presence was suspect. Here they
are checked with the assembly-independent, read-based calls (36), which
passed the negative control best (23-34 false hits vs 115 for union80).

Per gene (Panaroo cluster whose name contains the gene as a '~~~' token):
  * gene_table row (length, annotation vs DNA frequency, kept/redundant)
  * read breadth distribution (08_reads/summary/breadth.tsv.gz): full
    (>= 80%), partial (20-80%) or absent (< 20%) — partial means a deletion
    of part of the gene, not of the whole gene
  * association with its trait within C/non-B, on read calls (NA excluded):
    all genomes, MiSeq only, and inside each study with >= 5 genomes on each
    side (Fisher)
  * on the tree: Fitch changes of the presence pattern and the closest clade
    (the genomes that lack the gene: one clade = lineage marker, many
    independent losses = candidate repeated loss)
  * who lacks it: region, host, study, lineage

  python 47_candidate_genes.py [--pairs desA1_3:hostCattle,glgX:Region_Americas]
Outputs: 06_pangwas/candidates/<gene>.tsv (per genome) and report.txt
Environment: ml_env. Imports 44_homoplasy_site.py.
"""

import argparse
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

ROOT = Path("/cfs/earth/scratch/tabasran/MAP")


def load44(path):
    sp = importlib.util.spec_from_file_location("s44", path)
    m = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(m)
    return m


def fisher(pres, trait):
    """pres: 1/0 Series (no NA); trait: 1/0 Series aligned."""
    a = int(((pres == 1) & (trait == 1)).sum())
    b = int(((pres == 0) & (trait == 1)).sum())
    c = int(((pres == 1) & (trait == 0)).sum())
    d = int(((pres == 0) & (trait == 0)).sum())
    p = fisher_exact([[a, b], [c, d]])[1] if min(a + b, c + d) > 0 else np.nan
    fpos = a / (a + b) if a + b else np.nan
    fneg = c / (c + d) if c + d else np.nan
    return a + b, c + d, fpos, fneg, p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", default="desA1_3:hostCattle,glgX:Region_Americas")
    ap.add_argument("--min-side", type=int, default=5)
    ap.add_argument("--reads", default=str(ROOT / "08_reads/summary"))
    ap.add_argument("--genetable", default=str(ROOT / "04_pangenome/dna_presence/gene_table.tsv"))
    ap.add_argument("--union80", default=str(ROOT / "04_pangenome/dna_presence/presence_final.Rtab"))
    ap.add_argument("--tree", default=str(ROOT / "05_phylogeny/iqtree/tree_final.nwk"))
    ap.add_argument("--lineages", default=str(ROOT / "05_phylogeny/typing/lineage_assignments.tsv"))
    ap.add_argument("--traits", default=str(ROOT / "06_pangwas/trait.csv"))
    ap.add_argument("--ena", default=str(ROOT / "01_metadata/ena_raw.tsv"))
    ap.add_argument("--s44", default=str(Path(__file__).resolve().parent / "44_homoplasy_site.py"))
    ap.add_argument("--out", default=str(ROOT / "06_pangwas/candidates"))
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    s44 = load44(args.s44)
    R = []

    def say(s=""):
        print(s, flush=True)
        R.append(s)

    B = pd.read_csv(Path(args.reads) / "breadth.tsv.gz", sep="\t", index_col=0)
    T3 = pd.read_csv(Path(args.reads) / "presence_reads_ternary.tsv", sep="\t", index_col=0)
    GT = pd.read_csv(args.genetable, sep="\t").set_index("gene")
    U = pd.read_csv(args.union80, sep="\t", index_col=0) if Path(args.union80).exists() else None
    tr = pd.read_csv(args.traits, dtype=str).set_index("Name")
    lin = pd.read_csv(args.lineages, sep="\t", dtype=str).set_index("genome")
    lab = lin["lineage"] + "/" + lin["sublineage"]
    ena = pd.read_csv(args.ena, sep="\t", dtype=str).drop_duplicates("run_accession").set_index("run_accession")
    G = list(B.columns)
    plat = pd.Series([ena.at[g, "instrument_model"] if g in ena.index else "assembly" for g in G], index=G)
    study = pd.Series([ena.at[g, "study_accession"] if g in ena.index else "assembly" for g in G], index=G)
    reg = pd.Series("unknown", index=tr.index)
    for c in ("Region_Americas", "Region_Europe", "Region_Oceania"):
        reg[tr[c] == "1"] = c.replace("Region_", "")
    host = pd.Series("NA", index=tr.index)
    if "hostCattle" in tr.columns:
        host[tr["hostCattle"] == "1"] = "cattle"
        host[tr["hostCattle"] == "0"] = "non-cattle"
    cn = [g for g in G if lab.get(g) == "C/non-B"]

    adj, name, sup = s44.parse_newick(open(args.tree).read())
    tips = sorted(name.values())
    S = [t for t in tips if str(lab.get(t, "")).startswith("S/")]
    children, parent, post, root, edge = s44.root_on_split(adj, name, S)
    tipidx = {nd: tips.index(nm) for nd, nm in name.items()}

    for pair in args.pairs.split(","):
        gene, trait = pair.split(":")
        if trait not in tr.columns:
            sys.exit("[!] trait {} not in trait.csv; columns: {}".format(trait, list(tr.columns)))
        hits = [c for c in GT.index if gene in str(c).split("~~~")]
        say("\n" + "=" * 100 + "\n{}  vs  {}".format(gene, trait))
        if not hits:
            say("  [!] no Panaroo cluster contains the token {}".format(gene))
            continue
        for c in hits:
            say("  gene_table: " + "; ".join("{}={}".format(k, v) for k, v in GT.loc[c].items()))
        c = hits[0] if len(hits) == 1 else max(hits, key=lambda h: GT.at[h, "len_bp"] if "len_bp" in GT.columns else 0)
        if c not in B.index:
            say("  [!] {} not in the read matrices (removed as short/redundant: see gene_table)".format(c))
            continue
        say("  cluster used: {}".format(c))
        br = B.loc[c, G]
        t3 = T3.loc[c, G]
        cls = pd.cut(br, [-0.1, 20, 80, 101], labels=["absent(<20)", "partial(20-80)", "full(>=80)"], right=False)
        say("  read breadth classes, all genomes with reads: " + ", ".join(
            "{} {}".format(k, int(v)) for k, v in cls.value_counts().sort_index().items()))
        tv = pd.to_numeric(tr[trait].reindex(G).replace("NA", np.nan), errors="coerce")
        say("  breadth by trait within C/non-B (median; share full / partial / absent):")
        for val, nm_ in ((1, "trait=1"), (0, "trait=0")):
            m = [g for g in cn if tv[g] == val]
            cc = cls[m].value_counts(normalize=True)
            say("    {:<8} n {:>4}  median {:6.1f}  full {:.3f}  partial {:.3f}  absent {:.3f}".format(
                nm_, len(m), br[m].median(), cc.get("full(>=80)", 0), cc.get("partial(20-80)", 0),
                cc.get("absent(<20)", 0)))
        if U is not None and c in U.index:
            u = U.loc[c]
            say("  union80 (assembly) frequency, C/non-B: trait=1 {:.3f} | trait=0 {:.3f}".format(
                u[[g for g in cn if tv.get(g) == 1 and g in u.index]].mean(),
                u[[g for g in cn if tv.get(g) == 0 and g in u.index]].mean()))

        say("  association on read calls (NA excluded), C/non-B:")
        say("    {:<26} {:>6} {:>6} {:>9} {:>9} {:>10}".format("subset", "n_pos", "n_neg", "freq_pos", "freq_neg", "Fisher p"))

        def line(lbl, gs):
            x = t3[gs].dropna()
            y = tv[x.index].dropna()
            x = x[y.index]
            npos, nneg, fp, fn, p = fisher(x, y)
            say("    {:<26} {:>6} {:>6} {:>9.3f} {:>9.3f} {:>10.3g}".format(lbl, npos, nneg, fp, fn, p))
        line("all C/non-B", cn)
        line("MiSeq only", [g for g in cn if "miseq" in str(plat[g]).lower()])
        for s, n in study[cn].value_counts().items():
            gs = [g for g in cn if study[g] == s]
            k1, k0 = int((tv[gs] == 1).sum()), int((tv[gs] == 0).sum())
            if k1 >= args.min_side and k0 >= args.min_side:
                line("study " + s, gs)

        # tree: presence pattern (genomes without reads = missing)
        x = pd.Series(-1, index=tips, dtype=np.int8)
        known = t3.dropna()
        x[known.index.intersection(tips)] = known[known.index.intersection(tips)].astype(np.int8)
        absent = 0 if (x == 0).sum() <= (x == 1).sum() else 1
        Mx = x.values.reshape(-1, 1)
        fc = int(s44.fitch(children, post, tipidx, Mx)[0])
        lack = [t for t in tips if x[t] == absent]
        say("  tree: {} genomes {} the gene (minor state), Fitch changes {}".format(
            len(lack), "lack" if absent == 0 else "carry", fc))
        bit = {t: 1 << k for k, t in enumerate(tips)}
        ALL = (1 << len(tips)) - 1
        mask = {}
        for nd in post:
            kids = children.get(nd, [])
            mask[nd] = bit[name[nd]] if not kids else sum(mask[k] for k in kids)
        carm = sum(bit[t] for t in lack)
        knownm = sum(bit[t] for t in tips if x[t] != -1)
        best = min(((((m ^ carm) & knownm).bit_count(), m) for nd, m in mask.items() if nd != root),
                   key=lambda z: z[0])
        say("  closest clade: {} genomes, {} genomes differ from it".format(bin(best[1]).count("1"), best[0]))
        L = pd.DataFrame({"genome": lack, "region": reg.reindex(lack).values, "host": host.reindex(lack).values,
                          "study": [study.get(t, "assembly") for t in lack],
                          "lineage": [lab.get(t, "?") for t in lack],
                          "breadth": [br.get(t, np.nan) for t in lack]})
        for col in ("region", "host", "study", "lineage"):
            say("    minor state by {}: {}".format(col, "; ".join(
                "{}:{}".format(k, v) for k, v in L[col].value_counts().head(6).items())))
        P = pd.DataFrame({"breadth": br, "ternary": t3, "trait": tv, "region": reg.reindex(G).values,
                          "host": host.reindex(G).values, "study": study, "platform": plat,
                          "lineage": [lab.get(g, "?") for g in G]})
        P.to_csv(out / "{}.tsv".format(gene), sep="\t", index_label="genome")
    (out / "report.txt").write_text("\n".join(R) + "\n")
    say("\nwritten: {}/".format(out))


if __name__ == "__main__":
    main()
