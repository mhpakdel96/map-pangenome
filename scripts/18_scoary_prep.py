#!/usr/bin/env python3
"""
18_scoary_prep.py

Prepare one Scoary input set per (isolate subset x trait), so that every run
is exactly specified on disk and can be repeated. Each run gets its own:
  genes.csv  Panaroo gene_presence_absence_roary.csv restricted to the run's
             isolates, keeping only genes that VARY among them (genes present
             in all or none cannot associate and would only inflate the
             Bonferroni denominator)
  trait.csv  one trait; isolates with NA for it (unknown on that axis) are
             left out of the run instead of being counted as negatives
  tree.nwk   the IQ-TREE ML tree pruned to the run's isolates (unary nodes
             collapsed, branch lengths summed, support labels removed, root
             made bifurcating) for Scoary's pairwise comparisons

Subsets (why each exists)
  C_nonB         PRIMARY. The S/C split confounds sheep/Oceania with lineage
                 (98 of 104 sheep are S). Within C (non-Bison) that confound
                 is gone.
  all            all 1792, for comparison with the previous report only.
  C_nonB_MiSeq   region sensitivity: MiSeq is the only platform used in all
                 three regions; other platforms are region-specific (NovaSeq
                 6000 = Americas, NextSeq 500 = Europe, HiSeq 2500 = Americas),
                 so a region hit could be a sequencing/assembly batch effect.
  C_nonB_N50ge10k  assembly-quality sensitivity: drops genomes with N50 < 10 kb.
  study_<ACC>    within-study test, given with --study (e.g. PRJNA686527 holds
                 95% of Source_Environment and 52% of Source_Invasive, so the
                 source traits are tested inside that one study).

A run is written only if the trait has >= --min-n positives AND negatives in
that subset; skipped runs are listed with the reason.

Outputs: 06_pangwas/scoary_runs/<subset>__<trait>/{genes.csv,trait.csv,tree.nwk}
         06_pangwas/scoary_runs/manifest.tsv   (one line per run; array input)
         06_pangwas/scoary_runs/prep_log.txt

Usage
  python 18_scoary_prep.py --study PRJNA686527
"""

import argparse
import csv
import sys
from pathlib import Path

import pandas as pd

ROOT = Path("/cfs/earth/scratch/tabasran/MAP")


# --------------------------------------------------------------- Newick I/O
class Node(object):
    __slots__ = ("name", "length", "children")

    def __init__(self):
        self.name, self.length, self.children = "", 0.0, []


def read_newick(text):
    """Iterative parser (no recursion limit). Keeps tip names and lengths."""
    text = text.strip().rstrip(";")
    root = Node()
    stack, cur, i, n = [], root, 0, len(text)

    def read_label(i):
        j = i
        while j < n and text[j] not in ",():;":
            j += 1
        return text[i:j].strip(), j

    while i < n:
        c = text[i]
        if c == "(":
            child = Node()
            cur.children.append(child)
            stack.append(cur)
            cur = child
            i += 1
        elif c == ",":
            parent = stack[-1]
            child = Node()
            parent.children.append(child)
            cur = child
            i += 1
        elif c == ")":
            cur = stack.pop()
            i += 1
            lab, i = read_label(i)          # internal label = support; dropped
        elif c == ":":
            val, i = read_label(i + 1)
            cur.length = float(val) if val else 0.0
        else:
            lab, i = read_label(i)
            cur.name = lab.strip("'\"")
    return root


def prune(node, keep):
    """Return a pruned copy, or None if no kept tip below."""
    # iterative post-order
    out = {}
    stack = [(node, False)]
    while stack:
        nd, done = stack.pop()
        if not done:
            stack.append((nd, True))
            for ch in nd.children:
                stack.append((ch, False))
            continue
        if not nd.children:
            if nd.name in keep:
                m = Node()
                m.name, m.length = nd.name, nd.length
                out[id(nd)] = m
            else:
                out[id(nd)] = None
            continue
        kids = [out[id(ch)] for ch in nd.children if out[id(ch)] is not None]
        if not kids:
            out[id(nd)] = None
        elif len(kids) == 1:
            k = kids[0]
            k.length += nd.length
            out[id(nd)] = k
        else:
            m = Node()
            m.length, m.children = nd.length, kids
            out[id(nd)] = m
    return out[id(node)]


def binarise_root(root):
    if len(root.children) > 2:
        rest = Node()
        rest.children, rest.length = root.children[1:], 0.0
        root.children = [root.children[0], rest]
    return root


def write_newick(root):
    """Iterative writer; tip names and branch lengths only, no root length."""
    out, st = [], [(root, "open")]
    while st:
        x, kind = st.pop()
        if kind == "comma":
            out.append(",")
        elif kind == "close":
            out.append(")" if x is root else "):{:.10g}".format(x.length))
        elif x.children:
            out.append("(")
            st.append((x, "close"))
            for k in range(len(x.children) - 1, -1, -1):
                st.append((x.children[k], "open"))
                if k > 0:
                    st.append((None, "comma"))
        else:
            out.append("{}:{:.10g}".format(x.name, x.length))
    return "".join(out) + ";"


def tips(root):
    t, st = [], [root]
    while st:
        x = st.pop()
        if x.children:
            st.extend(x.children)
        else:
            t.append(x.name)
    return t


# ------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes", default=str(ROOT / "04_pangenome/panaroo/gene_presence_absence_roary.csv"))
    ap.add_argument("--traits", default=str(ROOT / "06_pangwas/trait.csv"))
    ap.add_argument("--tree", default=str(ROOT / "05_phylogeny/iqtree/tree_final.nwk"))
    ap.add_argument("--lineages", default=str(ROOT / "05_phylogeny/typing/lineage_assignments.tsv"))
    ap.add_argument("--qc", default=str(ROOT / "02_qc/genome_qc_master.tsv"))
    ap.add_argument("--ena", default=str(ROOT / "01_metadata/ena_raw.tsv"))
    ap.add_argument("--out", default=str(ROOT / "06_pangwas/scoary_runs"))
    ap.add_argument("--study", action="append", default=[],
                    help="add a within-study subset (repeatable)")
    ap.add_argument("--min-n", type=int, default=20)
    ap.add_argument("--n50-min", type=int, default=10000)
    args = ap.parse_args()

    out = Path(args.out)
    # logs/ (SLURM output of 19) may already exist; earlier runs may not.
    if out.exists() and ((out / "manifest.tsv").exists() or
                         any(d.is_dir() and "__" in d.name for d in out.iterdir())):
        sys.exit("[!] {} already holds runs — move it away first".format(out))
    out.mkdir(parents=True, exist_ok=True)
    log = []

    def say(s=""):
        print(s)
        log.append(s)

    # ---------------------------------------------------------- inputs
    genes = pd.read_csv(args.genes, dtype=str, keep_default_na=False)
    traits = pd.read_csv(args.traits, dtype=str, keep_default_na=False).set_index("Name")
    cohort = list(traits.index)
    iso_cols = [c for c in genes.columns if c in set(cohort)]
    ann_cols = list(genes.columns[: genes.columns.get_loc(iso_cols[0])])
    if set(iso_cols) != set(cohort) or len(ann_cols) + len(iso_cols) != genes.shape[1]:
        sys.exit("[!] gene matrix isolates do not match trait.csv, or isolate "
                 "columns are not one contiguous block")
    say("genes      : {} rows, {} annotation columns, {} isolates".format(
        len(genes), len(ann_cols), len(iso_cols)))
    pres = genes[iso_cols].apply(lambda c: c.str.strip() != "")

    tree = read_newick(open(args.tree).read())
    if set(tips(tree)) != set(cohort):
        sys.exit("[!] tree tips do not match the cohort")

    lin = pd.read_csv(args.lineages, sep="\t", dtype=str).set_index("genome").reindex(cohort)
    lab = lin["lineage"] + "/" + lin["sublineage"]
    qc = pd.read_csv(args.qc, sep="\t").set_index("genome").reindex(cohort)
    ena = pd.read_csv(args.ena, sep="\t", dtype=str).drop_duplicates("run_accession").set_index("run_accession")
    instr = pd.Series([ena.at[g, "instrument_model"] if g in ena.index else "assembly"
                       for g in cohort], index=cohort)
    study = pd.Series([ena.at[g, "study_accession"] if g in ena.index else "assembly"
                       for g in cohort], index=cohort)
    if lab.isna().any() or qc["n50"].isna().any():
        sys.exit("[!] lineage or N50 missing for some genomes")

    c_nonb = set(lab.index[lab == "C/non-B"])
    subsets = [
        ("C_nonB", c_nonb),
        ("all", set(cohort)),
        ("C_nonB_MiSeq", c_nonb & set(instr.index[instr == "Illumina MiSeq"])),
        ("C_nonB_N50ge{}k".format(args.n50_min // 1000),
         c_nonb & set(qc.index[qc["n50"] >= args.n50_min])),
    ]
    for s in args.study:
        subsets.append(("study_" + s, set(study.index[study == s])))
    say("\nsubsets:")
    for name, iso in subsets:
        say("  {:<20} {:>5} genomes".format(name, len(iso)))

    # ---------------------------------------------------------- runs
    man = []
    say("\n{:<20} {:<22} {:>5} {:>5} {:>5} {:>7}  status".format(
        "subset", "trait", "pos", "neg", "NA", "genes"))
    for sname, iso in subsets:
        for t in traits.columns:
            v = traits.loc[sorted(iso), t]
            known = v[v.isin(["0", "1"])]
            pos, neg = int((known == "1").sum()), int((known == "0").sum())
            na = len(v) - len(known)
            if min(pos, neg) < args.min_n:
                say("{:<20} {:<22} {:>5} {:>5} {:>5} {:>7}  skipped (<{} pos or neg)".format(
                    sname, t, pos, neg, na, "-", args.min_n))
                continue
            keep = list(known.index)
            p = pres.loc[:, keep]
            k = p.sum(axis=1)
            var = (k > 0) & (k < len(keep))
            d = out / "{}__{}".format(sname, t)
            d.mkdir()
            g = genes.loc[var, ann_cols + keep]
            g.to_csv(str(d / "genes.csv"), index=False, quoting=csv.QUOTE_MINIMAL)
            pd.DataFrame({"Name": keep, t: known.values}).to_csv(str(d / "trait.csv"), index=False)
            pt = binarise_root(prune(tree, set(keep)))
            if sorted(tips(pt)) != sorted(keep):
                sys.exit("[!] pruning failed for " + d.name)
            (d / "tree.nwk").write_text(write_newick(pt) + "\n")
            man.append((d.name, sname, t, pos, neg, na, int(var.sum())))
            say("{:<20} {:<22} {:>5} {:>5} {:>5} {:>7}  written".format(
                sname, t, pos, neg, na, int(var.sum())))

    with open(str(out / "manifest.tsv"), "w") as fh:
        fh.write("run\tsubset\ttrait\tpos\tneg\tNA\tgenes_tested\n")
        for r in man:
            fh.write("\t".join(str(x) for x in r) + "\n")
    say("\nruns written: {}  ->  {}/manifest.tsv".format(len(man), out))
    say("annotation columns kept: {} -> Scoary isolate start column (-s) = {}".format(
        len(ann_cols), len(ann_cols) + 1))
    (out / "prep_log.txt").write_text("\n".join(log) + "\n")


if __name__ == "__main__":
    main()
