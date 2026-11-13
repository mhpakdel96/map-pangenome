#!/usr/bin/env python3
"""
53_figures.py

Main and supplemental figures of the MAP manuscript (Fig. 1 = workflow
diagram, drawn separately in LaTeX). All panels are drawn from the cluster
outputs; values that exist only in logs of earlier steps are written below as
constants with the script that produced them.

  Fig2  circular cladogram of 1,792 genomes; rings: lineage/clade, region, host
  Fig3  (A) C-lineage cladogram, clade D marked, ring = region
        (B) regional composition of clade D and the basal grade
  Fig4  (A) changes per SNP allele before/after read verification
        (B) explanations of the 21 SNPs that first appeared homoplastic
        (C) genomes carrying icl T3M (4421235 T) on the cladogram
  Fig5  (A) spurious genes in the within-study negative control by method
        (B) root-to-tip distance vs date, C/non-B (52_rtt_data.R)
        (C) permuted vs observed rate: free / within region / within study
        (D) balanced accuracy of origin prediction by validation scheme, model
  FigS2 gene accumulation curves by presence method
  FigS3 batch artefacts of the Scoary hits (values from 21 and 22)
  FigS4 root-to-tip within single studies (52_rtt_data.R)

Size: 17.4 cm wide at most; PDF (vector) and PNG (300 dpi).
  python 53_figures.py [--only Fig2,Fig5]
Environment: ml_env + matplotlib (pip install matplotlib if missing).
Imports 44_homoplasy_site.py.
"""

import argparse
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Patch

ROOT = Path("/cfs/earth/scratch/tabasran/MAP")
W = 6.85                                  # 17.4 cm in inches
plt.rcParams.update({"font.family": "sans-serif", "font.size": 7, "axes.titlesize": 8,
                     "axes.labelsize": 7, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
                     "legend.fontsize": 6.5, "pdf.fonttype": 42, "axes.spines.top": False,
                     "axes.spines.right": False})
# Okabe-Ito palette
OI = {"orange": "#E69F00", "sky": "#56B4E9", "green": "#009E73", "yellow": "#F0E442",
      "blue": "#0072B2", "red": "#D55E00", "purple": "#CC79A7", "grey": "#999999", "black": "#000000"}
COL_GROUP = {"S/Type I": OI["orange"], "S/Type III": OI["yellow"], "Bison": OI["purple"],
             "C basal grade": OI["sky"], "C clade D": OI["blue"]}
COL_REGION = {"Americas": OI["red"], "Europe": OI["green"], "Oceania": OI["orange"], "other/unknown": "#DDDDDD"}
COL_HOST = {"cattle": "#555555", "sheep": OI["orange"], "other ruminants": OI["purple"], "other/unknown": "#DDDDDD"}

# ---- values recorded in earlier logs (not stored as tables) -----------------
# within-study negative control, significant genes (26_call_validation, 35, 36)
FP = pd.DataFrame({"method": ["Annotation", "DNA (union80)", "Reads, two states", "Reads, three states"],
                   "N50 split": [278, 115, 212, 23], "Depth split": [np.nan, 63, 260, 34]})
# explanation of the 21 SNPs first flagged as homoplastic (45_discordant_tips)
HOMOPLASY_EXPL = [("Mixed-strain samples", 16), ("Single terminal change", 3),
                  ("Assembly without reads", 1), ("Candidate parallel mutation", 1)]
# 21_batch_check: median |frequency difference| of hit genes, all genomes vs within one platform
BATCH = pd.DataFrame({"trait": ["Americas", "Europe", "Oceania", "Environment", "Excretory", "Invasive", "Cattle"],
                      "all genomes": [0.07, 0.09, 0.14, 0.38, 0.11, 0.13, 0.23],
                      "within platform": [0.01, 0.00, 0.09, 0.04, 0.05, 0.07, 0.04]})
# 22_contig_edge: % of gene occurrences within 100 bp of a contig end
EDGE = pd.DataFrame({"set": ["Core genes", "Accessory genes", "Hits, non-robust (enriched)",
                             "Hits, non-robust (depleted)", "Hits, robust"],
                     "pct": [2.6, 8.2, 15.8, 16.0, 56.9]})


def load44(path):
    sp = importlib.util.spec_from_file_location("s44", path)
    m = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(m)
    return m


def save(fig, out, name):
    fig.savefig(out / (name + ".pdf"), bbox_inches="tight")
    fig.savefig(out / (name + ".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  written", name)


def panel(ax, letter):
    ax.text(-0.08, 1.04, letter, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom")


# ---------------------------------------------------------------- tree tools
class Tree:
    def __init__(self, s44, tree_path, S):
        adj, name, sup = s44.parse_newick(open(tree_path).read())
        self.name = name
        self.children, parent, post, self.root, _ = s44.root_on_split(adj, name, S)
        self.tips = sorted(name.values())

    def subtree(self, keep):
        """children dict restricted to tips in keep, unary nodes collapsed."""
        keep = set(keep)
        ch = {}

        def rec(x):
            kids = self.children.get(x, [])
            if not kids:
                return x if self.name.get(x) in keep else None
            k2 = [y for y in (rec(k) for k in kids) if y is not None]
            if not k2:
                return None
            if len(k2) == 1:
                return k2[0]
            ch[x] = k2
            return x
        import sys
        sys.setrecursionlimit(100000)
        r = rec(self.root)
        return ch, r

    def layout(self, children=None, root=None):
        """circular cladogram: tip order by ladderised DFS, node height = max
        edges to a tip; returns angle, radius per node and tip order."""
        ch = children if children is not None else self.children
        r0 = root if root is not None else self.root
        size, height, order = {}, {}, []
        stack, post = [r0], []
        while stack:
            x = stack.pop()
            post.append(x)
            stack.extend(ch.get(x, []))
        for x in reversed(post):
            kids = ch.get(x, [])
            size[x] = 1 if not kids else sum(size[k] for k in kids)
            height[x] = 0 if not kids else 1 + max(height[k] for k in kids)
        stack = [r0]
        while stack:
            x = stack.pop()
            kids = sorted(ch.get(x, []), key=lambda k: size[k])
            if not kids:
                order.append(x)
            stack.extend(reversed(kids))
        n, H = len(order), height[r0]
        span = 2 * np.pi * 0.985
        ang = {t: span * i / n for i, t in enumerate(order)}
        for x in reversed(post):
            kids = ch.get(x, [])
            if kids:
                a = [ang[k] for k in kids]
                ang[x] = (min(a) + max(a)) / 2
        rad = {x: 1 - height[x] / H for x in post}
        return ang, rad, order, ch, r0

    def draw(self, ax, lay, lw=0.25, color="#333333"):
        ang, rad, order, ch, r0 = lay
        segs = []
        for x, kids in ch.items():
            if not kids:
                continue
            a = [ang[k] for k in kids]
            th = np.linspace(min(a), max(a), max(2, int((max(a) - min(a)) * 60)))
            segs.append(np.column_stack([th, np.full_like(th, rad[x])]))
            for k in kids:
                segs.append(np.array([[ang[k], rad[x]], [ang[k], rad[k]]]))
        ax.add_collection(LineCollection(segs, colors=color, linewidths=lw))
        ax.set_ylim(0, 1.32)
        ax.set_axis_off()

    def ring(self, ax, lay, colors, r0, h):
        ang, rad, order, ch, root = lay
        w = 2 * np.pi * 0.985 / len(order)
        th = [ang[t] for t in order]
        ax.bar(th, [h] * len(order), width=w * 1.02, bottom=r0,
               color=[colors.get(self.name[t], "#DDDDDD") for t in order], linewidth=0)


# ---------------------------------------------------------------- metadata
def metadata(root):
    lin = pd.read_csv(root / "05_phylogeny/typing/lineage_assignments.tsv", sep="\t", dtype=str).set_index("genome")
    sp = pd.read_csv(root / "05_phylogeny/typing/split_assignments.tsv", sep="\t", dtype=str).set_index("genome")
    tr = pd.read_csv(root / "06_pangwas/trait.csv", dtype=str).set_index("Name")
    g = pd.Series("C basal grade", index=lin.index)
    lab = lin["lineage"] + "/" + lin["sublineage"]
    g[lab == "S/non-III"] = "S/Type I"
    g[lab == "S/III"] = "S/Type III"
    g[lin["sublineage"] == "B"] = "Bison"
    g[sp["tree_side"].reindex(lin.index).fillna("").str.startswith("D ")] = "C clade D"
    reg = pd.Series("other/unknown", index=tr.index)
    for c in ("Region_Americas", "Region_Europe", "Region_Oceania"):
        if c in tr.columns:
            reg[tr[c] == "1"] = c.replace("Region_", "")
    host = pd.Series("other/unknown", index=tr.index)
    hc = {c.lower(): c for c in tr.columns if c.lower().startswith("host")}
    for key, val in (("hostother_ruminants", "other ruminants"), ("hostsheep", "sheep"), ("hostcattle", "cattle")):
        if key in hc:
            host[tr[hc[key]] == "1"] = val
    return g, reg, host


def legend(ax, colors, title, loc, anchor):
    ax.legend(handles=[Patch(color=c, label=k) for k, c in colors.items()], title=title, loc=loc,
              bbox_to_anchor=anchor, frameon=False, title_fontsize=7, handlelength=1, handleheight=1)


# ---------------------------------------------------------------- figures
def fig2(T, g, reg, host, out):
    fig = plt.figure(figsize=(W, W * 0.82))
    ax = fig.add_axes([0.0, 0.0, 0.78, 1.0], projection="polar")
    lay = T.layout()
    T.draw(ax, lay, lw=0.2)
    T.ring(ax, lay, {t: COL_GROUP[g.get(t, "C basal grade")] for t in T.tips}, 1.03, 0.07)
    T.ring(ax, lay, {t: COL_REGION[reg.get(t, "other/unknown")] for t in T.tips}, 1.12, 0.07)
    T.ring(ax, lay, {t: COL_HOST[host.get(t, "other/unknown")] for t in T.tips}, 1.21, 0.07)
    lg = fig.add_axes([0.78, 0.05, 0.22, 0.9])
    lg.set_axis_off()
    for i, (cols, title) in enumerate(((COL_GROUP, "Lineage (inner ring)"), (COL_REGION, "Region (middle ring)"),
                                       (COL_HOST, "Host (outer ring)"))):
        lg.legend(handles=[Patch(color=c, label=k) for k, c in cols.items()], title=title, loc="upper left",
                  bbox_to_anchor=(0, 1 - i * 0.3), frameon=False, title_fontsize=7, handlelength=1)
        lg.add_artist(lg.get_legend())
    save(fig, out, "Fig2")


def fig3(T, g, reg, out):
    keep = [t for t in T.tips if g.get(t) in ("C basal grade", "C clade D")]
    ch, r0 = T.subtree(keep)
    lay = T.layout(ch, r0)
    fig = plt.figure(figsize=(W, W * 0.55))
    ax = fig.add_axes([0.0, 0.0, 0.55, 1.0], projection="polar")
    ang, rad, order, _, _ = lay
    T.draw(ax, lay, lw=0.2)
    T.ring(ax, lay, {t: COL_REGION[reg.get(t, "other/unknown")] for t in keep}, 1.03, 0.08)
    T.ring(ax, lay, {t: COL_GROUP[g[t]] for t in keep}, 1.13, 0.05)
    ax.text(0.5, 1.0, "A", transform=ax.transAxes, fontsize=9, fontweight="bold")
    ax2 = fig.add_axes([0.66, 0.2, 0.32, 0.6])
    tab = pd.crosstab(g[keep], reg.reindex(keep).fillna("other/unknown"))
    tab = tab.reindex(index=["C clade D", "C basal grade"],
                      columns=[c for c in COL_REGION if c in tab.columns], fill_value=0)
    left = np.zeros(len(tab))
    for c in tab.columns:
        ax2.barh(range(len(tab)), tab[c].values, left=left, color=COL_REGION[c], label=c, height=0.6)
        left += tab[c].values
    ax2.set_yticks(range(len(tab)))
    ax2.set_yticklabels(["Clade D\n(n = {:,})".format(int(tab.loc["C clade D"].sum())),
                         "Basal grade\n(n = {:,})".format(int(tab.loc["C basal grade"].sum()))])
    ax2.invert_yaxis()
    ax2.set_xlabel("Genomes")
    ax2.legend(frameon=False, loc="lower right", fontsize=6)
    panel(ax2, "B")
    save(fig, out, "Fig3")
    print("    clade D x region:\n" + tab.to_string())


def fig4(T, g, root, out):
    A = pd.read_csv(root / "06_pangwas/convergence_v2/site_events.tsv", sep="\t")
    fig = plt.figure(figsize=(W, W * 0.45))
    ax = fig.add_axes([0.06, 0.15, 0.27, 0.75])
    bins = [1, 2, 3, 4, 5, 6, 8, 11]
    lab = ["1", "2", "3", "4", "5", "6-7", "8-10"]

    def counts(col, carr):
        x = A.loc[A[carr] > 0, col]
        return [int(((x >= lo) & (x < hi)).sum()) for lo, hi in zip(bins[:-1], bins[1:])]
    xx = np.arange(len(lab))
    ax.bar(xx - 0.2, counts("fitch_raw", "carriers_raw"), 0.4, color=OI["grey"], label="Assembly calls")
    ax.bar(xx + 0.2, counts("fitch", "carriers"), 0.4, color=OI["blue"], label="Read-verified calls")
    ax.set_yscale("log")
    ax.set_xticks(xx)
    ax.set_xticklabels(lab)
    ax.set_xlabel("Changes on the tree per SNP allele")
    ax.set_ylabel("SNP alleles")
    ax.legend(frameon=False)
    panel(ax, "A")
    ax = fig.add_axes([0.5, 0.15, 0.16, 0.75])
    ax.barh(range(len(HOMOPLASY_EXPL)), [v for _, v in HOMOPLASY_EXPL],
            color=[OI["red"], OI["sky"], OI["grey"], OI["green"]])
    ax.set_yticks(range(len(HOMOPLASY_EXPL)))
    ax.set_yticklabels([k for k, _ in HOMOPLASY_EXPL])
    ax.invert_yaxis()
    ax.set_xlabel("SNPs (of 21)")
    panel(ax, "B")
    ct = pd.read_csv(root / "05_phylogeny/core_masked/core.tab", sep="\t", dtype=str)
    row = ct[ct["POS"].astype(int) == 4421235].iloc[0]
    car = {t for t in T.tips if str(row.get(t, "")).upper() == "T"}
    ax = fig.add_axes([0.68, 0.0, 0.32, 1.0], projection="polar")
    lay = T.layout()
    T.draw(ax, lay, lw=0.15, color="#999999")
    T.ring(ax, lay, {t: COL_GROUP[g.get(t, "C basal grade")] for t in T.tips}, 1.03, 0.06)
    T.ring(ax, lay, {t: (OI["red"] if t in car else "#FFFFFF") for t in T.tips}, 1.11, 0.12)
    ax.text(0.0, 1.0, "C", transform=ax.transAxes, fontsize=9, fontweight="bold")
    ax.text(0.5, -0.02, "outer ring: genomes with icl 4421235 T (Thr3Met); n = {}".format(len(car)),
            transform=ax.transAxes, ha="center", fontsize=6)
    save(fig, out, "Fig4")


def fig5(root, out):
    fig = plt.figure(figsize=(W, W * 0.62))
    ax = fig.add_axes([0.07, 0.58, 0.38, 0.36])
    xx = np.arange(len(FP))
    ax.bar(xx - 0.2, FP["N50 split"], 0.4, color=OI["blue"], label="Split at study median N50")
    ax.bar(xx + 0.2, FP["Depth split"].fillna(0), 0.4, color=OI["sky"], label="Split at study median depth")
    ax.set_xticks(xx)
    ax.set_xticklabels(FP["method"], rotation=20, ha="right")
    ax.set_ylabel("Spurious significant genes")
    ax.legend(frameon=False)
    panel(ax, "A")
    R = pd.read_csv(root / "09_dating/fig/rtt_C_nonB.tsv", sep="\t")
    ax = fig.add_axes([0.55, 0.58, 0.32, 0.36])
    for r, c in COL_REGION.items():
        s = R[R["region"].replace({"unknown": "other/unknown"}) == r]
        ax.scatter(s["date"], s["rtt"], s=2, color=c, label=r, linewidths=0)
    b = np.polyfit(R["date"], R["rtt"], 1)
    xs = np.linspace(R["date"].min(), R["date"].max(), 10)
    ax.plot(xs, np.polyval(b, xs), color="black", lw=0.8)
    r2 = np.corrcoef(R["date"], R["rtt"])[0, 1] ** 2
    ax.text(0.02, 0.95, "R$^2$ = {:.3f}".format(r2), transform=ax.transAxes, va="top")
    ax.set_xlabel("Sampling date")
    ax.set_ylabel("Root-to-tip distance")
    ax.legend(frameon=False, markerscale=3, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    panel(ax, "B")
    P = pd.read_csv(root / "09_dating/fig/perm_C_nonB.tsv", sep="\t")
    obs = float(P.loc[P["scheme"] == "observed", "rate"].iloc[0])
    ax = fig.add_axes([0.07, 0.08, 0.38, 0.36])
    sch = [("free", "All dates\npermuted"), ("within_region", "Within\nregions"), ("within_study", "Within\nstudies")]
    data = [P.loc[P["scheme"] == s, "rate"].values for s, _ in sch]
    ax.boxplot(data, widths=0.5, showfliers=False, medianprops={"color": "black"})
    ax.axhline(obs, color=OI["red"], lw=0.8, ls="--")
    ax.text(3.45, obs, "observed", color=OI["red"], va="center", fontsize=6)
    ax.set_xticks([1, 2, 3])
    ax.set_xticklabels([l for _, l in sch])
    ax.set_ylabel("Root-to-tip slope")
    panel(ax, "C")
    M = pd.read_csv(root / "07_ml/C_nonB/metrics.tsv", sep="\t")
    M = M[M["label"] == "region"]
    ax = fig.add_axes([0.57, 0.08, 0.4, 0.36])
    schemes = [("random_cv", "Random\ncross-validation"), ("leave_study_out", "Leave one\nstudy out")]
    models = [("tree_1NN", "Nearest relative (tree)"), ("snp_1NN", "Nearest relative (SNPs)"),
              ("logreg_L1", "L1 logistic regression"), ("rf", "Random forest")]
    cols = [OI["grey"], "#BBBBBB", OI["blue"], OI["green"]]
    for i, (m, ml) in enumerate(models):
        v = [float(M.loc[(M["scheme"] == s) & (M["model"] == m), "balanced_accuracy"].iloc[0]) for s, _ in schemes]
        ax.bar(np.arange(len(schemes)) + (i - 1.5) * 0.18, v, 0.18, color=cols[i], label=ml)
    ax.axhline(1 / 3, color="black", lw=0.6, ls=":")
    ax.text(1.45, 1 / 3 + 0.02, "chance", fontsize=6)
    ax.set_xticks(range(len(schemes)))
    ax.set_xticklabels([l for _, l in schemes])
    ax.set_ylim(0, 1.3)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_ylabel("Balanced accuracy (region)")
    ax.legend(frameon=False, fontsize=5.5, loc="upper center", ncol=2)
    panel(ax, "D")
    save(fig, out, "Fig5")


def figS2(root, out, nperm=20):
    rng = np.random.default_rng(1)
    srcs = [("Annotation (Panaroo)", root / "04_pangenome/panaroo/gene_presence_absence.Rtab", "\t", OI["grey"]),
            ("DNA (union80)", root / "04_pangenome/dna_presence/presence_final.Rtab", "\t", OI["green"]),
            ("Reads, two states", root / "08_reads/summary/presence_reads.Rtab", "\t", OI["sky"]),
            ("Reads, three states (undetermined = absent)", root / "08_reads/summary/presence_reads_ternary.tsv", "\t", OI["blue"])]
    fig, ax = plt.subplots(figsize=(W * 0.55, W * 0.4))
    for lab, p, sep, c in srcs:
        if not p.exists():
            print("  [i] missing", p)
            continue
        X = pd.read_csv(p, sep=sep, index_col=0).fillna(0).values > 0.5
        n = X.shape[1]
        curves = []
        for _ in range(nperm):
            o = rng.permutation(n)
            Y = X[:, o]
            first = np.where(Y.any(axis=1), Y.argmax(axis=1), n + 1)
            curves.append(np.searchsorted(np.sort(first), np.arange(1, n + 1), side="left"))
        C = np.array(curves)
        xs = np.arange(1, n + 1)
        ax.plot(xs, C.mean(axis=0), color=c, lw=1, label=lab)
        ax.fill_between(xs, C.min(axis=0), C.max(axis=0), color=c, alpha=0.2, lw=0)
    ax.set_xlabel("Genomes sampled")
    ax.set_ylabel("Gene clusters (pangenome size)")
    ax.legend(frameon=False, fontsize=6)
    save(fig, out, "FigS2")


def figS3(out):
    fig, axs = plt.subplots(1, 2, figsize=(W, W * 0.33))
    ax = axs[0]
    xx = np.arange(len(BATCH))
    ax.bar(xx - 0.2, BATCH["all genomes"], 0.4, color=OI["grey"], label="All C genomes")
    ax.bar(xx + 0.2, BATCH["within platform"], 0.4, color=OI["blue"], label="Within one platform")
    ax.set_xticks(xx)
    ax.set_xticklabels(BATCH["trait"], rotation=30, ha="right")
    ax.set_ylabel("Median |frequency difference|\nof associated genes")
    ax.legend(frameon=False)
    panel(ax, "A")
    ax = axs[1]
    ax.barh(range(len(EDGE)), EDGE["pct"], color=[OI["grey"], "#BBBBBB", OI["sky"], OI["sky"], OI["red"]])
    ax.set_yticks(range(len(EDGE)))
    ax.set_yticklabels(EDGE["set"])
    ax.invert_yaxis()
    ax.set_xlabel("Gene copies within 100 bp of a contig end (%)")
    panel(ax, "B")
    fig.tight_layout()
    save(fig, out, "FigS3")


def figS4(root, out):
    B = pd.read_csv(root / "09_dating/fig/rtt_by_study.tsv", sep="\t")
    studies = list(B["study"].unique())
    fig, axs = plt.subplots(1, len(studies), figsize=(W, W * 0.3), squeeze=False)
    for ax, s in zip(axs[0], studies):
        d = B[B["study"] == s]
        ax.scatter(d["date"], d["rtt"], s=3, color=OI["blue"], linewidths=0)
        b = np.polyfit(d["date"], d["rtt"], 1)
        xs = np.linspace(d["date"].min(), d["date"].max(), 10)
        ax.plot(xs, np.polyval(b, xs), color="black", lw=0.8)
        r2 = np.corrcoef(d["date"], d["rtt"])[0, 1] ** 2
        ax.set_title("{} (n = {})".format(s, len(d)))
        ax.text(0.03, 0.95, "slope {:.2f}\nR$^2$ = {:.3f}".format(b[0], r2), transform=ax.transAxes, va="top", fontsize=6)
        ax.set_xlabel("Sampling date")
    axs[0][0].set_ylabel("Root-to-tip distance")
    fig.tight_layout()
    save(fig, out, "FigS4")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=str(ROOT))
    ap.add_argument("--only", default="")
    ap.add_argument("--s44", default=str(Path(__file__).resolve().parent / "44_homoplasy_site.py"))
    args = ap.parse_args()
    root = Path(args.root)
    out = root / "10_figures"
    out.mkdir(parents=True, exist_ok=True)
    only = set(x.strip() for x in args.only.split(",") if x.strip())
    want = lambda f: not only or f in only
    s44 = load44(args.s44)
    g, reg, host = metadata(root)
    T = None
    if any(want(f) for f in ("Fig2", "Fig3", "Fig4")):
        S = [t for t in g.index if g[t].startswith("S/")]
        T = Tree(s44, root / "05_phylogeny/iqtree/tree_final.nwk", S)
    if want("Fig2"):
        fig2(T, g, reg, host, out)
    if want("Fig3"):
        fig3(T, g, reg, out)
    if want("Fig4"):
        fig4(T, g, root, out)
    if want("Fig5"):
        fig5(root, out)
    if want("FigS2"):
        figS2(root, out)
    if want("FigS3"):
        figS3(out)
    if want("FigS4"):
        figS4(root, out)
    print("figures in", out)


if __name__ == "__main__":
    main()
