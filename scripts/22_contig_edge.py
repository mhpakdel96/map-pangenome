#!/usr/bin/env python3
"""
22_contig_edge.py

Mechanism test for the Scoary batch effect: are the non-robust hits genes
that sit at CONTIG ENDS (i.e. truncated by assembly breaks) or that Panaroo
had to RE-FIND because the annotation missed them?

Background (21_batch_check.py, N50 by study). Genes depleted in the Americas
are carried mainly by low-N50 assemblies (median ~13 kb; PRJNA841672, NextSeq
500, N50 14.5 kb) and Americas-enriched genes are missing from them. If the
hits are assembly artefacts, their copies should lie at contig edges or be
re-found far more often than ordinary core or accessory genes.

Categories compared (gene x genome occurrences)
  hit_nonrobust_enriched / hit_nonrobust_depleted / hit_robust
      from 06_pangwas/scoary_summary/batch_check.tsv (robust = supported
      inside every eligible platform and not batch-suspect)
  baseline_core       random --n-base genes present in >= 99% of genomes
  baseline_accessory  random --n-base non-hit genes present in 5-95%
Per occurrence:
  edge      CDS starts or ends within --edge-bp of its contig end
  refound   Panaroo ID with 'refound' (no CDS in the Prokka GFF)
  short_contig  contig shorter than --short-contig bp
Also printed: pangenome-graph position ('Genome Fragment', 'Order within
Fragment') of the robust hits, to see whether they form one locus.

Usage (reads ~1800 GFF headers; run as a job):
  python 22_contig_edge.py
"""

import argparse
import random
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd

ROOT = Path("/cfs/earth/scratch/tabasran/MAP")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes", default=str(ROOT / "04_pangenome/panaroo/gene_presence_absence_roary.csv"))
    ap.add_argument("--gff", default=str(ROOT / "03_genomes/gff"))
    ap.add_argument("--batch", default=str(ROOT / "06_pangwas/scoary_summary/batch_check.tsv"))
    ap.add_argument("--ena", default=str(ROOT / "01_metadata/ena_raw.tsv"))
    ap.add_argument("--out", default=str(ROOT / "06_pangwas/scoary_summary"))
    ap.add_argument("--edge-bp", type=int, default=100)
    ap.add_argument("--short-contig", type=int, default=5000)
    ap.add_argument("--n-base", type=int, default=300)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    gm = pd.read_csv(args.genes, dtype=str, keep_default_na=False)
    iso = list(gm.columns[14:])
    gm = gm.set_index("Gene")
    pres = (gm[iso] != "")
    freq = pres.mean(axis=1)

    b = pd.read_csv(args.batch, sep="\t")
    b["robust"] = (b["within_platform_supported"] == "yes") & (b["batch_suspect"] != "yes")
    cat = {}
    for _, r in b.iterrows():
        c = "hit_robust" if r["robust"] else "hit_nonrobust_" + r["direction"]
        # a gene that is a hit for several traits: robust wins, else first seen
        if cat.get(r["Gene"]) != "hit_robust":
            cat[r["Gene"]] = c
    hits = set(cat)
    rnd = random.Random(args.seed)
    core = sorted(g for g in freq.index[freq >= 0.99] if g not in hits)
    acc = sorted(g for g in freq.index[(freq >= 0.05) & (freq <= 0.95)] if g not in hits)
    for g in rnd.sample(core, min(args.n_base, len(core))):
        cat[g] = "baseline_core"
    for g in rnd.sample(acc, min(args.n_base, len(acc))):
        cat[g] = "baseline_accessory"
    cat = {g: c for g, c in cat.items() if g in gm.index}

    # per genome: which IDs to look up, and their category
    want = defaultdict(dict)
    refound = defaultdict(lambda: [0, 0])          # category -> [refound, total]
    occ_gene = {}
    for g, c in cat.items():
        for genome in iso:
            v = gm.at[g, genome]
            if not v:
                continue
            for i in re.split(r"[;\t]", v):
                i = i.strip()
                if not i:
                    continue
                refound[c][1] += 1
                if "refound" in i:
                    refound[c][0] += 1
                else:
                    want[genome][i] = (c, g)

    ex = [v for v in gm.loc[list(cat)[:5], iso[:3]].values.ravel() if v][:3]
    print("example ID cells: {}".format(ex))
    ena = pd.read_csv(args.ena, sep="\t", dtype=str).drop_duplicates("run_accession").set_index("run_accession")
    study = {x: (ena.at[x, "study_accession"] if x in ena.index else "assembly") for x in iso}

    stats = defaultdict(lambda: [0, 0, 0])           # category -> [edge, short, found]
    by_study = defaultdict(lambda: [0, 0])            # (category, study) -> [edge, found]
    missing_ids = 0
    for k, genome in enumerate(iso, 1):
        ids = want.get(genome)
        if not ids:
            continue
        clen, found = {}, 0
        with open(str(Path(args.gff) / (genome + ".gff"))) as fh:
            for line in fh:
                if line.startswith("##FASTA"):
                    break
                if line.startswith("##sequence-region"):
                    p = line.split()
                    clen[p[1]] = int(p[3])
                    continue
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9 or f[2] != "CDS":
                    continue
                m = re.search(r"(?:^|;)ID=([^;]+)", f[8])
                if not m or m.group(1) not in ids:
                    continue
                c, _ = ids[m.group(1)]
                L = clen.get(f[0])
                if L is None:
                    continue
                s, e = int(f[3]), int(f[4])
                edge = s <= args.edge_bp or e >= L - args.edge_bp
                st = stats[c]
                st[0] += edge
                st[1] += L < args.short_contig
                st[2] += 1
                bs = by_study[(c, study[genome])]
                bs[0] += edge
                bs[1] += 1
                found += 1
        missing_ids += len(ids) - found
        if k % 300 == 0:
            print("  GFFs read: {} / {}".format(k, len(iso)))

    order = ["baseline_core", "baseline_accessory", "hit_robust",
             "hit_nonrobust_enriched", "hit_nonrobust_depleted"]
    lines = ["edge = CDS within {} bp of a contig end; short contig < {} bp".format(
        args.edge_bp, args.short_contig), ""]
    lines.append("{:<24} {:>6} {:>10} {:>9} {:>11} {:>14}".format(
        "category", "genes", "occurr.", "edge %", "refound %", "short-contig %"))
    for c in order:
        ng = sum(1 for v in cat.values() if v == c)
        e, sh, n = stats[c]
        rf, tot = refound[c]
        lines.append("{:<24} {:>6} {:>10} {:>8.1f}% {:>10.1f}% {:>13.1f}%".format(
            c, ng, tot, 100.0 * e / n if n else float("nan"),
            100.0 * rf / tot if tot else float("nan"),
            100.0 * sh / n if n else float("nan")))
    lines.append("\nIDs not found in their GFF (should be ~0): {}".format(missing_ids))

    lines.append("\nedge % by study (studies with >= 200 occurrences in the category)")
    for c in order:
        rows = sorted(((s, v[0], v[1]) for (cc, s), v in by_study.items()
                       if cc == c and v[1] >= 200), key=lambda x: -x[1] / x[2])
        if rows:
            lines.append("  " + c + ": " + "  ".join(
                "{}={:.1f}%".format(s, 100.0 * e / n) for s, e, n in rows))

    rob = b[b["robust"]][["trait", "Gene", "direction", "delta_all"]].copy()
    if len(rob):
        gf = gm.reindex(rob["Gene"])
        rob["genome_fragment"] = gf["Genome Fragment"].values
        rob["order_in_fragment"] = gf["Order within Fragment"].values
        rob["freq_all"] = freq.reindex(rob["Gene"]).round(3).values
        lines.append("\n--- robust hits: position in the pangenome graph ---")
        lines.append(rob.sort_values(["genome_fragment", "order_in_fragment"]).to_string(index=False))

    txt = "\n".join(lines)
    Path(args.out).mkdir(parents=True, exist_ok=True)
    (Path(args.out) / "contig_edge_test.txt").write_text(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()
