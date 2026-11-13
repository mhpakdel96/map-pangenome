#!/usr/bin/env python3
"""
23_sequence_presence.py

Decisive test: is a Scoary hit a difference in SEQUENCE content, or only in
how the gene was assembled/annotated?

22_contig_edge.py showed that hit genes sit at contig ends far more often
than core genes (non-robust hits ~16%, the 12 'robust' hits 57%, core 2.6%).
A gene that ends at a contig break is often mis-called (truncated, split
into two clusters, or not called at all), so Panaroo's presence/absence can
differ between genomes that carry the same DNA. This script calls presence
from the DNA alone, independent of gene calling and contig breaks:

  each gene's representative sequence (pan_genome_reference.fa) is aligned
  with minimap2 to every assembly; aligned query intervals with identity
  >= --min-id are merged ACROSS contigs, so a gene split by a break still
  counts as covered. Present = covered fraction >= --min-cov.

Genes tested: every primary hit (batch_check.tsv) plus --n-base random core
and accessory genes as controls (same seed as 22).

Reported
  per category: how often Panaroo and sequence disagree
      ann0_seq1  'absent' in Panaroo but the DNA is there (annotation miss)
      ann1_seq0  'present' in Panaroo but < min-cov of the DNA is there
  per hit: delta (trait+ minus trait-) from Panaroo vs from sequence, over
      all isolates of the C_nonB run and inside MiSeq; 'retained' if the
      sequence delta keeps the sign and >= 50% of the Panaroo delta.

Outputs (06_pangwas/scoary_summary/)
  seq_coverage.tsv.gz     gene x genome covered fraction
  seq_presence_hits.tsv   per hit: deltas and verdict
  seq_presence_summary.txt
"""

import argparse
import random
import subprocess
import sys
from collections import defaultdict
from multiprocessing import Pool
from pathlib import Path

import pandas as pd

ROOT = Path("/cfs/earth/scratch/tabasran/MAP")
QUERY = None
MM2 = "minimap2"
MIN_ID = 0.9


def read_fasta(path, keep):
    seqs, name, buf = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if name in keep:
                seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name in keep:
        seqs[name] = "".join(buf)
    return seqs


def one(args):
    genome, fasta = args
    cmd = [MM2, "-c", "-x", "asm20", "-t", "1", "-N", "50", "--secondary=yes", fasta, QUERY]
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    if p.returncode != 0:
        return genome, None, p.stderr[-300:]
    iv, qlen = defaultdict(list), {}
    for line in p.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 12:
            continue
        q, ql, qs, qe = f[0], int(f[1]), int(f[2]), int(f[3])
        if int(f[10]) == 0 or int(f[9]) / int(f[10]) < MIN_ID:
            continue
        qlen[q] = ql
        iv[q].append((qs, qe))
    cov = {}
    for q, ints in iv.items():
        ints.sort()
        tot, cs, ce = 0, None, None
        for s, e in ints:
            if cs is None:
                cs, ce = s, e
            elif s <= ce:
                ce = max(ce, e)
            else:
                tot += ce - cs
                cs, ce = s, e
        tot += ce - cs
        cov[q] = round(tot / qlen[q], 4)
    return genome, cov, ""


def main():
    global QUERY, MM2, MIN_ID
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes", default=str(ROOT / "04_pangenome/panaroo/gene_presence_absence_roary.csv"))
    ap.add_argument("--ref", default=str(ROOT / "04_pangenome/panaroo/pan_genome_reference.fa"))
    ap.add_argument("--fasta", default=str(ROOT / "03_genomes/fasta"))
    ap.add_argument("--runs", default=str(ROOT / "06_pangwas/scoary_runs"))
    ap.add_argument("--batch", default=str(ROOT / "06_pangwas/scoary_summary/batch_check.tsv"))
    ap.add_argument("--ena", default=str(ROOT / "01_metadata/ena_raw.tsv"))
    ap.add_argument("--out", default=str(ROOT / "06_pangwas/scoary_summary"))
    ap.add_argument("--primary", default="C_nonB")
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--minimap2", default="minimap2")
    ap.add_argument("--min-id", type=float, default=0.9)
    ap.add_argument("--min-cov", type=float, default=0.8)
    ap.add_argument("--n-base", type=int, default=100)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()
    MM2, MIN_ID = args.minimap2, args.min_id
    out = Path(args.out)

    v = subprocess.run([MM2, "--version"], stdout=subprocess.PIPE, universal_newlines=True)
    if v.returncode != 0:
        sys.exit("[!] minimap2 not runnable: " + MM2)
    print("minimap2 " + v.stdout.strip())

    gm = pd.read_csv(args.genes, dtype=str, keep_default_na=False)
    iso = list(gm.columns[14:])
    gm = gm.set_index("Gene")
    ann = (gm[iso] != "").astype(int)
    freq = ann.mean(axis=1)

    b = pd.read_csv(args.batch, sep="\t")
    b["robust"] = (b["within_platform_supported"] == "yes") & (b["batch_suspect"] != "yes")
    cat = {}
    for _, r in b.iterrows():
        c = "hit_robust" if r["robust"] else "hit_nonrobust_" + r["direction"]
        if cat.get(r["Gene"]) != "hit_robust":
            cat[r["Gene"]] = c
    rnd = random.Random(args.seed)
    core = sorted(g for g in freq.index[freq >= 0.99] if g not in cat)
    acc = sorted(g for g in freq.index[(freq >= 0.05) & (freq <= 0.95)] if g not in cat)
    for g in rnd.sample(core, min(args.n_base, len(core))):
        cat[g] = "baseline_core"
    for g in rnd.sample(acc, min(args.n_base, len(acc))):
        cat[g] = "baseline_accessory"

    seqs = read_fasta(args.ref, set(cat))
    lost = sorted(set(cat) - set(seqs))
    print("genes: {} selected, {} with a reference sequence ({} missing: {})".format(
        len(cat), len(seqs), len(lost), lost[:10]))
    QUERY = str(out / "seq_presence_query.fa")
    with open(QUERY, "w") as fh:
        for g, s in seqs.items():
            fh.write(">{}\n{}\n".format(g, s))

    jobs = [(g, str(Path(args.fasta) / (g + ".fasta"))) for g in iso]
    covs, errs = {}, []
    pool = Pool(args.threads)
    for k, (g, cov, err) in enumerate(pool.imap_unordered(one, jobs, chunksize=4), 1):
        if cov is None:
            errs.append((g, err))
        else:
            covs[g] = cov
        if k % 300 == 0:
            print("  aligned {} / {}".format(k, len(jobs)))
    pool.close()
    if errs:
        sys.exit("[!] minimap2 failed for {} genomes, e.g. {}".format(len(errs), errs[:2]))
    C = pd.DataFrame({g: pd.Series(covs[g], dtype=float) for g in iso}).reindex(list(seqs)).fillna(0.0)
    C.to_csv(out / "seq_coverage.tsv.gz", sep="\t")
    S = (C >= args.min_cov).astype(int)

    # ------------------------------------------------ agreement by category
    order = ["baseline_core", "baseline_accessory", "hit_robust",
             "hit_nonrobust_enriched", "hit_nonrobust_depleted"]
    L = ["sequence presence: minimap2 asm20, identity >= {}, covered >= {} (merged across contigs)".format(
        args.min_id, args.min_cov), ""]
    L.append("{:<24} {:>6} {:>12} {:>14} {:>14}".format(
        "category", "genes", "agreement %", "ann0_seq1 %", "ann1_seq0 %"))
    for c in order:
        gs = [g for g in seqs if cat[g] == c]
        if not gs:
            continue
        a, s = ann.loc[gs, iso].values, S.loc[gs, iso].values
        n = a.size
        L.append("{:<24} {:>6} {:>11.1f}% {:>13.2f}% {:>13.2f}%".format(
            c, len(gs), 100.0 * (a == s).sum() / n,
            100.0 * ((a == 0) & (s == 1)).sum() / n, 100.0 * ((a == 1) & (s == 0)).sum() / n))

    # ------------------------------------------------ deltas per hit
    ena = pd.read_csv(args.ena, sep="\t", dtype=str).drop_duplicates("run_accession").set_index("run_accession")
    rows = []
    for trait, hs in b.groupby("trait"):
        tr = pd.read_csv(Path(args.runs) / "{}__{}".format(args.primary, trait) / "trait.csv",
                         dtype=str).set_index("Name").iloc[:, 0].astype(int)
        ids = list(tr.index)
        miseq = pd.Series([ena.at[i, "instrument_model"] == "Illumina MiSeq" if i in ena.index else False
                           for i in ids], index=ids)
        for _, h in hs.iterrows():
            g = h["Gene"]
            if g not in seqs:
                continue
            r = {"trait": trait, "Gene": g, "Annotation": h["Annotation"], "category": cat[g],
                 "len_bp": len(seqs[g])}
            for lab, M in (("ann", ann), ("seq", S)):
                x = M.loc[g, ids]
                r["delta_" + lab] = round(x[tr == 1].mean() - x[tr == 0].mean(), 3)
                xm, tm = x[miseq], tr[miseq]
                r["delta_" + lab + "_MiSeq"] = round(xm[tm == 1].mean() - xm[tm == 0].mean(), 3) \
                    if (tm == 1).sum() >= 20 and (tm == 0).sum() >= 20 else float("nan")
            r["freq_ann"] = round(ann.loc[g, ids].mean(), 3)
            r["freq_seq"] = round(S.loc[g, ids].mean(), 3)
            da, ds = r["delta_ann"], r["delta_seq"]
            r["retained"] = "yes" if (da * ds > 0 and abs(ds) >= 0.5 * abs(da)) else "no"
            rows.append(r)
    H = pd.DataFrame(rows)
    H.to_csv(out / "seq_presence_hits.tsv", sep="\t", index=False)

    L.append("\nper trait: hits whose association survives sequence-based presence "
             "(same sign, >= 50% of the Panaroo delta)")
    for t, g in H.groupby("trait"):
        L.append("  {:<22} hits {:>4}   retained {:>4}   median |delta| ann {:.3f} -> seq {:.3f}"
                 "   (MiSeq: ann {:.3f} -> seq {:.3f})".format(
                     t, len(g), int((g["retained"] == "yes").sum()),
                     g["delta_ann"].abs().median(), g["delta_seq"].abs().median(),
                     g["delta_ann_MiSeq"].abs().median(), g["delta_seq_MiSeq"].abs().median()))
    keep = H[(H["retained"] == "yes")].sort_values("delta_seq", key=abs, ascending=False) \
        if len(H) else H
    pd.set_option("display.width", 250)
    pd.set_option("display.max_colwidth", 45)
    L.append("\n--- hits retained on sequence (top 40 by |delta_seq|) ---")
    L.append(keep.head(40)[["trait", "Gene", "Annotation", "category", "len_bp", "freq_ann", "freq_seq",
                            "delta_ann", "delta_seq", "delta_ann_MiSeq", "delta_seq_MiSeq"]].to_string(index=False))
    txt = "\n".join(L)
    (out / "seq_presence_summary.txt").write_text(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()
