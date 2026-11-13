#!/usr/bin/env python3
"""
24_dna_pangenome.py

Pangenome with gene presence called from DNA, not from gene annotation.

Why. 23_sequence_presence.py showed that Panaroo's presence/absence (from
Prokka gene calls on short-read assemblies) is wrong for ~25% of the
'absences' of ordinary accessory genes and for >50% of those of the Scoary
hits: the DNA is in the assembly, but the gene was not called (typically
broken at a contig end). Accessory size, the core/shell/cloud split and the
open/closed verdict all depend on these calls, so they are recomputed here.

Method (every rule is a parameter, printed in the report)
  1. Query = every Panaroo cluster's representative sequence
     (pan_genome_reference.fa).
  2. REDUNDANT FRAGMENTS: a cluster whose sequence lies >= --min-cov inside a
     LONGER cluster at >= --min-id (minimap2 self-alignment) is a fragment
     of that gene (Panaroo keeps split genes as separate clusters). Once
     presence is called from DNA it would just duplicate its parent, so it is
     excluded from the DNA pangenome and listed with its parent.
  3. DNA PRESENCE: minimap2 (-x asm20) of all queries against every assembly;
     aligned query intervals with identity >= --min-id are merged ACROSS
     contigs; present = covered fraction >= --min-cov.
  4. FLAG 'dna_below_annotation': DNA frequency lower than annotation
     frequency by more than --rep-gap. In a clonal species this means the
     representative sequence is unusual (chimera, junction, divergent
     allele), not that the DNA is missing; such clusters are reported and
     shown with and without in the counts.
  Categories as Panaroo: core >= 99%, soft core 95-99%, shell 15-95%,
  cloud < 15%.

Also reported
  * transitions annotation category -> DNA category
  * accessory genes per genome vs N50 and platform (Spearman), annotation vs
    DNA: if annotation counts track assembly quality and DNA counts do not,
    the batch effect is removed
  * pangenome growth (Heaps' law P(n) = k n^gamma, gamma > 0 = open), on the
    same gene set for annotation and DNA, mean of --perm permutations

Outputs (04_pangenome/dna_presence/)
  coverage.tsv.gz        cluster x genome covered fraction (all clusters)
  presence_dna.Rtab      0/1, clusters kept in the DNA pangenome
  gene_table.tsv         per cluster: length, frequencies, categories, flags
  dna_pangenome_report.txt
"""

import argparse
import subprocess
import sys
from collections import defaultdict
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/cfs/earth/scratch/tabasran/MAP")
QUERY, MM2, MIN_ID = None, "minimap2", 0.9


def read_fasta(path):
    seqs, name, buf = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def merge_cov(ints, qlen):
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
    if cs is not None:
        tot += ce - cs
    return tot / float(qlen)


def paf_cov(text, min_id):
    iv, ql = defaultdict(list), {}
    for line in text.splitlines():
        f = line.split("\t")
        if len(f) < 12 or int(f[10]) == 0 or int(f[9]) / float(f[10]) < min_id:
            continue
        iv[f[0]].append((int(f[2]), int(f[3])))
        ql[f[0]] = int(f[1])
    return {q: round(merge_cov(v, ql[q]), 4) for q, v in iv.items()}


def one(job):
    genome, fasta = job
    p = subprocess.run([MM2, "-c", "-x", "asm20", "-t", "1", "-N", "50", "--secondary=yes",
                        fasta, QUERY], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                       universal_newlines=True)
    if p.returncode != 0:
        return genome, None, p.stderr[-300:]
    return genome, paf_cov(p.stdout, MIN_ID), ""


def category(f):
    return np.select([f >= 0.99, f >= 0.95, f >= 0.15, f > 0],
                     ["core", "soft_core", "shell", "cloud"], "absent")


def heaps(M, perms, rng):
    """M: genes x genomes bool. Returns n grid, mean pangenome size, gamma."""
    n = M.shape[1]
    grid = sorted(set([1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, n]) & set(range(1, n + 1)))
    sizes = np.zeros((perms, len(grid)))
    for p in range(perms):
        cum = np.zeros(M.shape[0], dtype=bool)
        order = rng.permutation(n)
        k = 0
        for j, c in enumerate(order, 1):
            cum |= M[:, c]
            if j == grid[k]:
                sizes[p, k] = cum.sum()
                k += 1
                if k == len(grid):
                    break
    mean = sizes.mean(axis=0)
    g = np.array(grid, dtype=float)
    sel = g >= 10
    gamma = float(np.polyfit(np.log(g[sel]), np.log(mean[sel]), 1)[0])
    return grid, mean, gamma


def main():
    global QUERY, MM2, MIN_ID
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes", default=str(ROOT / "04_pangenome/panaroo/gene_presence_absence_roary.csv"))
    ap.add_argument("--ref", default=str(ROOT / "04_pangenome/panaroo/pan_genome_reference.fa"))
    ap.add_argument("--fasta", default=str(ROOT / "03_genomes/fasta"))
    ap.add_argument("--qc", default=str(ROOT / "02_qc/genome_qc_master.tsv"))
    ap.add_argument("--ena", default=str(ROOT / "01_metadata/ena_raw.tsv"))
    ap.add_argument("--out", default=str(ROOT / "04_pangenome/dna_presence"))
    ap.add_argument("--minimap2", default="minimap2")
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--min-id", type=float, default=0.9)
    ap.add_argument("--min-cov", type=float, default=0.8)
    ap.add_argument("--rep-gap", type=float, default=0.2)
    ap.add_argument("--perm", type=int, default=50)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()
    MM2, MIN_ID = args.minimap2, args.min_id
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    R = []

    def say(s=""):
        print(s)
        R.append(s)

    v = subprocess.run([MM2, "--version"], stdout=subprocess.PIPE, universal_newlines=True)
    if v.returncode != 0:
        sys.exit("[!] minimap2 not runnable: " + MM2)
    say("minimap2 {} | identity >= {} | covered >= {} | rep-gap {}".format(
        v.stdout.strip(), args.min_id, args.min_cov, args.rep_gap))

    gm = pd.read_csv(args.genes, dtype=str, keep_default_na=False)
    iso = list(gm.columns[14:])
    gm = gm.set_index("Gene")
    ann = (gm[iso] != "").values
    genes = list(gm.index)
    f_ann = pd.Series(ann.mean(axis=1), index=genes)

    seqs = read_fasta(args.ref)
    noref = [g for g in genes if g not in seqs]
    extra = [g for g in seqs if g not in gm.index]
    say("clusters: {} | with reference sequence: {} | no reference: {} | in reference only: {}".format(
        len(genes), len(genes) - len(noref), len(noref), len(extra)))
    QUERY = str(out / "query.fa")
    with open(QUERY, "w") as fh:
        for g in genes:
            if g in seqs:
                fh.write(">{}\n{}\n".format(g, seqs[g]))

    # ---------------------------------------------- 2. redundant fragments
    p = subprocess.run([MM2, "-c", "-x", "asm20", "-X", "-t", str(args.threads), "-N", "50",
                        QUERY, QUERY], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                       universal_newlines=True)
    if p.returncode != 0:
        sys.exit("[!] self-alignment failed: " + p.stderr[-300:])
    iv = defaultdict(list)
    for line in p.stdout.splitlines():
        f = line.split("\t")
        q, t = f[0], f[5]
        if q == t or int(f[10]) == 0 or int(f[9]) / float(f[10]) < args.min_id:
            continue
        if len(seqs[t]) <= len(seqs[q]):
            continue                      # only containment in a LONGER cluster
        iv[(q, t)].append((int(f[2]), int(f[3])))
    parent = {}
    for (q, t), ints in iv.items():
        c = merge_cov(ints, len(seqs[q]))
        if c >= args.min_cov and (q not in parent or len(seqs[t]) > len(seqs[parent[q]])):
            parent[q] = t
    say("redundant fragments (inside a longer cluster): {}".format(len(parent)))

    # ---------------------------------------------- 3. DNA presence
    jobs = [(g, str(Path(args.fasta) / (g + ".fasta"))) for g in iso]
    covs, errs = {}, []
    pool = Pool(args.threads)
    for k, (g, cov, err) in enumerate(pool.imap_unordered(one, jobs, chunksize=2), 1):
        if cov is None:
            errs.append((g, err))
        else:
            covs[g] = cov
        if k % 200 == 0:
            print("  aligned {} / {}".format(k, len(jobs)))
    pool.close()
    if errs:
        sys.exit("[!] minimap2 failed for {} genomes, e.g. {}".format(len(errs), errs[:2]))
    C = pd.DataFrame({g: pd.Series(covs[g], dtype=float) for g in iso}).reindex(genes).fillna(0.0)
    C.to_csv(out / "coverage.tsv.gz", sep="\t", float_format="%.3f")
    dna = (C.values >= args.min_cov)
    f_dna = pd.Series(dna.mean(axis=1), index=genes)

    # ---------------------------------------------- 4. table + flags
    T = pd.DataFrame({"gene": genes, "annotation": gm["Annotation"].values,
                      "len_bp": [len(seqs[g]) if g in seqs else np.nan for g in genes],
                      "freq_ann": f_ann.round(4).values, "freq_dna": f_dna.round(4).values})
    T["cat_ann"] = category(T["freq_ann"].values)
    T["cat_dna"] = category(T["freq_dna"].values)
    T["redundant_of"] = [parent.get(g, "") for g in genes]
    T["no_reference"] = [g in set(noref) for g in genes]
    T["dna_below_annotation"] = (T["freq_ann"] - T["freq_dna"]) > args.rep_gap
    T["kept"] = (~T["no_reference"]) & (T["redundant_of"] == "")
    T.to_csv(out / "gene_table.tsv", sep="\t", index=False)
    keep = T["kept"].values
    say("kept in DNA pangenome: {} | flagged dna_below_annotation among kept: {}".format(
        int(keep.sum()), int((keep & T["dna_below_annotation"].values).sum())))
    pd.DataFrame(dna[keep].astype(int), index=np.array(genes)[keep], columns=iso) \
        .to_csv(out / "presence_dna.Rtab", sep="\t", index_label="Gene")

    # ---------------------------------------------- category counts
    cats = ["core", "soft_core", "shell", "cloud", "absent"]

    def counts(col, mask):
        c = T.loc[mask, col].value_counts()
        return [int(c.get(x, 0)) for x in cats]
    flag = T["dna_below_annotation"].values
    say("\n{:<46} {:>6} {:>9} {:>6} {:>6} {:>7} {:>7}".format(
        "gene set / presence call", "core", "soft_core", "shell", "cloud", "absent", "total*"))
    for lab, col, mask in (
            ("Panaroo, all clusters (annotation)", "cat_ann", np.ones(len(T), bool)),
            ("kept clusters, annotation", "cat_ann", keep),
            ("kept clusters, DNA", "cat_dna", keep),
            ("kept & not flagged, annotation", "cat_ann", keep & ~flag),
            ("kept & not flagged, DNA", "cat_dna", keep & ~flag)):
        c = counts(col, mask)
        say("{:<46} {:>6} {:>9} {:>6} {:>6} {:>7} {:>7}".format(lab, *(c + [sum(c[:4])])))
    say("* total excludes 'absent' (cluster present in no genome by that call)")

    say("\ntransitions, kept clusters (rows annotation, columns DNA):")
    tr = pd.crosstab(T.loc[keep, "cat_ann"], T.loc[keep, "cat_dna"]).reindex(
        index=cats, columns=cats, fill_value=0)
    say(tr.to_string())

    say("\nlength of kept clusters by annotation category -> DNA category (median bp):")
    say(T[keep].groupby(["cat_ann", "cat_dna"])["len_bp"].agg(["count", "median"]).to_string())

    # ---------------------------------------------- per genome vs quality
    try:
        from scipy.stats import spearmanr
    except ImportError:
        spearmanr = None
    qc = pd.read_csv(args.qc, sep="\t").set_index("genome").reindex(iso)
    ena = pd.read_csv(args.ena, sep="\t", dtype=str).drop_duplicates("run_accession").set_index("run_accession")
    plat = pd.Series([ena.at[i, "instrument_model"] if i in ena.index else "assembly" for i in iso], index=iso)
    acc_ann = (T.loc[keep, "cat_ann"] != "core").values
    acc_dna = (T.loc[keep, "cat_dna"] != "core").values
    n_ann = pd.Series(ann[keep][acc_ann].sum(axis=0), index=iso)
    n_dna = pd.Series(dna[keep][acc_dna].sum(axis=0), index=iso)
    tot_ann = pd.Series(ann[keep].sum(axis=0), index=iso)
    tot_dna = pd.Series(dna[keep].sum(axis=0), index=iso)
    say("\nper genome (kept clusters): genes and accessory genes vs assembly quality")
    if spearmanr is not None:
        for lab, s in (("genes, annotation", tot_ann), ("genes, DNA", tot_dna),
                       ("accessory, annotation", n_ann), ("accessory, DNA", n_dna)):
            r1 = spearmanr(s, qc["n50"])
            r2 = spearmanr(s, qc["contigs"])
            say("  {:<24} Spearman vs N50 {:+.3f} (p={:.1e})   vs contigs {:+.3f} (p={:.1e})".format(
                lab, r1[0], r1[1], r2[0], r2[1]))
    say("  median per platform (genes ann / DNA | accessory ann / DNA):")
    d = pd.DataFrame({"plat": plat, "ga": tot_ann, "gd": tot_dna, "aa": n_ann, "ad": n_dna,
                      "n50": qc["n50"]})
    say(d.groupby("plat").agg(n=("ga", "size"), n50=("n50", "median"), genes_ann=("ga", "median"),
                              genes_dna=("gd", "median"), acc_ann=("aa", "median"),
                              acc_dna=("ad", "median")).sort_values("n50").round(0).to_string())

    # ---------------------------------------------- openness
    rng = np.random.RandomState(args.seed)
    say("\npangenome growth, Heaps' law P(n) = k n^gamma (n >= 10), mean of {} permutations".format(args.perm))
    for lab, M in (("kept, annotation", ann[keep]), ("kept, DNA", dna[keep]),
                   ("kept & not flagged, annotation", ann[keep & ~flag]),
                   ("kept & not flagged, DNA", dna[keep & ~flag])):
        grid, mean, gamma = heaps(M, args.perm, rng)
        say("  {:<32} gamma = {:.4f}   P(n): ".format(lab, gamma) + "  ".join(
            "{}:{:.0f}".format(a, b) for a, b in zip(grid, mean)))
    say("  gamma > 0 open, ~0 closed; compare annotation vs DNA on the SAME clusters")

    (out / "dna_pangenome_report.txt").write_text("\n".join(R) + "\n")
    say("\nwritten: {}/coverage.tsv.gz, presence_dna.Rtab, gene_table.tsv, dna_pangenome_report.txt".format(out))


if __name__ == "__main__":
    main()
