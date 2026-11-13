#!/usr/bin/env python3
"""
46_split_markers.py

A fourth typing level: SNP markers for the major split inside C/non-B.

44/45 showed that ~1161 C/non-B genomes (96% of the Americas, 95% of the
Oceania C isolates, incl. K-10) form one DERIVED clade 'D', while ~508 mostly
European C/non-B genomes keep the ancestral alleles shared with S and Bison
(a paraphyletic basal grade). The split is carried by >= 10 clade-compatible
SNPs. This script finds every core SNP that marks clade D and ranks them as
typing markers, in the same spirit as the 8-SNP gyrA/gyrB scheme (16).

Steps
  1. clade D = the tree bipartition closest to the carrier set of --anchor
     (default AE016958.1_374739_T_C), taken on the side WITHOUT lineage S
  2. for every core SNP allele: mismatches with D (missing calls ignored; an
     allele may mark D or its complement), missing calls, Fitch changes
     (05_phylogeny/homoplasy/fitch_all_sites.tsv), annotation in K-10 (CDS,
     codon change), other core SNPs within +-50 bp (primer design), distance
     to the nearest masked repeat; mismatching genomes flagged if they are
     possible mixed samples (08_reads/summary/genome_qc_reads.tsv)
  3. two best markers (fewest mismatches outside mixed samples, no missing,
     no neighbouring SNP, not near a repeat) combined into a two-SNP call and
     compared with the tree clade; composition of D and of the basal grade by
     region, host and lineage

  python 46_split_markers.py [--anchor V] [--max-mismatch 5]
Outputs: 05_phylogeny/typing/split_markers.tsv, split_markers_report.txt,
         split_assignments.tsv (genome, tree side, two-SNP call)
Environment: ml_env (Python >= 3.10). Imports 44_homoplasy_site.py.
"""

import argparse
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/cfs/earth/scratch/tabasran/MAP")


def load44(path):
    sp = importlib.util.spec_from_file_location("s44", path)
    m = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(m)
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--anchor", default="AE016958.1_374739_T_C")
    ap.add_argument("--max-mismatch", type=int, default=5)
    ap.add_argument("--repeat-min", type=int, default=1000, help="min bp to a masked repeat for a good marker")
    ap.add_argument("--mixed", type=int, default=20)
    ap.add_argument("--tree", default=str(ROOT / "05_phylogeny/iqtree/tree_final.nwk"))
    ap.add_argument("--coretab", default=str(ROOT / "05_phylogeny/core_masked/core.tab"))
    ap.add_argument("--fitch", default=str(ROOT / "05_phylogeny/homoplasy/fitch_all_sites.tsv"))
    ap.add_argument("--mask", default=str(ROOT / "05_phylogeny/core_masked/mask_genes.tsv"))
    ap.add_argument("--ref", default=str(ROOT / "03_genomes/fasta/GCA_000007865.1.fasta"))
    ap.add_argument("--gff", default=str(ROOT / "03_genomes/gff/GCA_000007865.1.gff"))
    ap.add_argument("--lineages", default=str(ROOT / "05_phylogeny/typing/lineage_assignments.tsv"))
    ap.add_argument("--traits", default=str(ROOT / "06_pangwas/trait.csv"))
    ap.add_argument("--ena", default=str(ROOT / "01_metadata/ena_raw.tsv"))
    ap.add_argument("--readqc", default=str(ROOT / "08_reads/summary/genome_qc_reads.tsv"))
    ap.add_argument("--s44", default=str(Path(__file__).resolve().parent / "44_homoplasy_site.py"))
    ap.add_argument("--out", default=str(ROOT / "05_phylogeny/typing"))
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    s44 = load44(args.s44)
    R = []

    def say(s=""):
        print(s, flush=True)
        R.append(s)

    # ------------------------------------------------ tree bipartitions
    adj, name, sup = s44.parse_newick(open(args.tree).read())
    tips = sorted(name.values())
    bit = {t: 1 << k for k, t in enumerate(tips)}
    ALL = (1 << len(tips)) - 1
    start = next(iter(adj))
    parent, order, st = {start: None}, [], [start]
    while st:
        x = st.pop()
        order.append(x)
        for y in adj[x]:
            if y != parent[x]:
                parent[y] = x
                st.append(y)
    mask = {}
    for x in reversed(order):
        mask[x] = bit[name[x]] if x in name else 0
        for y in adj[x]:
            if y != parent[x]:
                mask[x] |= mask[y]
    edges = [(mask[x], sup.get(frozenset((x, parent[x])), "")) for x in order if parent[x] is not None]

    lin = pd.read_csv(args.lineages, sep="\t", dtype=str).set_index("genome")
    lab = (lin["lineage"] + "/" + lin["sublineage"]).reindex(tips).fillna("?")
    Smask = sum(bit[t] for t in tips if lab[t].startswith("S/"))
    tr = pd.read_csv(args.traits, dtype=str).set_index("Name")
    ena = pd.read_csv(args.ena, sep="\t", dtype=str).drop_duplicates("run_accession").set_index("run_accession")
    q = pd.read_csv(args.readqc, sep="\t", index_col=0) if Path(args.readqc).exists() else pd.DataFrame()
    mixed = set(q.index[q["mixed_sites"] >= args.mixed]) if "mixed_sites" in q.columns else set()

    # ------------------------------------------------ genotypes
    ct = pd.read_csv(args.coretab, sep="\t", dtype=str)
    refcol = "REF" if "REF" in ct.columns else "Reference"
    ct["POS"] = ct["POS"].astype(int)
    Bm = ct[tips].apply(lambda c: c.str.upper())
    BA = Bm.fillna("N").values                      # sites x tips, same order as `tips`
    allpos = np.sort(ct["POS"].values)

    def to_int(arr):
        """bool array over tips -> Python int with bit k = tips[k]"""
        return int.from_bytes(np.packbits(arr, bitorder="little").tobytes(), "little")

    def allele_masks(k, alt):
        row = BA[k]
        return to_int(row == alt), to_int(~np.isin(row, ["N", "-", "?"]))

    chrom, apos, aref, aalt = args.anchor.rsplit("_", 3)
    ka = ct.index[ct["POS"] == int(apos)]
    if not len(ka):
        sys.exit("[!] anchor {} not in core.tab".format(args.anchor))
    car, known = allele_masks(ka[0], aalt)
    best = None
    for m, s in edges:
        for side in (m, ALL ^ m):
            d = ((side ^ car) & known).bit_count()
            if best is None or d < best[0]:
                best = (d, side, s)
    _, side, sD = best
    D = side if not (side & Smask) else ALL ^ side
    if D & Smask:
        sys.exit("[!] neither side of the closest bipartition is free of lineage S")
    Dset = [t for t in tips if D & bit[t]]
    say("clade D: {} genomes (UFBoot {}), K-10 inside: {}".format(
        len(Dset), sD, "GCA_000007865.1" in Dset))

    # composition
    reg = pd.Series("unknown", index=tips)
    for c in ("Region_Americas", "Region_Europe", "Region_Oceania"):
        if c in tr.columns:
            reg[[t for t in tips if tr[c].get(t) == "1"]] = c.replace("Region_", "")
    host = pd.Series("NA", index=tips)
    if "hostCattle" in tr.columns:
        host[[t for t in tips if tr["hostCattle"].get(t) == "1"]] = "cattle"
        host[[t for t in tips if tr["hostCattle"].get(t) == "0"]] = "non-cattle"
    side_lab = pd.Series(np.where([D & bit[t] for t in tips], "D (derived clade)",
                                  np.where(lab.str.startswith("C/non-B"), "C/non-B basal grade", lab)), index=tips)
    say("\ncomposition by region:\n" + pd.crosstab(side_lab, reg).to_string())
    say("\ncomposition by host:\n" + pd.crosstab(side_lab, host).to_string())
    study = pd.Series([ena.at[t, "study_accession"] if t in ena.index else "assembly" for t in tips], index=tips)
    say("\nstudies represented: D {} | basal grade {}".format(
        study[side_lab == "D (derived clade)"].nunique(), study[side_lab == "C/non-B basal grade"].nunique()))

    # ------------------------------------------------ marker scan
    seq = {}
    cur = None
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
    mk = pd.read_csv(args.mask, sep="\t", header=None, names=["chr", "s0", "e", "product"]) \
        if Path(args.mask).exists() else None
    F = pd.read_csv(args.fitch, sep="\t").set_index("variant") if Path(args.fitch).exists() else None

    def annotate(pos, ref, alt):
        for s, e, strand, lt, gene, prod in cds:
            if s <= pos <= e:
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
                eff = "syn" if a1 == a2 else ("stop" if a2 == "*" else "missense")
                return "{} {} {}{}{} {}".format(lt, gene or "-", a1, off // 3 + 1, a2, eff), prod
        return "intergenic", ""

    rows = []
    for k in range(len(ct)):
        pos = int(ct.at[k, "POS"])
        ref = str(ct.at[k, refcol]).upper()
        for alt in sorted(set(BA[k]) - {ref, "N", "-", "?"}):
            car, known = allele_masks(k, alt)
            d1 = ((car ^ D) & known).bit_count()
            d2 = ((car ^ (ALL ^ D)) & known).bit_count()
            if min(d1, d2) > args.max_mismatch:
                continue
            marks = alt if d1 <= d2 else ref
            diff = ((car ^ D) if d1 <= d2 else (car ^ (ALL ^ D))) & known
            odd = [t for t in tips if diff & bit[t]]
            v = "{}_{}_{}_{}".format(ct.at[k, "CHR"], pos, ref, alt)
            locus, prod = annotate(pos, ref, alt)
            near = int(((allpos >= pos - 50) & (allpos <= pos + 50)).sum()) - 1
            if mk is not None:
                dd = np.where(pos < mk["s0"] + 1, mk["s0"] + 1 - pos, np.where(pos > mk["e"], pos - mk["e"], 0))
                drep = int(dd.min())
            else:
                drep = -1
            rows.append({"variant": v, "D_allele": marks, "other_allele": ref if marks == alt else alt,
                         "mismatches": len(odd), "mismatches_not_mixed": sum(t not in mixed for t in odd),
                         "missing": len(tips) - known.bit_count(),
                         "fitch": int(F.at[v, "fitch_changes"]) if F is not None and v in F.index else np.nan,
                         "snps_within_50bp": near, "bp_to_repeat": drep, "locus": locus, "product": prod,
                         "mismatching_genomes": ",".join(odd)})
    M = pd.DataFrame(rows)
    if M.empty:
        sys.exit("[!] no SNP within {} mismatches of clade D".format(args.max_mismatch))
    M["good"] = (M["mismatches_not_mixed"] == 0) & (M["missing"] == 0) & (M["snps_within_50bp"] == 0) & \
                (M["bp_to_repeat"] >= args.repeat_min)
    M = M.sort_values(["good", "mismatches_not_mixed", "mismatches", "missing", "snps_within_50bp"],
                      ascending=[False, True, True, True, True])
    M.to_csv(out / "split_markers.tsv", sep="\t", index=False)
    pd.set_option("display.width", 260)
    pd.set_option("display.max_colwidth", 50)
    say("\nSNPs marking clade D with <= {} mismatching genomes: {} (good markers: {})".format(
        args.max_mismatch, len(M), int(M["good"].sum())))
    say(M.drop(columns="mismatching_genomes").to_string(index=False))

    # ------------------------------------------------ two-SNP call
    top = M[M["good"]].head(2)
    if len(top) < 2:
        top = M.head(2)
        say("\n[i] fewer than two 'good' markers; using the two best ranked")
    calls = {}
    for _, r in top.iterrows():
        pos = int(r["variant"].rsplit("_", 3)[1])
        b = Bm.iloc[int(ct.index[ct["POS"] == pos][0])]
        calls[r["variant"]] = b.map(lambda x: "D" if x == r["D_allele"] else ("?" if x in ("N", "-") or not isinstance(x, str) else "not-D"))
    C = pd.DataFrame(calls)
    two = np.where((C == "D").all(axis=1), "D", np.where((C == "not-D").all(axis=1), "not-D", "conflict/missing"))
    A = pd.DataFrame({"genome": tips, "tree_side": side_lab.values, "lineage": lab.values,
                      "two_snp_call": two})
    A.to_csv(out / "split_assignments.tsv", sep="\t", index=False)
    say("\ntwo-SNP call ({}) vs tree:".format(" + ".join(top["variant"])))
    say(pd.crosstab(A["tree_side"], A["two_snp_call"]).to_string())
    (out / "split_markers_report.txt").write_text("\n".join(R) + "\n")
    say("\nwritten: {}/split_markers.tsv, split_assignments.tsv, split_markers_report.txt".format(out))


if __name__ == "__main__":
    main()
