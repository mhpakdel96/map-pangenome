#!/usr/bin/env python3
"""
45_discordant_tips.py

Which genomes make a SNP look homoplastic?

44 showed that the 21 'homoplastic' SNPs of 30/31 need only 1-2 extra
changes on the ML tree (Fitch, root-independent; genome-wide CI 0.966).
The earlier 'origins' count was inflated because, on a rooted tree, carriers
of an ANCESTRAL allele always fall into several groups. For 374739 the two
extra genomes were one possible mixed sample (SRR11839074, reads 60% ALT)
and one assembly (GCA_041877225.1). This script does the same for every
site, root-free:

  for each SNP allele, the tree bipartition (edge) closest to its carrier
  set is found (smallest symmetric difference, missing calls ignored); the
  genomes in that difference are the 'discordant tips' that the SNP would
  have to move to become a clean clade marker. For each discordant tip:
  lineage, study, read-based mixed-site count (35), N50, and the ALT allele
  fraction at that site from reads (08_reads/k10/<ID>.ad.tsv).

Default variants: those in 06_pangwas/pyseer/homoplasy_check.tsv.
  python 45_discordant_tips.py [--variants V1,V2,...] [--min-fitch 2 --all]
Outputs: 05_phylogeny/homoplasy/discordant_tips.tsv, discordant_summary.txt
Environment: ml_env (Python >= 3.10).
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
    ap.add_argument("--variants", default="")
    ap.add_argument("--all", action="store_true", help="all alleles with fitch_changes >= --min-fitch")
    ap.add_argument("--min-fitch", type=int, default=2)
    ap.add_argument("--hom", default=str(ROOT / "06_pangwas/pyseer/homoplasy_check.tsv"))
    ap.add_argument("--fitch", default=str(ROOT / "05_phylogeny/homoplasy/fitch_all_sites.tsv"))
    ap.add_argument("--tree", default=str(ROOT / "05_phylogeny/iqtree/tree_final.nwk"))
    ap.add_argument("--coretab", default=str(ROOT / "05_phylogeny/core_masked/core.tab"))
    ap.add_argument("--lineages", default=str(ROOT / "05_phylogeny/typing/lineage_assignments.tsv"))
    ap.add_argument("--ena", default=str(ROOT / "01_metadata/ena_raw.tsv"))
    ap.add_argument("--readqc", default=str(ROOT / "08_reads/summary/genome_qc_reads.tsv"))
    ap.add_argument("--reads", default=str(ROOT / "08_reads/k10"))
    ap.add_argument("--mixed", type=int, default=20)
    ap.add_argument("--s44", default=str(Path(__file__).resolve().parent / "44_homoplasy_site.py"))
    ap.add_argument("--out", default=str(ROOT / "05_phylogeny/homoplasy"))
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    s44 = load44(args.s44)
    R = []

    def say(s=""):
        print(s, flush=True)
        R.append(s)

    F = pd.read_csv(args.fitch, sep="\t").set_index("variant") if Path(args.fitch).exists() else None
    if args.variants:
        V = [v.strip() for v in args.variants.split(",") if v.strip()]
    elif args.all:
        if F is None:
            sys.exit("[!] --all needs {} (run 44 first)".format(args.fitch))
        V = F.index[(F["fitch_changes"] >= args.min_fitch) & (F["minor_count"] >= 2)].tolist()
    else:
        V = pd.read_csv(args.hom, sep="\t")["variant"].unique().tolist()
    say("variants: {}".format(len(V)))

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
        kids = [y for y in adj[x] if y != parent[x]]
        mask[x] = bit[name[x]] if x in name else 0
        for y in kids:
            mask[x] |= mask[y]
    edges = [(x, mask[x], sup.get(frozenset((x, parent[x])), "")) for x in order if parent[x] is not None]
    say("tree: {} tips, {} edges".format(len(tips), len(edges)))

    # ------------------------------------------------ genotypes
    ct = pd.read_csv(args.coretab, sep="\t", dtype=str)
    ct["POS"] = ct["POS"].astype(int)
    lin = pd.read_csv(args.lineages, sep="\t", dtype=str).set_index("genome")
    lab = lin["lineage"] + "/" + lin["sublineage"]
    ena = pd.read_csv(args.ena, sep="\t", dtype=str).drop_duplicates("run_accession").set_index("run_accession")
    q = pd.read_csv(args.readqc, sep="\t", index_col=0) if Path(args.readqc).exists() else pd.DataFrame()

    def read_af(g, pos, alt):
        f = Path(args.reads) / (g + ".ad.tsv")
        if not f.exists():
            return np.nan
        for line in open(f):
            p = line.split("\t")
            if len(p) >= 4 and int(p[0]) == pos:
                ad = [int(v) for v in p[3].strip().split(",") if v not in (".", "")]
                alts = p[2].split(",")
                tot = sum(ad)
                if not tot:
                    return np.nan
                return round((ad[1 + alts.index(alt)] if alt in alts and 1 + alts.index(alt) < len(ad) else 0) / tot, 2)
        return np.nan

    rows, per_site = [], []
    for v in V:
        chrom, pos, ref, alt = v.rsplit("_", 3)
        pos = int(pos)
        r = ct[ct["POS"] == pos]
        if r.empty:
            say("[!] {} not in core.tab".format(v))
            continue
        b = r.iloc[0][tips].str.upper()
        car = sum(bit[t] for t in tips if b[t] == alt)
        known = sum(bit[t] for t in tips if b[t] not in ("N", "-", "?") and isinstance(b[t], str))
        best = None
        for x, m, s in edges:
            for side in (m, ALL ^ m):
                d = ((side ^ car) & known).bit_count()
                if best is None or d < best[0]:
                    best = (d, side, s)
        d, side, s = best
        diff = (side ^ car) & known
        odd = [t for t in tips if diff & bit[t]]
        nmix = 0
        for g in odd:
            mixed = q["mixed_sites"].get(g, np.nan) if "mixed_sites" in q.columns else np.nan
            nmix += int(mixed >= args.mixed) if mixed == mixed else 0
            rows.append({"variant": v, "genome": g,
                         "state": "carrier outside clade" if car & bit[g] else "non-carrier inside clade",
                         "lineage": lab.get(g, "?"),
                         "study": ena.at[g, "study_accession"] if g in ena.index else "assembly",
                         "mixed_sites": mixed,
                         "n50": q["n50"].get(g, np.nan) if "n50" in q.columns else np.nan,
                         "reads_alt_fraction": read_af(g, pos, alt)})
        per_site.append({"variant": v, "carriers": bin(car).count("1"),
                         "fitch": int(F.at[v, "fitch_changes"]) if F is not None and v in F.index else np.nan,
                         "closest_clade_size": bin(side).count("1"), "clade_support": s,
                         "discordant_tips": d, "of_them_mixed": nmix,
                         "tips": ",".join(odd[:8]) + (",..." if len(odd) > 8 else "")})

    D = pd.DataFrame(rows)
    P = pd.DataFrame(per_site)
    D.to_csv(out / "discordant_tips.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_colwidth", 90)
    say("\nper site: closest tree clade and the genomes that break it")
    say(P.to_string(index=False))
    if len(D):
        say("\ngenomes breaking at least one site:")
        G = D.groupby("genome").agg(sites=("variant", "nunique"), lineage=("lineage", "first"),
                                    study=("study", "first"), mixed_sites=("mixed_sites", "first"),
                                    n50=("n50", "first"),
                                    reads_alt_fraction=("reads_alt_fraction",
                                                        lambda s: ",".join(str(x) for x in s))).sort_values(
            "sites", ascending=False)
        say(G.to_string())
        say("\nsites fully explained by possible mixed samples (>= {} mixed sites): {} of {}".format(
            args.mixed, int((P["discordant_tips"] == P["of_them_mixed"]).sum()), len(P)))
        say("discordant genomes: {} distinct | possible mixed: {} | without reads: {}".format(
            len(G), int((G["mixed_sites"] >= args.mixed).sum()), int(G["mixed_sites"].isna().sum())))
    (out / "discordant_summary.txt").write_text("\n".join(R) + "\n")
    say("\nwritten: {}/discordant_tips.tsv, discordant_summary.txt".format(out))


if __name__ == "__main__":
    main()
