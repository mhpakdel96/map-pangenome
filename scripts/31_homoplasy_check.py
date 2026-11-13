#!/usr/bin/env python3
"""
31_homoplasy_check.py

The SNP hits that conflict with well-supported clades (30_snp_hits.py,
'homoplastic') need an explanation before any of them is called convergent.
Seven Region_Europe SNPs at seven distant positions share one identical
carrier set (506-507 genomes, 19-20 origins, 7-8 conflicts with UFBoot>=95
clades). Identical patterns at unlinked sites are not expected from
independent convergent mutations. Candidate explanations tested here:

  1. RECOMBINATION: the SNP lies in a Gubbins recombination block that
     affects (some of) its carriers. The ML tree was built on the
     Gubbins-filtered alignment, so recombined sites are expected to
     conflict with it.
  2. TECHNICAL BATCH: carriers concentrated in one platform/study beyond what
     region explains (e.g. a systematic miscall on one chemistry).
  3. SHARED PATTERN: SNPs with identical carrier sets are grouped, so one
     event is not counted many times.

Outputs: 06_pangwas/pyseer/homoplasy_check.tsv and a printed summary.
"""

import argparse
import re
from pathlib import Path

import pandas as pd

ROOT = Path("/cfs/earth/scratch/tabasran/MAP")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=str(ROOT / "06_pangwas/pyseer"))
    ap.add_argument("--gubbins", default=str(ROOT / "05_phylogeny/gubbins/gubbins.recombination_predictions.gff"))
    ap.add_argument("--ena", default=str(ROOT / "01_metadata/ena_raw.tsv"))
    ap.add_argument("--traits", default=str(ROOT / "06_pangwas/trait.csv"))
    args = ap.parse_args()
    d = Path(args.dir)

    H = pd.read_csv(d / "snp_hits_annotated.tsv", sep="\t")
    H = H[(H["subset"] == "C_nonB") & (H["verdict"] == "homoplastic")]
    V = pd.read_csv(d / "variants_snps.Rtab", sep="\t", index_col=0)
    tips = list(V.columns)
    ena = pd.read_csv(args.ena, sep="\t", dtype=str).drop_duplicates("run_accession").set_index("run_accession")
    plat = pd.Series([ena.at[g, "instrument_model"] if g in ena.index else "assembly" for g in tips], index=tips)
    study = pd.Series([ena.at[g, "study_accession"] if g in ena.index else "assembly" for g in tips], index=tips)
    tr = pd.read_csv(args.traits, dtype=str).set_index("Name")
    reg = pd.Series("other", index=tips)
    for c in ("Region_Americas", "Region_Europe", "Region_Oceania"):
        reg[tr[c] == "1"] = c.replace("Region_", "")
    sub = list(pd.read_csv(d / "pheno" / "C_nonB__Region_Europe.tsv", sep="\t")["sample"])

    # Gubbins blocks: (start, end, set of taxa)
    blocks = []
    for line in open(args.gubbins):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        m = re.search(r'taxa="([^"]*)"', f[8])
        taxa = set(m.group(1).split()) if m else set()
        blocks.append((int(f[3]), int(f[4]), taxa))
    print("Gubbins recombination blocks: {}".format(len(blocks)))

    rows = []
    for v in H["variant"].unique():
        pos = int(v.split("_")[-3])
        x = V.loc[v, sub]
        if x.mean() > 0.5:
            x = 1 - x
        car = set(x.index[x == 1])
        inb = [(s, e, t) for s, e, t in blocks if s <= pos <= e]
        in_rec = set().union(*[t for _, _, t in inb]) if inb else set()
        rows.append({"variant": v, "pos": pos, "carriers": len(car),
                     "pattern": hash(frozenset(car)),
                     "gubbins_blocks_here": len(inb),
                     "carriers_in_recomb": len(car & in_rec),
                     "frac_carriers_in_recomb": round(len(car & in_rec) / len(car), 3) if car else 0,
                     "platform_mix": "; ".join("{}:{}".format(k, n) for k, n in
                                               plat[list(car)].value_counts().head(4).items()),
                     "region_mix": "; ".join("{}:{}".format(k, n) for k, n in
                                             reg[list(car)].value_counts().items()),
                     "n_studies": int(study[list(car)].nunique()),
                     "locus": H.loc[H["variant"] == v, "locus"].iloc[0]})
    R = pd.DataFrame(rows)
    # identical carrier sets
    R["pattern_group"] = R.groupby("pattern").ngroup() + 1
    R["sites_in_pattern"] = R.groupby("pattern")["variant"].transform("count")
    R = R.drop(columns="pattern").sort_values(["pattern_group", "pos"])
    R.to_csv(d / "homoplasy_check.tsv", sep="\t", index=False)

    # background: platform mix of the whole C_nonB subset, by region
    print("\nC_nonB background, platform x region:")
    bg = pd.crosstab(plat[sub], reg[sub])
    print(bg.to_string())
    pd.set_option("display.width", 260)
    pd.set_option("display.max_colwidth", 80)
    print("\nhomoplastic C_nonB SNPs: {} sites in {} distinct carrier patterns".format(
        len(R), R["pattern_group"].nunique()))
    print(R[["pattern_group", "sites_in_pattern", "variant", "carriers", "gubbins_blocks_here",
             "frac_carriers_in_recomb", "n_studies", "region_mix", "platform_mix",
             "locus"]].to_string(index=False))
    print("\nwritten: {}/homoplasy_check.tsv".format(d))


if __name__ == "__main__":
    main()
