#!/usr/bin/env python3
"""
07_build_traits.py

Build the Scoary trait matrix from the cleaned metadata, with every recoding
rule written down in this file rather than applied by hand.

WHY THE PREVIOUS MATRIX IS BEING REPLACED
-----------------------------------------
Reverse-engineering the old trait.csv against the metadata showed four defects
that would each distort a pan-GWAS on their own:

  1. 406 isolates carried no isolation-source label at all, including 294 whose
     metadata says 'Feces' (21% of all faecal isolates), 8 'Farm Environment'
     and 2 'Tissue & Gut'. Conversely 12 isolates with NO recorded source were
     given one.
  2. ProductTypeDiary matched hostCattle in 1755 of 1780 rows (98.6%); the only
     difference was 25 Goat isolates. It was a relabelling of host, not an
     independent variable, so the same genes could surface twice under two
     names.
  3. hostEnviroment and Source_Environment covered the identical 79 isolates —
     one trait under two names.
  4. Oryx, an actual ruminant, was excluded from hostOther_Ruminants, while
     Goat, Deer and Camel were included. Egypt and Saudi Arabia fell into no
     region at all.

Every isolate now receives a label from each axis whenever its metadata
supports one, and no label when it does not. Nothing is inferred.

Outputs
-------
    06_pangwas/trait.csv                Scoary input, one row per isolate
    01_metadata/trait_rules.txt         the mapping, the counts, the exclusions
    01_metadata/trait_vs_old.tsv        per-isolate diff against the old matrix

COHORT
------
Only genomes in 01_metadata/cohort_final.txt are written, named exactly as
in the Panaroo gene_presence_absence.Rtab header (Scoary matches by name).
Metadata IDs without a version suffix (GCA_000835225) are matched to cohort
IDs with or without one. Any cohort genome without metadata stops the run.
If 05_phylogeny/typing/lineage_assignments.tsv exists, every trait is
cross-tabulated against lineage (S/I, S/III, C/B, C/non-B) to show how far
it is confounded with population structure.

Usage
-----
    python 07_build_traits.py
    python 07_build_traits.py --old-traits 00_archive/panaroo_run2/panaroo_Rana/trait.csv
    python 07_build_traits.py --min-n 20
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/cfs/earth/scratch/tabasran/MAP")

# ---------------------------------------------------------------------------
# MAPPINGS. Anything not listed receives no label on that axis.
# ---------------------------------------------------------------------------

REGION = {
    "Region_Americas": ["USA", "Canada", "Argentina", "Chile"],
    "Region_Europe":   ["Ireland", "France", "Italy", "Germany",
                        "Portugal", "Spain", "Russia"],
    "Region_Oceania":  ["Australia", "New Zealand"],
    "Region_Asia":     ["China", "South Korea", "India", "Japan"],
    # Egypt and Saudi Arabia were unlabelled in the old matrix. They are given
    # their own region here; with n=4 it will fall below --min-n and be dropped
    # with a note, which is honest, unlike silently leaving them out.
    "Region_Africa_MiddleEast": ["Egypt", "Saudi Arabia"],
}

HOST = {
    "hostCattle":          ["Cattle"],
    "hostSheep":           ["Sheep"],
    "hostHuman":           ["Human"],
    # Oryx added: it is a ruminant, and its exclusion from the old matrix was
    # inconsistent with Goat, Deer and Camel being included.
    "hostOther_Ruminants": ["Goat", "Deer", "Camel", "Oryx", "Bison"],
    # 'Environment' is not a host; those isolates are captured by
    # Source_Environment instead. 'Pig' is not a ruminant and is the only
    # member of its category. 'Microbe' means the submitter recorded no host
    # (checked in SRA for SRR14863282 / SRR14863283, PRJNA739165).
}

SOURCE = {
    "Source_Invasive":    ["Lymph Node", "Ileum", "Tissue & Gut"],
    "Source_Excretory":   ["Feces", "Human Breast Milk"],
    "Source_Environment": ["Farm Environment"],
}

# Deliberately NOT built:
#   ProductTypeDiary / ProductTypeNonDiary — collinear with hostCattle.
#   hostEnviroment — identical to Source_Environment (same 79 isolates).

AXES = [("country", REGION), ("host", HOST), ("isolation_source", SOURCE)]

# Values meaning "not known" on an axis. Isolates with these values are written
# as NA for every trait of that axis, so that Scoary does not count them as
# negatives (e.g. an isolate with no recorded source is NOT a non-invasive one).
# 'Environment' and 'Microbe' on the host axis record no host animal.
# Any other unmapped value (e.g. Pig, or an Asian country) is a KNOWN value
# outside the trait and stays a negative (0).
UNKNOWN = {"country": {""}, "host": {"", "Environment", "Microbe"},
           "isolation_source": {""}}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--metadata", default=str(ROOT / "01_metadata/metadata_clean.tsv"))
    ap.add_argument("--old-traits", default=None,
                    help="previous trait.csv, for the comparison report")
    ap.add_argument("--outdir", default=str(ROOT / "06_pangwas"))
    ap.add_argument("--metadir", default=str(ROOT / "01_metadata"))
    ap.add_argument("--cohort", default=str(ROOT / "01_metadata/cohort_final.txt"))
    ap.add_argument("--rtab", default=str(ROOT / "04_pangenome/panaroo/gene_presence_absence.Rtab"))
    ap.add_argument("--lineages", default=str(ROOT / "05_phylogeny/typing/lineage_assignments.tsv"))
    ap.add_argument("--min-n", type=int, default=20,
                    help="traits with fewer positives than this are dropped")
    args = ap.parse_args()

    meta_path = Path(args.metadata)
    if not meta_path.exists():
        sys.exit(f"[!] not found: {meta_path}\n    Run 06_clean_metadata.py first.")

    outdir, metadir = Path(args.outdir), Path(args.metadir)
    outdir.mkdir(parents=True, exist_ok=True)
    metadir.mkdir(parents=True, exist_ok=True)

    md = pd.read_csv(meta_path, sep="\t", dtype=str)
    if "accession" not in md.columns:
        sys.exit("[!] metadata has no 'accession' column")
    md["accession"] = md["accession"].str.strip()
    if md["accession"].duplicated().any():
        sys.exit(f"[!] duplicated accessions in metadata: "
                 f"{sorted(md.loc[md['accession'].duplicated(), 'accession'])[:10]}")
    md = md.set_index("accession")

    # ---- restrict to the locked cohort, named as in the Panaroo matrix
    cohort = [l.strip() for l in open(args.cohort) if l.strip()]
    base = lambda x: x.rsplit(".", 1)[0] if x.rsplit(".", 1)[-1].isdigit() else x
    by_base = {}
    for a in md.index:
        by_base.setdefault(base(a), []).append(a)
    rows, unmatched, ambiguous = {}, [], []
    for g in cohort:
        if g in md.index:
            rows[g] = g
        elif len(by_base.get(base(g), [])) == 1:
            rows[g] = by_base[base(g)][0]
        elif len(by_base.get(base(g), [])) > 1:
            ambiguous.append(g)
        else:
            unmatched.append(g)
    if unmatched or ambiguous:
        sys.exit(f"[!] cohort genomes without a unique metadata row: "
                 f"unmatched {unmatched[:20]} ambiguous {ambiguous[:20]}")
    n_meta_all = len(md)
    md = md.loc[[rows[g] for g in cohort]].copy()
    md.index = pd.Index(cohort, name="accession")

    rtab = Path(args.rtab)
    if rtab.exists():
        names = open(rtab).readline().rstrip("\n").split("\t")[1:]
        if set(names) != set(cohort):
            sys.exit(f"[!] Rtab genomes != cohort: only in Rtab "
                     f"{sorted(set(names) - set(cohort))[:10]}, only in cohort "
                     f"{sorted(set(cohort) - set(names))[:10]}")
    else:
        print(f"[i] {rtab} not found — name check against Panaroo skipped")

    lines, add = [], lambda s: lines.append(s)
    add("=" * 70)
    add(f" metadata : {meta_path}")
    add(f" metadata rows   : {n_meta_all}")
    add(f" cohort          : {args.cohort}")
    add(f" isolates        : {len(md)}  (cohort, named as in {Path(args.rtab).name})")
    add(f" min-n    : {args.min_n}")
    add("=" * 70)

    traits = pd.DataFrame(index=md.index)
    for col, mapping in AXES:
        if col not in md.columns:
            add(f"\n[!] column absent, axis skipped: {col}")
            continue
        values = md[col].fillna("")
        assigned = pd.Series(0, index=md.index)

        add(f"\n--- axis: {col} ---")
        for trait, members in mapping.items():
            hit = values.isin(members).astype(int)
            traits[trait] = hit.astype(float).where(~values.isin(UNKNOWN[col]))
            assigned += hit
            add(f"  {trait:<26} n={int(hit.sum()):<5} <- {', '.join(members)}")

        # Anything on this axis that received no label
        unl = values[assigned == 0]
        add(f"  unlabelled on this axis: {len(unl)}")
        if len(unl):
            vc = unl.replace("", "(missing)").value_counts()
            for k, n in vc.items():
                add(f"      {k:<24} {n}")
        add(f"  NA (unknown: {', '.join(sorted(v or '(missing)' for v in UNKNOWN[col]))}): "
            f"{int(values.isin(UNKNOWN[col]).sum())}")
        multi = int((assigned > 1).sum())
        if multi:
            add(f"  [!] {multi} isolates got more than one label on this axis")

    # ------------------------------------------------------- collinearity
    add("\n--- collinearity check (identical or near-identical traits) ---")
    cols = list(traits.columns)
    flagged = []
    for i in range(len(cols)):
        for j in range(i + 1, len(cols)):
            a, b = traits[cols[i]].fillna(0).values, traits[cols[j]].fillna(0).values
            inter = int(np.sum((a == 1) & (b == 1)))
            union = int(np.sum((a == 1) | (b == 1)))
            if union == 0:
                continue
            jac = inter / union
            if jac >= 0.90:
                flagged.append((cols[i], cols[j], jac, inter, union))
    if flagged:
        for a, b, jac, inter, union in flagged:
            add(f"  [!] {a} vs {b}: Jaccard={jac:.3f} ({inter}/{union})")
        add("      Traits this similar should not both enter the analysis.")
    else:
        add("  none above Jaccard 0.90 — OK")

    # ------------------------------------------------------- size filter
    add(f"\n--- traits below min-n ({args.min_n}) ---")
    npos = {c: int((traits[c] == 1).sum()) for c in traits.columns}
    nneg = {c: int((traits[c] == 0).sum()) for c in traits.columns}
    dropped = [c for c in traits.columns if min(npos[c], nneg[c]) < args.min_n]
    if dropped:
        for c in dropped:
            add(f"  dropped {c} (pos={npos[c]}, neg={nneg[c]}) — too few for a "
                f"meaningful association test")
        traits = traits.drop(columns=dropped)
    else:
        add("  none")

    # ------------------------------------------------ lineage confounding
    lp = Path(args.lineages)
    if lp.exists():
        lin = pd.read_csv(lp, sep="\t", dtype=str).set_index("genome")
        lab = (lin["lineage"] + "/" + lin["sublineage"]).replace(
            {"S/non-III": "S/I(non-III)"}).reindex(traits.index).fillna("NA")
        order = sorted(lab.unique())
        add("\n--- traits x lineage (positives per lineage; lineage sizes in header) ---")
        add("  " + f"{'trait':<26}" + "".join(f"{o:>14}" for o in order))
        add("  " + f"{'(lineage n)':<26}" + "".join(f"{int((lab == o).sum()):>14}" for o in order))
        for c in traits.columns:
            add("  " + f"{c:<26}" + "".join(
                f"{int((traits.loc[lab == o, c] == 1).sum()):>14}" for o in order))
        add("  A trait concentrated in one lineage cannot be separated from that\n"
            "  lineage by association alone; Scoary's pairwise test addresses this.")
    else:
        add(f"\n[i] {lp} not found — lineage cross-tab skipped")

    # ------------------------------------------------------------- write
    traits.index.name = "Name"
    trait_path = outdir / "trait.csv"
    traits.astype("Int64").to_csv(trait_path, na_rep="NA")

    add("\n" + "=" * 70)
    add(f" traits retained : {traits.shape[1]}")
    add(f" isolates        : {traits.shape[0]}")
    add("=" * 70)
    add("\n--- final trait counts ---")
    for c in traits.columns:
        n = int((traits[c] == 1).sum())
        n0 = int((traits[c] == 0).sum())
        na = int(traits[c].isna().sum())
        add(f"  {c:<26} pos {n:>5}   neg {n0:>5}   NA {na:>4}")

    # -------------------------------------------------------- comparison
    if args.old_traits:
        old_path = Path(args.old_traits)
        if not old_path.exists():
            add(f"\n[!] old trait file not found: {old_path}")
        else:
            old = pd.read_csv(old_path).set_index("Name")
            shared_i = traits.index.intersection(old.index)
            shared_c = [c for c in traits.columns if c in old.columns]
            add("\n" + "=" * 70)
            add(" COMPARISON WITH PREVIOUS MATRIX")
            add(f"  isolates: old {len(old)}, new {len(traits)}, "
                f"shared {len(shared_i)}")
            add(f"  traits removed : "
                f"{sorted(set(old.columns) - set(traits.columns))}")
            add(f"  traits added   : "
                f"{sorted(set(traits.columns) - set(old.columns))}")
            add("=" * 70)

            rows = []
            for c in shared_c:
                known = shared_i[traits.loc[shared_i, c].notna()]
                add(f"  {c:<26} NA in new (unknown on axis): {len(shared_i) - len(known)}")
                o = old.loc[known, c].astype(int)
                n = traits.loc[known, c].astype(int)
                gained = int(((o == 0) & (n == 1)).sum())
                lost = int(((o == 1) & (n == 0)).sum())
                add(f"  {c:<26} old={int(o.sum()):<5} new={int(n.sum()):<5} "
                    f"+{gained} / -{lost}")
                for acc in known[(o != n).values]:
                    rows.append({"accession": acc, "trait": c,
                                 "old": int(o.at[acc]), "new": int(n.at[acc])})
            diff = pd.DataFrame(rows, columns=["accession", "trait", "old", "new"])
            diff_path = metadir / "trait_vs_old.tsv"
            diff.to_csv(diff_path, sep="\t", index=False)
            add(f"\n  per-isolate differences: {len(diff)} -> {diff_path}")

    report = "\n".join(lines)
    (metadir / "trait_rules.txt").write_text(report, encoding="utf-8")
    print(report)
    print(f"\ntrait matrix : {trait_path}")
    print(f"rules/report : {metadir / 'trait_rules.txt'}")


if __name__ == "__main__":
    main()
