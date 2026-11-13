#!/usr/bin/env python3
"""
09_apply_filter.py

Apply the assembly-quality thresholds to the cohort and lock the final genome
set. This is the step after which the pan-genome, phylogeny and pan-GWAS are
run, so any change to it later means re-running everything downstream.

Design
------
  * Thresholds are given explicitly on the command line and recorded in the
    log, so the Methods section can quote them exactly.
  * DRY RUN BY DEFAULT. Nothing is moved unless --apply is given.
  * Nothing is ever deleted. Failing genomes are MOVED (GFF and the FASTA
    symlink) into 03_genomes/excluded/step3_assembly_qc/, where they stay
    inspectable.
  * Protected genomes (default: the K-10 reference, GCA_000007865.1) must
    pass. If one fails, the script stops: a failing reference means the
    thresholds or the QC data are wrong, not the reference.

Usage
-----
    # 1) dry run — shows exactly what would be removed and why
    python 09_apply_filter.py --n50 10000 --completeness 95 --contamination 5

    # 2) once the list has been read and agreed
    python 09_apply_filter.py --n50 10000 --completeness 95 --contamination 5 --apply

Outputs
-------
    01_metadata/cohort_final.txt          IDs of every retained genome
    02_qc/excluded_step3_assembly_qc.tsv  every removed genome, metrics, reasons
    02_qc/filter_log.txt                  thresholds, counts, date
"""

import argparse
import datetime
import shutil
import sys
from pathlib import Path

import pandas as pd

ROOT = Path("/cfs/earth/scratch/tabasran/MAP")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--master", default=str(ROOT / "02_qc/genome_qc_master.tsv"))
    ap.add_argument("--n50", type=float, default=None, help="minimum N50")
    ap.add_argument("--contigs", type=float, default=None, help="maximum contigs")
    ap.add_argument("--completeness", type=float, default=None,
                    help="minimum CheckM2 completeness (%%)")
    ap.add_argument("--contamination", type=float, default=None,
                    help="maximum CheckM2 contamination (%%)")
    ap.add_argument("--exclude", default=None,
                    help="TSV of genomes to exclude regardless of QC: two "
                         "columns, genome<TAB>reason (e.g. not_MAP). These go "
                         "to excluded/step4_manual/, separate from QC failures")
    ap.add_argument("--protect", nargs="*", default=["GCA_000007865.1"],
                    help="genomes that must pass (default: K-10 reference)")
    ap.add_argument("--apply", action="store_true",
                    help="actually move failing genomes (default: dry run)")
    ap.add_argument("--root", default=str(ROOT))
    args = ap.parse_args()

    root = Path(args.root)
    gff_dir = root / "03_genomes/gff"
    fa_dir = root / "03_genomes/fasta"
    ex_dir = root / "03_genomes/excluded/step3_assembly_qc"
    ex_dir_manual = root / "03_genomes/excluded/step4_manual"
    qc_dir = root / "02_qc"
    meta_dir = root / "01_metadata"

    if all(v is None for v in (args.n50, args.contigs,
                               args.completeness, args.contamination)) \
            and args.exclude is None:
        sys.exit("[!] give at least one threshold (--n50, --contigs, "
                 "--completeness, --contamination) or --exclude")

    mpath = Path(args.master)
    if not mpath.exists():
        sys.exit(f"[!] {mpath} not found. Run 05_filter_decision.py first.")
    df = pd.read_csv(mpath, sep="\t")

    need = []
    if args.n50 is not None:
        need.append("n50")
    if args.contigs is not None:
        need.append("contigs")
    if args.completeness is not None:
        need.append("completeness")
    if args.contamination is not None:
        need.append("contamination")
    missing_cols = [c for c in need if c not in df.columns]
    if missing_cols:
        sys.exit(f"[!] master table lacks columns: {missing_cols}")

    # A genome with a missing metric cannot be judged; it is reported, not
    # silently passed or failed.
    na = df[need].isna().any(axis=1)
    if na.any():
        print(f"[!] {int(na.sum())} genomes have missing QC values:")
        print(df.loc[na, ["genome"] + need].to_string(index=False))
        sys.exit("    Resolve these before filtering.")

    reasons = {}
    if args.n50 is not None:
        reasons["n50"] = df["n50"] < args.n50
    if args.contigs is not None:
        reasons["contigs"] = df["contigs"] > args.contigs
    if args.completeness is not None:
        reasons["completeness"] = df["completeness"] < args.completeness
    if args.contamination is not None:
        reasons["contamination"] = df["contamination"] > args.contamination

    manual_reason = pd.Series("", index=df.index)
    if args.exclude is not None:
        xpath = Path(args.exclude)
        if not xpath.exists():
            sys.exit(f"[!] exclude file not found: {xpath}")
        ex = pd.read_csv(xpath, sep="\t", header=None, names=["genome", "reason"],
                         dtype=str, comment="#").dropna(subset=["genome"])
        ex["genome"] = ex["genome"].str.strip()
        if ex["reason"].isna().any() or (ex["reason"].str.strip() == "").any():
            sys.exit("[!] every line in the exclude file needs a reason (column 2)")
        unknown = sorted(set(ex["genome"]) - set(df["genome"]))
        if unknown:
            sys.exit(f"[!] exclude file lists genomes not in the QC table: {unknown}")
        rmap = dict(zip(ex["genome"], ex["reason"].str.strip()))
        manual_reason = df["genome"].map(rmap).fillna("")
        reasons["manual"] = manual_reason != ""

    fail = pd.Series(False, index=df.index)
    for r in reasons.values():
        fail |= r
    df["fail_reason"] = [
        ",".join((f"manual:{manual_reason.iloc[i]}" if k == "manual" else k)
                 for k, v in reasons.items() if bool(v.iloc[i]))
        for i in range(len(df))
    ]
    df["passes"] = ~fail
    df["is_manual"] = manual_reason != ""

    # ------------------------------------------------------ protected check
    bad_protect = [g for g in args.protect
                   if g in set(df["genome"]) and not
                   bool(df.loc[df["genome"] == g, "passes"].iloc[0])]
    absent_protect = [g for g in args.protect if g not in set(df["genome"])]
    if absent_protect:
        print(f"[i] protected genomes not in cohort: {absent_protect}")
    if bad_protect:
        print("[!] PROTECTED GENOME(S) FAIL THE THRESHOLDS — stopping.")
        print(df.loc[df["genome"].isin(bad_protect),
                     ["genome"] + need + ["fail_reason"]].to_string(index=False))
        sys.exit(1)

    # ------------------------------------------------------------- report
    removed = df[fail].copy()
    kept = df[~fail].copy()
    stamp = datetime.datetime.now().isoformat(timespec="seconds")

    lines, add = [], lambda s: lines.append(s)
    add("=" * 66)
    add(f" date       : {stamp}")
    add(f" mode       : {'APPLY' if args.apply else 'DRY RUN'}")
    add(" thresholds :")
    for k in ("n50", "contigs", "completeness", "contamination"):
        v = getattr(args, k)
        if v is not None:
            op = {"n50": ">=", "contigs": "<=",
                  "completeness": ">=", "contamination": "<="}[k]
            add(f"     {k:<14} {op} {v:g}")
    if args.exclude is not None:
        add(f"     exclude list   {args.exclude} ({int(df['is_manual'].sum())} genomes)")
    add(f" input      : {len(df)}")
    add(f" removed    : {len(removed)}")
    add(f" retained   : {len(kept)}")
    add("=" * 66)

    add("\n--- removed, by reason (a genome may fail several) ---")
    for k, v in reasons.items():
        add(f"   {k:<14} {int(v.sum())}")
    add("\n--- removed, by exact reason combination ---")
    for combo, n in removed["fail_reason"].value_counts().items():
        add(f"   {combo:<40} {n}")
    if "source" in df.columns:
        add("\n--- by source ---")
        add(f"   {'source':<12} {'input':>6} {'removed':>8} {'retained':>9}")
        for s, g in df.groupby("source"):
            r = int((~g["passes"]).sum())
            add(f"   {s:<12} {len(g):>6} {r:>8} {len(g) - r:>9}")

    show = ["genome"] + [c for c in ("source", "contigs", "n50", "length",
                                     "completeness", "contamination")
                         if c in df.columns] + ["fail_reason"]
    add("\n--- removed genomes ---")
    add(removed.sort_values("fail_reason")[show].to_string(index=False)
        if len(removed) else "   none")

    report = "\n".join(lines)
    print(report)

    if not args.apply:
        print("\nDRY RUN — nothing was moved. Re-run with --apply once the list "
              "above has been reviewed.")
        return

    # --------------------------------------------------------------- apply
    for d in (ex_dir, ex_dir_manual):
        (d / "gff").mkdir(parents=True, exist_ok=True)
        (d / "fasta").mkdir(parents=True, exist_ok=True)

    moved, already, problems = 0, 0, []
    for g, is_manual in zip(removed["genome"], removed["is_manual"]):
        # QC failures keep going to step3 even if also listed manually, so a
        # re-run never moves a genome that is already in its original place.
        qc_fail = any(bool(v[df["genome"] == g].iloc[0])
                      for k, v in reasons.items() if k != "manual")
        tgt = ex_dir if (qc_fail or not is_manual) else ex_dir_manual
        for src, dst in ((gff_dir / f"{g}.gff", tgt / "gff" / f"{g}.gff"),
                         (fa_dir / f"{g}.fasta", tgt / "fasta" / f"{g}.fasta")):
            if src.exists() or src.is_symlink():
                if dst.exists() or dst.is_symlink():
                    problems.append(f"destination exists: {dst}")
                    continue
                shutil.move(str(src), str(dst))
                moved += 1
            elif dst.exists() or dst.is_symlink():
                already += 1
            else:
                problems.append(f"not found anywhere: {src.name}")

    meta_dir.mkdir(parents=True, exist_ok=True)
    (meta_dir / "cohort_final.txt").write_text(
        "\n".join(sorted(kept["genome"])) + "\n")
    removed[show].to_csv(qc_dir / "excluded_step3_assembly_qc.tsv",
                         sep="\t", index=False)
    (qc_dir / "filter_log.txt").write_text(report + "\n", encoding="utf-8")

    n_gff = len(list(gff_dir.glob("*.gff")))
    n_fa = len(list(fa_dir.glob("*.fasta")))

    print("\n" + "=" * 66)
    print(f" files moved             : {moved}")
    print(f" already in excluded/    : {already}")
    print(f" GFFs now in gff/        : {n_gff}")
    print(f" FASTAs now in fasta/    : {n_fa}")
    print(f" expected (retained)     : {len(kept)}")
    print("=" * 66)
    for p in problems:
        print(f" [!] {p}")

    if n_gff == n_fa == len(kept) and not problems:
        print("\n OK — final cohort locked.")
        print(f"   {meta_dir / 'cohort_final.txt'}")
        print(f"   {qc_dir / 'excluded_step3_assembly_qc.tsv'}")
        print(f"   {qc_dir / 'filter_log.txt'}")
        sys.exit(0)
    print("\n NOT OK — counts do not match. Do not start the pan-genome run.")
    sys.exit(1)


if __name__ == "__main__":
    main()
