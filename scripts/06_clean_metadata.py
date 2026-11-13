#!/usr/bin/env python3
"""
06_clean_metadata.py

Normalise the isolate metadata table and recompute the descriptive statistics
that go into the manuscript.

WHY THIS IS NEEDED
------------------
The published Supplementary Table S1 still carries the defects the original
sanitisation step was supposed to remove:

  * trailing whitespace splitting one category into two
        'Ileum' (75) vs 'Ileum ' (1)
        'Lymph Node' (92) vs 'Lymph Node ' (1)
        'Illumina HiSeq 2500' (318) vs 'Illumina HiSeq 2500 ' (6)
  * non-breaking spaces (U+00A0) at the end of 17 study_accession values
  * misspelled categories
        'Camle'      -> 'Camel'
        'Enviroment' -> 'Environment'
  * one numeric value in an otherwise textual column (454 -> '454')

Every change is recorded in a change log. Nothing is guessed: values this
script cannot resolve unambiguously are reported for a human decision rather
than silently recoded.

Outputs
-------
    01_metadata/metadata_clean.tsv        cleaned table, one row per isolate
    01_metadata/metadata_changelog.tsv    every value that was altered
    01_metadata/metadata_summary.txt      recomputed descriptive statistics

Usage
-----
    python 06_clean_metadata.py --input <path to Supplementary_Table_S1.xlsx>
    python 06_clean_metadata.py --input S1.xlsx --sheet Filter_Data
"""

import argparse
import re
import sys
import unicodedata
from pathlib import Path

import pandas as pd

# --- Unambiguous spelling corrections. Keys are compared after whitespace
# --- normalisation. Anything not listed here is left exactly as found.
SPELLING = {
    "host": {
        "Camle": "Camel",
        "Enviroment": "Environment",
    },
    "isolation_source": {},
    "country": {},
    "sequencing_model": {},
}

# Values that cannot be resolved without domain input. Reported, never recoded.
NEEDS_DECISION = {
    "host": ["Microbe"],
}

TEXT_COLS = ["country", "host", "isolation_source", "study_accession",
             "sequencing_model"]


def normalise_text(value):
    """Strip whitespace incl. NBSP, collapse internal runs, keep case."""
    if pd.isna(value):
        return pd.NA
    s = str(value)
    s = unicodedata.normalize("NFKC", s)      # turns U+00A0 into a plain space
    s = s.replace("\u00a0", " ").replace("\r", " ").replace("\n", " ")
    s = re.sub(r"\s+", " ", s).strip()
    return s if s else pd.NA


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="Supplementary_Table_S1.xlsx")
    ap.add_argument("--sheet", default="Filter_Data")
    ap.add_argument("--outdir", default="/cfs/earth/scratch/tabasran/MAP/01_metadata")
    args = ap.parse_args()

    src = Path(args.input)
    if not src.exists():
        sys.exit(f"[!] not found: {src}")
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_excel(src, sheet_name=args.sheet)
    df = df.loc[:, ~df.columns.astype(str).str.startswith("Unnamed")]
    df.columns = [str(c).strip() for c in df.columns]

    # The source spells this column 'Sequencign Model'; fix the header too.
    df = df.rename(columns={"Sequencign Model": "sequencing_model",
                            "Sequencing Model": "sequencing_model"})

    n_in = len(df)
    changes = []

    for col in TEXT_COLS:
        if col not in df.columns:
            print(f"[i] column absent, skipped: {col}")
            continue
        before = df[col].copy()
        df[col] = df[col].map(normalise_text)

        mapping = SPELLING.get(col, {})
        if mapping:
            df[col] = df[col].map(lambda v: mapping.get(v, v) if pd.notna(v) else v)

        for i in df.index:
            b, a = before.at[i], df[col].at[i]
            if pd.isna(b) and pd.isna(a):
                continue
            if str(b) != str(a):
                changes.append({
                    "accession": df.at[i, "accession"],
                    "column": col,
                    "before": repr(b),
                    "after": repr(a),
                })

    # ---------------------------------------------------------------- checks
    dup = df["accession"].duplicated().sum()
    if dup:
        print(f"[!] {dup} duplicated accessions — not removed, listed below:")
        print(df[df["accession"].duplicated(keep=False)]["accession"].to_string())

    flagged = []
    for col, vals in NEEDS_DECISION.items():
        if col not in df.columns:
            continue
        for v in vals:
            hits = df[df[col] == v]
            if len(hits):
                flagged.append((col, v, len(hits), list(hits["accession"])))

    # --------------------------------------------------------------- outputs
    clean_path = outdir / "metadata_clean.tsv"
    df.to_csv(clean_path, sep="\t", index=False)

    chg = pd.DataFrame(changes,
                       columns=["accession", "column", "before", "after"])
    chg_path = outdir / "metadata_changelog.tsv"
    chg.to_csv(chg_path, sep="\t", index=False)

    lines = []
    add = lines.append
    add("=" * 66)
    add(f" source            : {src.name}  [sheet: {args.sheet}]")
    add(f" isolates          : {n_in}")
    add(f" duplicated ids    : {dup}")
    add(f" values corrected  : {len(chg)}")
    add("=" * 66)

    if len(chg):
        add("")
        add("--- corrections by column ---")
        for col, grp in chg.groupby("column"):
            add(f"  {col}: {len(grp)}")
            for (b, a), g in grp.groupby(["before", "after"]):
                add(f"      {b} -> {a}   (n={len(g)})")

    if flagged:
        add("")
        add("--- values needing a decision (NOT recoded) ---")
        for col, v, n, accs in flagged:
            add(f"  {col} = {v!r}  n={n}")
            add(f"      {', '.join(map(str, accs))}")

    for col in TEXT_COLS:
        if col not in df.columns:
            continue
        vc = df[col].value_counts(dropna=False)
        pct = (vc / len(df) * 100).round(2)
        add("")
        add(f"--- {col}  ({df[col].nunique(dropna=True)} categories, "
            f"{int(df[col].isna().sum())} missing) ---")
        for k, n in vc.items():
            label = "missing" if pd.isna(k) else str(k)
            add(f"  {label:<34} {n:>5}  {pct[k] if pd.notna(k) else pct.get(k, 0):>6}%")

    summary = "\n".join(lines)
    sum_path = outdir / "metadata_summary.txt"
    sum_path.write_text(summary, encoding="utf-8")

    print(summary)
    print()
    print(f"clean table : {clean_path}")
    print(f"change log  : {chg_path}")
    print(f"summary     : {sum_path}")
    print()
    print("The percentages above are what the manuscript should report. Replace")
    print("the figures currently in the cohort-description section with these.")


if __name__ == "__main__":
    main()
