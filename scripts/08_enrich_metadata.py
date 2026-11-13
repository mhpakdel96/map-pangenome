#!/usr/bin/env python3
"""
08_enrich_metadata.py

Fetch the full BioSample-level attributes for every sequencing run in the
cohort from the ENA Portal API, and report what they would change.

This script is deliberately NON-DESTRUCTIVE. It writes a proposed table and a
field-by-field diff; it never overwrites metadata_clean.tsv. Nothing enters the
analysis until the diff has been read.

WHAT IT CAN FILL
----------------
    isolation_source   125 isolates currently missing
    country              2 isolates currently missing
    host                 2 isolates recorded as 'Microbe' (SRA lists no host)

WHAT IT CAN ADD
---------------
    collection_date    enables a temporal / dated-phylogeny analysis
    strain / isolate   the real strain name, for linking to earlier studies
    lat_lon            finer geography than a continental bucket
    serotype           if present, the classical MAP C-type / S-type split —
                       a far better balanced and more biologically meaningful
                       trait than Source_Excretory, which sits on 77% of the
                       cohort

RUN ON THE LOGIN NODE (needs internet).

Usage
-----
    python 08_enrich_metadata.py
    python 08_enrich_metadata.py --chunk 50 --sleep 1.0
    python 08_enrich_metadata.py --offline   # re-parse a previous download
"""

import argparse
import io
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path("/cfs/earth/scratch/tabasran/MAP")
ENA = "https://www.ebi.ac.uk/ena/portal/api/search"

FIELDS = [
    "run_accession", "sample_accession", "secondary_sample_accession",
    "study_accession", "experiment_accession",
    "scientific_name", "tax_id", "sub_species",
    "instrument_model", "library_layout", "library_strategy",
    "read_count", "base_count",
    "country", "location", "collection_date", "first_public",
    "host", "host_scientific_name", "host_status",
    "isolation_source", "environment_material",
    "strain", "isolate", "serotype", "sample_title", "sample_alias",
]


def fetch_chunk(accessions, timeout=120, tries=3, sleep=2.0):
    """POST one chunk of run accessions to the ENA Portal API -> DataFrame."""
    data = urllib.parse.urlencode({
        "result": "read_run",
        "includeAccessions": ",".join(accessions),
        "fields": ",".join(FIELDS),
        "format": "tsv",
        "limit": "0",
    }).encode()

    last = None
    for attempt in range(1, tries + 1):
        try:
            req = urllib.request.Request(
                ENA, data=data,
                headers={"Content-Type": "application/x-www-form-urlencoded",
                         "Accept": "text/plain"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                text = r.read().decode("utf-8", errors="replace")
            if not text.strip():
                return pd.DataFrame(columns=FIELDS)
            return pd.read_csv(io.StringIO(text), sep="\t", dtype=str)
        except Exception as exc:                      # network, 5xx, timeout
            last = exc
            if attempt < tries:
                time.sleep(sleep * attempt)
    print(f"  [!] chunk failed after {tries} tries: {last}")
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--metadata", default=str(ROOT / "01_metadata/metadata_clean.tsv"))
    ap.add_argument("--outdir", default=str(ROOT / "01_metadata"))
    ap.add_argument("--chunk", type=int, default=100)
    ap.add_argument("--sleep", type=float, default=0.5,
                    help="pause between chunks, to stay polite to the API")
    ap.add_argument("--offline", action="store_true",
                    help="skip the download and re-parse ena_raw.tsv")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    raw_path = outdir / "ena_raw.tsv"

    md = pd.read_csv(args.metadata, sep="\t", dtype=str)
    if "accession" not in md.columns:
        sys.exit("[!] metadata has no 'accession' column")

    runs = [a for a in md["accession"] if str(a)[:3] in ("SRR", "ERR", "DRR")]
    refs = [a for a in md["accession"] if str(a)[:3] not in ("SRR", "ERR", "DRR")]
    print(f"run accessions to query : {len(runs)}")
    print(f"assembly accessions     : {len(refs)}  (not in read_run; skipped)")

    # ------------------------------------------------------------ download
    if args.offline:
        if not raw_path.exists():
            sys.exit(f"[!] --offline given but {raw_path} does not exist")
        ena = pd.read_csv(raw_path, sep="\t", dtype=str)
        print(f"re-using {raw_path}: {len(ena)} rows")
    else:
        frames, failed = [], []
        chunks = [runs[i:i + args.chunk] for i in range(0, len(runs), args.chunk)]
        for i, ch in enumerate(chunks, 1):
            print(f"  chunk {i}/{len(chunks)} ({len(ch)} accessions)...", flush=True)
            df = fetch_chunk(ch)
            if df is None:
                failed.extend(ch)
            elif len(df):
                frames.append(df)
            time.sleep(args.sleep)

        if failed:
            print(f"[!] {len(failed)} accessions in failed chunks; retrying singly")
            for a in failed:
                df = fetch_chunk([a], tries=2)
                if df is not None and len(df):
                    frames.append(df)
                time.sleep(args.sleep)

        if not frames:
            sys.exit("[!] nothing retrieved. Check network access from this node.")
        ena = pd.concat(frames, ignore_index=True)
        ena.to_csv(raw_path, sep="\t", index=False)
        print(f"\nraw ENA table: {raw_path}  ({len(ena)} rows)")

    ena = ena.drop_duplicates(subset="run_accession", keep="first")
    got = set(ena["run_accession"])
    missing = [a for a in runs if a not in got]
    print(f"retrieved {len(got)} / {len(runs)} runs; {len(missing)} not found")
    if missing:
        (outdir / "ena_not_found.txt").write_text("\n".join(missing))

    # -------------------------------------------------------------- report
    lines, add = [], lambda s: lines.append(s)
    add("=" * 70)
    add(f" cohort runs      : {len(runs)}")
    add(f" retrieved        : {len(got)}")
    add(f" not found in ENA : {len(missing)}")
    add("=" * 70)

    add("\n--- coverage of each ENA field (non-empty / retrieved) ---")
    for f in FIELDS:
        if f not in ena.columns:
            add(f"  {f:<28} field absent from response")
            continue
        n = int((ena[f].fillna("").astype(str).str.strip() != "").sum())
        add(f"  {f:<28} {n:>5}  ({n / max(len(ena), 1) * 100:5.1f}%)")

    merged = md.merge(ena, left_on="accession", right_on="run_accession",
                      how="left", suffixes=("", "_ena"))

    # Gaps the enrichment could close
    add("\n--- gaps in the current metadata that ENA can fill ---")
    pairs = [("isolation_source", "isolation_source_ena"),
             ("country", "country_ena"),
             ("host", "host_ena")]
    for ours, theirs in pairs:
        if ours not in merged.columns or theirs not in merged.columns:
            continue
        blank = merged[ours].isna() | (merged[ours].fillna("") == "")
        fillable = blank & merged[theirs].notna() & (merged[theirs].fillna("") != "")
        add(f"  {ours:<20} missing {int(blank.sum()):>4}  -> "
            f"fillable from ENA {int(fillable.sum()):>4}")
        if int(fillable.sum()):
            vc = merged.loc[fillable, theirs].value_counts().head(12)
            for k, n in vc.items():
                add(f"        {k:<34} {n}")

    # Disagreements between our value and ENA's
    add("\n--- disagreements (both present, values differ) ---")
    diff_rows = []
    for ours, theirs in pairs:
        if ours not in merged.columns or theirs not in merged.columns:
            continue
        both = (merged[ours].fillna("") != "") & (merged[theirs].fillna("") != "")
        neq = both & (merged[ours].astype(str).str.strip().str.lower()
                      != merged[theirs].astype(str).str.strip().str.lower())
        add(f"  {ours:<20} {int(neq.sum())}")
        for _, r in merged[neq].iterrows():
            diff_rows.append({"accession": r["accession"], "field": ours,
                              "current": r[ours], "ena": r[theirs]})

    pd.DataFrame(diff_rows, columns=["accession", "field", "current", "ena"]) \
        .to_csv(outdir / "ena_disagreements.tsv", sep="\t", index=False)

    # New fields worth having
    add("\n--- new fields available ---")
    for f in ["collection_date", "strain", "isolate", "serotype",
              "sub_species", "location", "host_status"]:
        c = f if f in merged.columns else f"{f}_ena"
        if c not in merged.columns:
            continue
        n = int((merged[c].fillna("") != "").sum())
        add(f"  {f:<20} present for {n:>5} / {len(merged)}")
        if n and f in ("serotype", "sub_species", "host_status"):
            for k, v in merged[c].value_counts().head(8).items():
                add(f"        {k:<30} {v}")

    proposed = outdir / "metadata_enriched_proposed.tsv"
    merged.to_csv(proposed, sep="\t", index=False)

    report = "\n".join(lines)
    (outdir / "ena_enrichment_report.txt").write_text(report, encoding="utf-8")
    print("\n" + report)
    print(f"\nproposed table  : {proposed}")
    print(f"disagreements   : {outdir / 'ena_disagreements.tsv'}")
    print(f"report          : {outdir / 'ena_enrichment_report.txt'}")
    print("\nNothing was overwritten. Read the disagreements file before")
    print("deciding which ENA values to adopt.")


if __name__ == "__main__":
    main()
