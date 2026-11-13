#!/usr/bin/env python3
"""
38_dates.py

Collection dates for the dated phylogeny (BactDating), with every rule
explicit and every excluded genome listed with its reason.

Sources
  runs (SRR/ERR/DRR)  ENA read_run field collection_date (01_metadata/ena_raw.tsv)
  NCBI assemblies     BioSample attribute 'collection date' / 'collection_date'
                      fetched from ENA (sample accession from
                      02_qc/integrity/ena_vs_file.tsv); needs internet ->
                      run on the login node
Parsing (decimal years; BactDating takes a date range for imprecise dates)
  YYYY-MM-DD -> exact day (lower = upper)
  YYYY-MM    -> that month   [first day, last day]
  YYYY       -> that year    [Y, Y+1)
  'A/B' ranges (YYYY/YYYY or full dates) -> [A, B]
  anything else (missing, not collected, restricted, ...) -> no date
Excluded from dating (kept in all other analyses), each with a reason:
  culture-collection / lab-passaged strains whose isolation date is not the
  sampling date: ATCC 19698, DSM 44135 (by genome ID below)
  possible mixed infections: >= --max-mixed mixed core-SNP sites from reads
  (08_reads/summary/genome_qc_reads.tsv)
  dates later than the run's first_public date (impossible) or before 1900

Outputs (09_dating/)
  dates.tsv            genome, raw, lower, upper, mid, precision, used, reason
  dates_report.txt     counts by precision, by lineage and region, range
"""

import argparse
import calendar
import re
import sys
import time
import urllib.request
import xml.etree.ElementTree as ET
from datetime import date
from pathlib import Path

import pandas as pd

ROOT = Path("/cfs/earth/scratch/tabasran/MAP")
LAB_STRAINS = {"GCA_000240525.2": "culture collection ATCC 19698",
               "GCA_013357385.1": "culture collection DSM 44135"}


def dec(y, m=1, d=1):
    start = date(y, 1, 1).toordinal()
    n = 366 if calendar.isleap(y) else 365
    return y + (date(y, m, d).toordinal() - start) / n


def parse(raw):
    """Return (lower, upper, precision) in decimal years or None."""
    if raw is None:
        return None
    s = str(raw).strip()
    if not s or s.lower() in ("nan", "none"):
        return None
    if "/" in s:
        a, b = s.split("/", 1)
        pa, pb = parse(a), parse(b)
        if pa and pb:
            return (pa[0], pb[1], "range")
        return None
    m = re.fullmatch(r"(\d{4})-(\d{2})-(\d{2})(?:T.*)?", s)
    if m:
        y, mo, d = map(int, m.groups())
        try:
            v = dec(y, mo, d)
        except ValueError:
            return None
        return (v, v, "day")
    m = re.fullmatch(r"(\d{4})-(\d{2})", s)
    if m:
        y, mo = map(int, m.groups())
        if not 1 <= mo <= 12:
            return None
        last = calendar.monthrange(y, mo)[1]
        return (dec(y, mo, 1), dec(y, mo, last), "month")
    m = re.fullmatch(r"(\d{4})", s)
    if m:
        y = int(m.group(1))
        return (float(y), y + 0.9999, "year")
    return None


def fetch_sample_date(acc, tries=3):
    url = "https://www.ebi.ac.uk/ena/browser/api/xml/" + acc
    for t in range(tries):
        try:
            x = urllib.request.urlopen(url, timeout=60).read()
            root = ET.fromstring(x)
            for a in root.iter("SAMPLE_ATTRIBUTE"):
                tag = (a.findtext("TAG") or "").lower().replace("_", " ")
                if tag == "collection date":
                    return a.findtext("VALUE")
            return None
        except Exception:
            time.sleep(2 * (t + 1))
    return "FETCH_ERROR"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", default=str(ROOT / "01_metadata/cohort_final.txt"))
    ap.add_argument("--ena", default=str(ROOT / "01_metadata/ena_raw.tsv"))
    ap.add_argument("--gca", default=str(ROOT / "02_qc/integrity/ena_vs_file.tsv"))
    ap.add_argument("--readqc", default=str(ROOT / "08_reads/summary/genome_qc_reads.tsv"))
    ap.add_argument("--lineages", default=str(ROOT / "05_phylogeny/typing/lineage_assignments.tsv"))
    ap.add_argument("--traits", default=str(ROOT / "06_pangwas/trait.csv"))
    ap.add_argument("--out", default=str(ROOT / "09_dating"))
    ap.add_argument("--max-mixed", type=int, default=20)
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    R = []

    def say(s=""):
        print(s, flush=True)
        R.append(s)

    cohort = [l.strip() for l in open(args.cohort) if l.strip()]
    ena = pd.read_csv(args.ena, sep="\t", dtype=str).drop_duplicates("run_accession").set_index("run_accession")
    gca = pd.read_csv(args.gca, sep="\t", dtype=str).set_index("genome") if Path(args.gca).exists() else None
    mix = pd.read_csv(args.readqc, sep="\t", index_col=0) if Path(args.readqc).exists() else None
    if mix is None:
        say("[!] {} not found — mixed-infection exclusion skipped".format(args.readqc))

    rows = []
    for g in cohort:
        raw, pub = None, None
        if g in ena.index:
            raw = ena.at[g, "collection_date"] if "collection_date" in ena.columns else None
            pub = ena.at[g, "first_public"] if "first_public" in ena.columns else None
        elif gca is not None and g in gca.index and isinstance(gca.at[g, "ena_sample"], str):
            raw = fetch_sample_date(gca.at[g, "ena_sample"])
            time.sleep(0.3)
        p = parse(raw)
        r = {"genome": g, "raw": raw, "lower": None, "upper": None, "mid": None,
             "precision": None, "used": False, "reason": ""}
        if p:
            r["lower"], r["upper"], r["precision"] = round(p[0], 4), round(p[1], 4), p[2]
            r["mid"] = round((p[0] + p[1]) / 2, 4)
        reasons = []
        if p is None:
            reasons.append("no parseable date ({})".format(raw))
        if g in LAB_STRAINS:
            reasons.append(LAB_STRAINS[g])
        if mix is not None and g in mix.index and mix.at[g, "mixed_sites"] >= args.max_mixed:
            reasons.append("possible mixed infection ({} mixed sites)".format(int(mix.at[g, "mixed_sites"])))
        if p:
            pp = parse(pub)
            if pp and p[0] > pp[1] + 0.01:
                reasons.append("collection after first_public ({})".format(pub))
            if p[0] < 1900:
                reasons.append("before 1900")
        r["used"] = not reasons
        r["reason"] = "; ".join(reasons)
        rows.append(r)

    D = pd.DataFrame(rows)
    D.to_csv(out / "dates.tsv", sep="\t", index=False)
    lin = pd.read_csv(args.lineages, sep="\t", dtype=str).set_index("genome")
    tr = pd.read_csv(args.traits, dtype=str).set_index("Name")
    reg = pd.Series("unknown", index=tr.index)
    for c in ("Region_Americas", "Region_Europe", "Region_Oceania"):
        reg[tr[c] == "1"] = c.replace("Region_", "")
    D["lineage"] = (lin["lineage"] + "/" + lin["sublineage"]).reindex(D["genome"]).values
    D["region"] = reg.reindex(D["genome"]).values

    say("genomes: {} | with a parseable date: {} | used for dating: {}".format(
        len(D), int(D["precision"].notna().sum()), int(D["used"].sum())))
    say("\nprecision of parseable dates:\n" + D["precision"].value_counts().to_string())
    say("\nexclusion reasons (a genome may have several):")
    for k, n in pd.Series([x for r in D["reason"] for x in r.split("; ") if x]).str.replace(
            r"\(.*\)", "(...)", regex=True).value_counts().items():
        say("  {:>5}  {}".format(n, k))
    U = D[D["used"]]
    say("\nused genomes: date range {:.1f} - {:.1f} (median {:.1f})".format(
        U["mid"].min(), U["mid"].max(), U["mid"].median()))
    say("\nused, by lineage:\n" + U.groupby("lineage")["mid"].agg(["count", "min", "median", "max"]).round(1).to_string())
    say("\nused, by region:\n" + U.groupby("region")["mid"].agg(["count", "min", "median", "max"]).round(1).to_string())
    say("\nused, by decade:\n" + (U["mid"] // 10 * 10).astype(int).value_counts().sort_index().to_string())
    fe = int((D["raw"] == "FETCH_ERROR").sum())
    if fe:
        say("\n[!] {} assembly dates could not be fetched (network) — rerun on the login node".format(fe))
    (out / "dates_report.txt").write_text("\n".join(R) + "\n")
    say("\nwritten: {}/dates.tsv, dates_report.txt".format(out))


if __name__ == "__main__":
    main()
