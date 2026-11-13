#!/usr/bin/env python3
"""
01_audit_gca.py

Inventory MAP/GCA/ and report every way it overlaps or duplicates the cohort,
BEFORE anything is added to the pan-genome set.

Three kinds of redundancy are known to be present:
  1. GenBank/RefSeq pairs of the same assembly  (GCA_000007865.1 / GCF_000007865.1)
  2. Versioned and unversioned copies           (GCA_000835225 / GCA_000835225.1)
  3. Strain-named directories alongside accessions (ATCC19698, DT3, JTC1281, ...)

There is also a large GCA_9091xxxxx / GCA_9055xxxxx block that is absent from
Final_gff. Those are ENA-side assemblies and are very likely the SAME isolates
already present as ERR6415xxx de novo assemblies. Adding both would put
duplicate genomes into the pan-genome, which would distort Scoary and the tree.
This script flags them; it does not decide.

Nothing is written except report files. Nothing is moved or deleted.

Usage:  python 01_audit_gca.py
"""

import re
from collections import defaultdict
from pathlib import Path

ROOT = Path("/cfs/earth/scratch/tabasran/MAP")
GCA_DIR = ROOT / "GCA"
GFF_DIR = ROOT / "03_genomes/gff"          # after 00_archive_and_restructure.sh
GFF_DIR_FALLBACK = ROOT / "Final_gff"      # before it
META = ROOT / "01_metadata"


def base_acc(name: str) -> str:
    """GCF_000007865.1 -> 000007865 ; GCA_000835225 -> 000835225 ; else name."""
    m = re.match(r"^GC[AF]_(\d+)", name)
    return m.group(1) if m else name


def main() -> None:
    META.mkdir(parents=True, exist_ok=True)

    if not GCA_DIR.is_dir():
        raise SystemExit(f"[!] {GCA_DIR} not found.")

    gff_dir = GFF_DIR if GFF_DIR.is_dir() else GFF_DIR_FALLBACK
    if not gff_dir.is_dir():
        raise SystemExit("[!] No GFF directory found.")

    entries = sorted(p.name for p in GCA_DIR.iterdir() if p.is_dir())

    # Which entries actually have a usable FASTA?
    with_fasta, without_fasta = [], []
    for e in entries:
        hits = list((GCA_DIR / e).glob("*.fasta")) + \
               list((GCA_DIR / e).glob("*.fna")) + \
               list((GCA_DIR / e).glob("*.fa"))
        (with_fasta if hits else without_fasta).append(e)

    cohort = {p.stem for p in gff_dir.glob("*.gff")}
    cohort_base = {base_acc(c) for c in cohort}

    groups: dict[str, list[str]] = defaultdict(list)
    for e in entries:
        groups[base_acc(e)].append(e)

    redundant = {k: v for k, v in groups.items() if len(v) > 1}

    in_cohort = [e for e in entries if e in cohort]
    same_assembly_diff_name = [
        e for e in entries if e not in cohort and base_acc(e) in cohort_base
    ]
    novel = [
        e for e in entries if e not in cohort and base_acc(e) not in cohort_base
    ]

    ena_block = [e for e in novel if re.match(r"^GC[AF]_9\d{8}", e)]
    other_novel = [e for e in novel if e not in ena_block]

    print("=" * 66)
    print(f" GCA/ sub-directories        : {len(entries)}")
    print(f"   with a FASTA              : {len(with_fasta)}")
    print(f"   WITHOUT a FASTA           : {len(without_fasta)}")
    print(f" GFFs in cohort ({gff_dir.name}) : {len(cohort)}")
    print("-" * 66)
    print(f" already in the cohort       : {len(in_cohort)}")
    print(f" same assembly, other naming : {len(same_assembly_diff_name)}")
    print(f" genuinely not in cohort     : {len(novel)}")
    print(f"   of which ENA 9xxxxxxxx    : {len(ena_block)}   <-- likely duplicates of ERR6415xxx")
    print(f"   of which other            : {len(other_novel)}")
    print("=" * 66)

    if redundant:
        print(f"\n--- {len(redundant)} accessions present under more than one name ---")
        for k in sorted(redundant)[:15]:
            print(f"   {k}: {', '.join(sorted(redundant[k]))}")
        if len(redundant) > 15:
            print(f"   ... and {len(redundant) - 15} more (see report)")

    if same_assembly_diff_name:
        print("\n--- in cohort under a different name (do NOT add again) ---")
        for e in sorted(same_assembly_diff_name)[:15]:
            match = sorted(c for c in cohort if base_acc(c) == base_acc(e))
            print(f"   {e}  ==  {', '.join(match)}")
        if len(same_assembly_diff_name) > 15:
            print(f"   ... and {len(same_assembly_diff_name) - 15} more")

    if other_novel:
        print(f"\n--- {len(other_novel)} non-ENA entries absent from the cohort ---")
        print("   " + ", ".join(sorted(other_novel)[:30]))

    if without_fasta:
        print(f"\n--- {len(without_fasta)} directories with no FASTA ---")
        print("   " + ", ".join(sorted(without_fasta)[:30]))

    # ----------------------------------------------------------------- write
    out = META / "gca_audit.tsv"
    with out.open("w") as fh:
        fh.write("entry\tbase_accession\thas_fasta\tstatus\n")
        for e in entries:
            if e in cohort:
                st = "in_cohort"
            elif base_acc(e) in cohort_base:
                st = "duplicate_of_cohort_member"
            elif e in ena_block:
                st = "not_in_cohort_ena_block"
            else:
                st = "not_in_cohort_other"
            fh.write(f"{e}\t{base_acc(e)}\t{int(e in with_fasta)}\t{st}\n")

    print(f"\nFull table: {out}")
    print("\nNEXT: decide whether the ENA block is added. If those isolates are")
    print("already in the cohort as ERR6415xxx de novo assemblies, adding them")
    print("would duplicate genomes. Cross-check a few accessions on ENA before")
    print("going further.")


if __name__ == "__main__":
    main()
