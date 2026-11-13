#!/usr/bin/env python3
"""
17_integrity_check.py

Check that every genome FILE in the cohort is what its name says it is.
Prompted by GCA_000218155.2: its file holds the Env 210 assembly (WGS AGAO01,
byte-identical to GCA_000240465.2), while ENA says GCA_000218155.2 is strain
Pt164 (WGS AFPI). QC, CheckM2 and the species check cannot catch this.

Two independent modes (Python 3.5 compatible, standard library only):

  hash   for every cohort genome: number of contigs, length, and a content
         hash (md5 of the upper-cased contig sequences, sorted, so contig
         order and names do not matter). Genomes with the same hash are the
         same assembly under two names. Reads ~9 GB: run it as a SLURM job.

  ena    for every cohort genome named GCA_*: fetch the ENA assembly record
         and compare it with the accessions written in the FASTA headers
         (WGS prefix for drafts, chromosome accession for complete genomes).
         Needs internet: run it on the login node. ~1 s per genome.

Outputs (02_qc/integrity/)
  genome_hashes.tsv        genome, contigs, bp, hash, header accessions
  duplicate_hashes.tsv     groups of genomes with identical content
  ena_vs_file.tsv          genome, ENA strain/title/WGS/chromosomes/sample,
                           accessions in the file, verdict

Usage
  python 17_integrity_check.py hash --threads 16
  python 17_integrity_check.py ena
"""

import argparse
import hashlib
import re
import sys
import time
import urllib.request
import xml.etree.ElementTree as ET
from collections import defaultdict
from multiprocessing import Pool
from pathlib import Path

ROOT = Path("/cfs/earth/scratch/tabasran/MAP")
FASTA = ROOT / "03_genomes/fasta"
COHORT = ROOT / "01_metadata/cohort_final.txt"
OUT = ROOT / "02_qc/integrity"

# accession-like tokens in FASTA headers: WGS contigs (AGAO01000001,
# JAHWGM010000001, NZ_LKUS01000001) and complete records (CP010113, AE016958)
RE_WGS = re.compile(r"(?:NZ_)?([A-Z]{4,6})(\d{2})\d{6,9}(?:\.\d+)?")
RE_NUC = re.compile(r"(?:NZ_)?([A-Z]{1,2}\d{5,8})(?:\.\d+)?")


def header_accessions(headers):
    wgs, nuc = set(), set()
    for h in headers:
        for tok in re.split(r"[|\s]+", h):
            m = RE_WGS.fullmatch(tok)
            if m:
                wgs.add(m.group(1))
                continue
            m = RE_NUC.fullmatch(tok)
            if m:
                nuc.add(m.group(1))
    return wgs, nuc


def one_genome(g):
    path = FASTA / (g + ".fasta")
    seqs, headers, cur = [], [], []
    with open(str(path)) as fh:
        for line in fh:
            if line.startswith(">"):
                if headers:
                    seqs.append("".join(cur).upper())
                headers.append(line[1:].strip())
                cur = []
            else:
                cur.append(line.strip())
    if headers:
        seqs.append("".join(cur).upper())
    h = hashlib.md5()
    for s in sorted(seqs):
        h.update(s.encode())
        h.update(b"\n")
    wgs, nuc = header_accessions(headers)
    return (g, len(seqs), sum(len(s) for s in seqs), h.hexdigest()[:16],
            ",".join(sorted(wgs)), ",".join(sorted(nuc)),
            str(path.resolve()))


def mode_hash(cohort, threads):
    OUT.mkdir(parents=True, exist_ok=True)
    pool = Pool(threads)
    rows = []
    for k, r in enumerate(pool.imap_unordered(one_genome, cohort, chunksize=4)):
        rows.append(r)
        if (k + 1) % 200 == 0:
            print("  hashed {} / {}".format(k + 1, len(cohort)))
    pool.close()
    rows.sort()
    with open(str(OUT / "genome_hashes.tsv"), "w") as fh:
        fh.write("genome\tcontigs\tbp\thash\theader_wgs\theader_acc\tfile\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    by = defaultdict(list)
    for r in rows:
        by[r[3]].append(r)
    dups = [v for v in by.values() if len(v) > 1]
    with open(str(OUT / "duplicate_hashes.tsv"), "w") as fh:
        fh.write("group\tgenome\tcontigs\tbp\thash\tfile\n")
        for i, grp in enumerate(dups, 1):
            for r in grp:
                fh.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(i, r[0], r[1], r[2], r[3], r[6]))
    print("\ngenomes hashed : {}".format(len(rows)))
    print("identical-content groups: {}".format(len(dups)))
    for i, grp in enumerate(dups, 1):
        print("  group {}: {}".format(i, "  ".join("{} ({} contigs, {} bp)".format(
            r[0], r[1], r[2]) for r in grp)))
    print("\nwritten: {}/genome_hashes.tsv, duplicate_hashes.tsv".format(OUT))


def fetch_xml(acc, tries=3):
    url = "https://www.ebi.ac.uk/ena/browser/api/xml/" + acc
    for t in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as e:          # network hiccup: retry, then report
            err = str(e)
            time.sleep(2 * (t + 1))
    return "ERROR: " + err


def mode_ena(cohort):
    OUT.mkdir(parents=True, exist_ok=True)
    gca = [g for g in cohort if g.startswith("GCA_")]
    print("GCA-named genomes in cohort: {}".format(len(gca)))
    rows = []
    for k, g in enumerate(gca, 1):
        x = fetch_xml(g)
        headers = [l[1:].strip() for l in open(str(FASTA / (g + ".fasta"))) if l.startswith(">")]
        fw, fn = header_accessions(headers)
        title = strain = prefix = sample = ""
        chrom = set()
        if x.startswith("ERROR"):
            verdict = "ena_error"
        else:
            try:
                root = ET.fromstring(x)
            except ET.ParseError:
                root = None
            if root is None or root.find(".//ASSEMBLY") is None:
                verdict = "ena_no_record"
            else:
                a = root.find(".//ASSEMBLY")
                title = (a.findtext("TITLE") or "").strip()
                strain = (a.findtext(".//STRAIN") or "").strip()
                prefix = (a.findtext(".//WGS_SET/PREFIX") or "").strip()
                sample = (a.findtext(".//SAMPLE_REF/IDENTIFIERS/PRIMARY_ID") or "").strip()
                for c in a.iter("CHROMOSOME"):
                    ac = c.get("accession", "")
                    if ac:
                        chrom.add(ac.split(".")[0])
                if prefix and fw:
                    verdict = "match" if fw == {prefix} else "MISMATCH"
                elif chrom and fn:
                    verdict = "match" if (fn & chrom) else "MISMATCH"
                else:
                    verdict = "undetermined"
        rows.append((g, verdict, strain, title, prefix, ",".join(sorted(chrom)), sample,
                     ",".join(sorted(fw)), ",".join(sorted(fn)), headers[0][:90] if headers else ""))
        print("  [{:>2}/{}] {:<18} {:<13} ENA strain={:<22} WGS={:<7} file_wgs={:<8} file_acc={}".format(
            k, len(gca), g, verdict, strain[:22], prefix, ",".join(sorted(fw)),
            ",".join(sorted(fn))[:30]))
        time.sleep(0.4)
    with open(str(OUT / "ena_vs_file.tsv"), "w") as fh:
        fh.write("genome\tverdict\tena_strain\tena_title\tena_wgs_prefix\tena_chromosomes\t"
                 "ena_sample\tfile_wgs\tfile_acc\tfile_first_header\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    bad = [r for r in rows if r[1] != "match"]
    print("\nmatch: {}   other: {}".format(len(rows) - len(bad), len(bad)))
    for r in bad:
        print("  {:<18} {:<13} ENA: {} | file: {}".format(r[0], r[1], r[3][:60], r[9][:70]))
    s = defaultdict(list)
    for r in rows:
        if r[6]:
            s[r[6]].append(r[0])
    shared = {k: v for k, v in s.items() if len(v) > 1}
    print("\nGCA genomes sharing one BioSample: {}".format(len(shared)))
    for k, v in shared.items():
        print("  {}: {}".format(k, ", ".join(v)))
    print("\nwritten: {}/ena_vs_file.tsv".format(OUT))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["hash", "ena"])
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()
    cohort = [l.strip() for l in open(str(COHORT)) if l.strip()]
    if args.mode == "hash":
        mode_hash(cohort, args.threads)
    else:
        mode_ena(cohort)


if __name__ == "__main__":
    main()
