#!/usr/bin/env python3
"""
59_external_identity.py

What are the genomes that 57_panel_validation.py part B genotyped?

Part B takes every FASTA in --external-dir that is not in the cohort list.  It
does not check that those files are MAP.  minimap2 will happily align a related
Mycobacterium to K-10 across conserved genes such as gyrA and gyrB, so a
non-MAP genome still receives a panel call, and that call is meaningless.

This script applies the same species criterion the cohort definition used:
Mash distance to K-10.  In the cohort, the most divergent genuine MAP lineage
sat far below the one genome that was excluded as not-MAP (distance 0.0128).
Genomes far from K-10 here are not MAP and their panel calls must be dropped
before any concordance figure is computed.

It reports the distance distribution per panel call, which is the direct test
of whether the surprising Bison and S/Type III calls are real MAP or contaminants.

Outputs (12_validation/)
  external_identity.tsv   per genome: path, panel call, mash distance, verdict
  external_identity_summary.tsv
  identity_report.txt

Usage
  python 59_external_identity.py \
      --mash "$(ls /cfs/earth/scratch/tabasran/.conda/envs/*/bin/mash | head -1)"

NOTE  This script decides nothing about the manuscript.  It tells you what the
      714 genomes are.  Read the verdict column before using any part B number.
"""

import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path("/cfs/earth/scratch/tabasran/MAP")
# the genome excluded from the cohort as not-MAP sat at this distance from K-10
NOT_MAP_EXEMPLAR = 0.0128


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=str(ROOT))
    ap.add_argument("--mash", default="mash")
    ap.add_argument("--calls", default="", help="panel_external_calls.tsv from script 57")
    ap.add_argument("--k10", default="")
    ap.add_argument("--max-map-dist", type=float, default=0.005,
                    help="distance above which a genome is not treated as MAP")
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()

    root = Path(args.root)
    out = root / "12_validation"
    calls_path = Path(args.calls) if args.calls else out / "panel_external_calls.tsv"
    k10 = Path(args.k10) if args.k10 else root / "03_genomes/fasta/GCA_000007865.1.fasta"
    if not calls_path.exists():
        sys.exit("[!] missing " + str(calls_path) + " -- run 57 part B first")
    if not k10.exists():
        sys.exit("[!] missing K-10 reference " + str(k10))

    rep = open(out / "identity_report.txt", "w")

    def log(*a):
        t = " ".join(str(x) for x in a)
        print(t)
        rep.write(t + "\n")
        rep.flush()

    E = pd.read_csv(calls_path, sep="\t", dtype=str)
    log("genomes genotyped by part B: {}".format(len(E)))
    log("panel calls: " + ", ".join(
        "{} {}".format(k, v) for k, v in E["panel_call"].value_counts().items()))

    # where do the files come from, and what are they named?
    E["dir"] = E["path"].map(lambda p: str(Path(p).parent))
    log("\nsource directories:")
    for d, n in E["dir"].value_counts().items():
        log("  {:>5}  {}".format(n, d))
    log("\nfilename prefixes (first 12 characters):")
    pref = E["genome"].str.slice(0, 12).value_counts()
    for p, n in pref.head(15).items():
        log("  {:>5}  {}".format(n, p))
    if len(pref) > 15:
        log("  ... {} further prefixes".format(len(pref) - 15))

    # mash distance to K-10
    sketch = out / "tmp" / "k10"
    sketch.parent.mkdir(parents=True, exist_ok=True)
    log("\nsketching K-10 ...")
    r = subprocess.run([args.mash, "sketch", "-o", str(sketch), str(k10)],
                       capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit("[!] mash sketch failed: " + r.stderr.strip())

    paths = E["path"].tolist()
    log("running mash dist on {} genomes ...".format(len(paths)))
    listfile = out / "tmp" / "external_paths.txt"
    listfile.write_text("\n".join(paths) + "\n")
    r = subprocess.run([args.mash, "dist", "-p", str(args.threads),
                        str(sketch) + ".msh", "-l", str(listfile)],
                       capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit("[!] mash dist failed: " + r.stderr.strip())

    rows = []
    for line in r.stdout.strip().split("\n"):
        f = line.split("\t")
        if len(f) >= 5:
            rows.append({"path": f[1], "mash_dist": float(f[2]),
                         "mash_p": float(f[3]), "shared_hashes": f[4]})
    D = pd.DataFrame(rows)
    M = E.merge(D, on="path", how="left")
    M["verdict"] = ["no mash result" if pd.isna(d)
                    else ("MAP" if d <= args.max_map_dist else "NOT MAP")
                    for d in M["mash_dist"]]

    keep = [c for c in ("genome", "path", "panel_call", "n_sites_called",
                        "mash_dist", "shared_hashes", "verdict") if c in M.columns]
    M[keep].sort_values("mash_dist").to_csv(out / "external_identity.tsv",
                                            sep="\t", index=False)

    log("\nverdict at a Mash distance cutoff of {} (the genome excluded from the "
        "cohort as not-MAP was at {}):".format(args.max_map_dist, NOT_MAP_EXEMPLAR))
    for v, n in M["verdict"].value_counts().items():
        log("  {:<14} {}".format(v, n))

    log("\nMash distance to K-10, by panel call:")
    log("  {:<16} {:>5} {:>9} {:>9} {:>9} {:>7}".format(
        "panel call", "n", "min", "median", "max", "n MAP"))
    g = M.groupby("panel_call")["mash_dist"]
    S = pd.DataFrame({"n": g.size(), "min": g.min(), "median": g.median(),
                      "max": g.max()}).reset_index()
    S["n_MAP"] = M[M["verdict"] == "MAP"].groupby("panel_call").size().reindex(
        S["panel_call"]).fillna(0).astype(int).values
    for r_ in S.itertuples():
        log("  {:<16} {:>5} {:>9.5f} {:>9.5f} {:>9.5f} {:>7}".format(
            r_.panel_call, r_.n, r_.min, r_.median, r_.max, r_.n_MAP))
    S.to_csv(out / "external_identity_summary.tsv", sep="\t", index=False)

    n_map = int((M["verdict"] == "MAP").sum())
    log("\n{} of {} external genomes pass the MAP criterion.".format(n_map, len(M)))
    if n_map < len(M):
        log("The rest are not MAP. Their panel calls are artefacts of mapping a related")
        log("genome to K-10 across conserved genes and must be excluded from part B.")
    log("\nEven among the genomes that pass, the panel call is still not ground truth.")
    log("A concordance figure needs an independent group assignment (SNP distance to")
    log("the cohort, or a tree that includes these genomes).")
    log("\nwritten: " + str(out / "external_identity.tsv"))
    rep.close()


if __name__ == "__main__":
    main()
