#!/usr/bin/env python3
"""
16_lineage_typing.py

Assign every genome of the cohort to a MAP lineage (S / C, Type III within S,
Bison-type within C) from the ML tree, using ONLY strains whose type is
stated in a citable publication as anchors. Then test the S sub-typing
independently with the gyrA/gyrB alleles that define Type I vs Type III.

Standard library only; Python 3.5 compatible (no f-strings).

Anchors (literature-verified; see ANCHORS below for the source of each)
  S / Type III   S397, JIII-386
  C              K-10, JII-1961, MAP4, E1, MAP/TANUVAS/TN/India/2008
  C-type (non-B) MAPK_CN7/15, MAPK_CN9/15        (IS1311 PCR-REA, MIRU-VNTR)
  B (Bison) type MAPK_CN4/13, MAPK_JB16/15, MAPK_JJ1/13, S5
  external       Telford (CP033688.1) = S / Type I, not in the cohort; used
                 only as a gyrA/gyrB positive control (--external)

Tree logic (rooting-free). Every internal branch of the unrooted tree splits
the tips into two sides; a node's leaf-set and its complement are those two
sides. A lineage = a side that contains all its anchors and none of the
anchors of the lineage it must be separated from.
  S      largest side with all S anchors and no C/B anchor
  S/III  smallest side with all Type III anchors and no C/B anchor
  C/B    smallest side with all B anchors and no S or C-type(non-B) anchor
Everything else is C (non-B). Labels that cannot be resolved are reported,
never forced.

gyrA/gyrB. Price-Carter et al. 2021 typed sheep strains by these two genes:
Type III differs from all Type I by 2 SNPs in gyrA and 2-3 in gyrB. Alleles
are read from each genome's Snippy alignment (K-10 coordinates) over the
gyrA/gyrB CDS of the K-10 Prokka GFF.

Usage
  python 16_lineage_typing.py [--external Telford=path/to/snps.aligned.fa]

Outputs (05_phylogeny/typing/)
  anchors.tsv               the anchors, their published type and source
  lineage_assignments.tsv   one row per cohort genome
  gyr_haplotypes.tsv        one row per genome (cohort + external)
"""

import argparse
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path("/cfs/earth/scratch/tabasran/MAP")

SRC_JIII = ("Complete genome sequence of ovine MAP strain JIII-386 (MAP-S/type III) "
            "and its comparison to MAP-S/type I, MAP-C and M. avium complex genomes; "
            "Microorganisms 2021; PMC7823733")
SRC_LEE = ("Genomic diversity of MAP: pangenomic approach ... newly constructed complete "
           "genomes; Vet Res 2021; doi:10.1186/s13567-021-00905-1; PMC7977185")
SRC_S5 = ("Singh SV et al. Genome sequence of the 'Indian Bison Type' biotype of MAP "
          "strain S5; Genome Announc 2013;1:e00005; PMC3587920")
SRC_PC = ("Price-Carter M et al. Comparative genomics of MAP sheep strains; "
          "Front Vet Sci 2021; doi:10.3389/fvets.2021.637637")

# genome, strain, published type, group used here, method in source, source
ANCHORS = [
    ("GCA_000219085.3", "S397",           "S / Type III", "S_III", "genome comparison; gyrA/gyrB", SRC_JIII + " | " + SRC_PC),
    ("GCA_002232005.2", "JIII-386",       "S / Type III", "S_III", "genome comparison",            SRC_JIII),
    ("GCA_000007865.1", "K-10",           "C",            "C",     "listed as MAP-C",              SRC_JIII),
    ("GCA_002211525",   "JII-1961",       "C",            "C",     "listed as MAP-C",              SRC_JIII),
    ("GCA_000390085.1", "MAP4",           "C",            "C",     "listed as MAP-C",              SRC_JIII),
    ("GCA_000835225",   "E1",             "C",            "C",     "listed as MAP-C",              SRC_JIII),
    ("GCA_001653355.1", "MAP/TANUVAS/TN/India/2008", "C", "C",     "listed as MAP-C",              SRC_JIII),
    ("GCA_003713025.1", "MAPK_CN7/15",    "C-type",       "C_nonB", "IS1311 PCR-REA; MIRU-VNTR INMV2", SRC_LEE),
    ("GCA_003713045.1", "MAPK_CN9/15",    "C-type",       "C_nonB", "IS1311 PCR-REA; MIRU-VNTR INMV2", SRC_LEE),
    ("GCA_003815795.1", "MAPK_CN4/13",    "B-type",       "B",     "IS1311 PCR-REA; MIRU-VNTR INMV68", SRC_LEE),
    ("GCA_003815815.1", "MAPK_JB16/15",   "B-type",       "B",     "IS1311 PCR-REA; MIRU-VNTR INMV68", SRC_LEE),
    ("GCA_003816035.1", "MAPK_JJ1/13",    "B-type",       "B",     "IS1311 PCR-REA; MIRU-VNTR INMV68", SRC_LEE),
    ("GCA_000330785.1", "S5",             "Indian Bison type", "B", "as stated in title",          SRC_S5),
]
# Published label that the genome contradicts; reported, never used as anchor.
CONFLICTS = [
    ("GCA_021696875.1", "MAP-Gt-9", "S-type",
     "Singh VK et al. 3 Biotech 2024;14:38 (title: 'Complete genome sequence of S-type MAP ...')"),
]


# ------------------------------------------------------------------ tree
def parse(nwk):
    """Return (tip names, list of (leaf-set of names, label)) per internal node."""
    nwk = nwk.strip().rstrip(";")
    stack, nodes, tips = [[]], [], []
    i, n = 0, len(nwk)
    while i < n:
        c = nwk[i]
        if c == "(":
            stack.append([])
            i += 1
        elif c == ",":
            i += 1
        elif c == ")":
            s = set()
            for x in stack.pop():
                s |= x
            i += 1
            j = i
            while j < n and nwk[j] not in ",():;":
                j += 1
            nodes.append((s, nwk[i:j].strip()))
            stack[-1].append(s)
            while j < n and nwk[j] not in ",()":
                j += 1
            i = j
        else:
            j = i
            while j < n and nwk[j] not in ",():":
                j += 1
            name = nwk[i:j].strip().strip("'\"")
            tips.append(name)
            stack[-1].append({name})
            while j < n and nwk[j] not in ",()":
                j += 1
            i = j
    return tips, nodes


def sides_of(tips, nodes):
    """Both sides of every internal branch, as (frozenset, support label)."""
    allt = frozenset(tips)
    out, seen = [], set()
    for s, lab in nodes:
        a = frozenset(s)
        for side in (a, allt - a):
            if 0 < len(side) < len(allt) and side not in seen:
                seen.add(side)
                out.append((side, lab))
    return out


def pick(sides, include, exclude, largest=False):
    cand = [(s, l) for s, l in sides if include <= s and not (s & exclude)]
    if not cand:
        return None, None
    return (max if largest else min)(cand, key=lambda x: len(x[0]))


# ------------------------------------------------------------ gyrA/gyrB
def gyr_regions(gff):
    """{gene: [(start, end)]} for gyrA/gyrB CDS, 1-based inclusive.

    The 1795-run showed that Prokka's gene=gyrA tag in K-10 sits on a 666-bp
    ORF at 4.67 Mb; the real gyrA (7303-9822, 2520 bp, directly downstream of
    gyrB) is tagged gene=gyrA1. Candidates are therefore collected by gene
    name OR by product 'DNA gyrase subunit A/B', all are printed, and the
    LONGEST candidate per subunit is used.
    """
    cand = defaultdict(list)
    for line in open(str(gff)):
        if line.startswith("##FASTA"):
            break
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS":
            continue
        m = re.search(r"(?:^|;)gene=gyr([AB])(?:_?\d+)?(?:;|$)", f[8])
        p = re.search(r"product=DNA gyrase subunit ([AB])(?:;|$)", f[8], re.I)
        sub = (m.group(1) if m else None) or (p.group(1).upper() if p else None)
        if sub:
            prod = re.search(r"product=([^;]*)", f[8])
            cand["gyr" + sub].append((int(f[3]), int(f[4]), f[6],
                                      prod.group(1) if prod else ""))
    reg = {}
    for g in ("gyrA", "gyrB"):
        for a, b, s, pr in cand.get(g, []):
            print("   candidate {}: {}-{} ({} bp, strand {}) '{}'".format(
                g, a, b, b - a + 1, s, pr))
        if cand.get(g):
            a, b, _, _ = max(cand[g], key=lambda x: x[1] - x[0])
            reg[g] = [(a, b)]
    return reg


def read_single_seq(path):
    parts = []
    with open(str(path)) as fh:
        for line in fh:
            if not line.startswith(">"):
                parts.append(line.strip())
    return "".join(parts).upper()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree", default=str(ROOT / "05_phylogeny/iqtree/tree_final.nwk"))
    ap.add_argument("--snps", default=str(ROOT / "05_phylogeny/core_masked/per_genome_aligned.tsv"))
    ap.add_argument("--gff", default=str(ROOT / "03_genomes/gff/GCA_000007865.1.gff"))
    ap.add_argument("--snippy", default=str(ROOT / "05_phylogeny/snippy"))
    ap.add_argument("--cohort", default=str(ROOT / "01_metadata/cohort_final.txt"))
    ap.add_argument("--out", default=str(ROOT / "05_phylogeny/typing"))
    ap.add_argument("--external", action="append", default=[],
                    help="NAME=path/to/snps.aligned.fa (gyr check only)")
    ap.add_argument("--s-min", type=int, default=1000,
                    help="SNPs to K-10 used only for the cross-check (default 1000)")
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    cohort = [l.strip() for l in open(args.cohort) if l.strip()]
    cset = set(cohort)

    snps = {}
    for line in open(args.snps):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 3 and f[2].isdigit():
            snps[f[0]] = int(f[2])

    # ---------------------------------------------------------- anchors
    missing = [a[0] for a in ANCHORS if a[0] not in cset]
    if missing:
        sys.exit("[!] anchors not in cohort: {}".format(missing))
    with open(str(out / "anchors.tsv"), "w") as fh:
        fh.write("genome\tstrain\tpublished_type\tgroup\tmethod\tsource\tsnps_to_K10\n")
        for a in ANCHORS:
            fh.write("\t".join(a) + "\t{}\n".format(snps.get(a[0], "NA")))
    grp = defaultdict(set)
    for a in ANCHORS:
        grp[a[3]].add(a[0])
    S_anc = grp["S_III"]
    C_anc = grp["C"] | grp["C_nonB"] | grp["B"]
    B_anc = grp["B"]
    nonB_anc = grp["C_nonB"]

    # ------------------------------------------------------------- tree
    tips, nodes = parse(open(args.tree).read())
    tset = set(tips)
    if tset != cset or len(tips) != len(tset):
        sys.exit("[!] tree tips do not match the cohort (run 15_check_tree.py)")
    sides = sides_of(tips, nodes)
    print("tree: {} tips, {} distinct splits".format(len(tips), len(sides)))

    S, lS = pick(sides, S_anc, C_anc, largest=True)
    III, lIII = pick(sides, S_anc, C_anc)
    B, lB = pick(sides, B_anc, S_anc | nonB_anc)
    print("\n=== lineages from the tree (anchors only) ===")
    if S is None:
        sys.exit("[!] no branch separates the S anchors from the C anchors")
    print(" S       : {:5d} genomes   support {}".format(len(S), lS or "n/a"))
    print(" S/III   : {:5d} genomes   support {}".format(len(III), lIII or "n/a"))
    rest = frozenset(S - III)
    if rest:
        mono = any(s == rest for s, _ in sides)
        print(" S/non-III: {:5d} genomes   monophyletic: {}".format(
            len(rest), "yes" if mono else "NO (paraphyletic w.r.t. III)"))
    if B is None:
        print(" C/B     : NOT RESOLVED — no branch separates the B anchors from "
              "the S and C-type(non-B) anchors")
        Bx, _ = pick(sides, B_anc, S_anc)
        if Bx is not None:
            print("           smallest side with all B anchors and no S anchor has "
                  "{} genomes; it contains: {}".format(
                      len(Bx), sorted((Bx & nonB_anc) | (Bx & grp["C"]))))
        B = frozenset()
    else:
        print(" C/B     : {:5d} genomes   support {}".format(len(B), lB or "n/a"))
        if "GCA_000007865.1" in B:
            print(" [!] K-10 falls inside the B clade — anchors are inconsistent")
    C = cset - S
    print(" C (all) : {:5d} genomes".format(len(C)))

    # cross-check with the SNP distance to K-10
    hi = set(g for g in cohort if snps.get(g, 0) > args.s_min)
    print("\n cross-check: genomes > {} SNPs to K-10: {} ; identical to tree S: {}".format(
        args.s_min, len(hi), "yes" if hi == set(S) else "NO"))
    for g in sorted(hi ^ set(S))[:20]:
        print("   discordant: {}  snps={}  in_tree_S={}".format(g, snps.get(g), g in S))

    def dist(sel):
        v = sorted(snps[g] for g in sel if g in snps)
        if not v:
            return "n/a"
        return "min {} | median {} | max {}".format(v[0], v[len(v) // 2], v[-1])
    print("\n SNPs to K-10 by lineage")
    print("   S/III      " + dist(III))
    print("   S/non-III  " + dist(rest))
    print("   C/B        " + dist(B))
    print("   C/non-B    " + dist(C - B))

    for g, strain, lab, src in CONFLICTS:
        if g in cset:
            where = "S" if g in S else ("C/B" if g in B else "C/non-B")
            print("\n [conflict] {} ({}) published as '{}' but tree places it in {} "
                  "with {} SNPs to K-10\n            source: {}".format(
                      g, strain, lab, where, snps.get(g, "NA"), src))

    # -------------------------------------------------------- gyrA/gyrB
    print("\n=== gyrA / gyrB (K-10 Prokka GFF) ===")
    reg = gyr_regions(args.gff)
    for gname in ("gyrA", "gyrB"):
        print(" {}: {}".format(gname, reg.get(gname, "NOT FOUND")))
    gyr = {}
    if not reg.get("gyrA") or not reg.get("gyrB"):
        print(" [!] gyrA/gyrB not found in the GFF — gyr test skipped")
    else:
        pos = []
        for gname in ("gyrA", "gyrB"):
            for a, b in reg[gname]:
                pos.extend((gname, p) for p in range(a, b + 1))
        files = [(g, Path(args.snippy) / g / "snps.aligned.fa") for g in cohort]
        for e in args.external:
            name, _, path = e.partition("=")
            files.append((name, Path(path)))
        alle = {}
        for k, (g, f) in enumerate(files):
            if not f.exists():
                sys.exit("[!] missing {}".format(f))
            s = read_single_seq(f)
            alle[g] = "".join(s[p - 1] if s[p - 1] in "ACGT" else "N" for _, p in pos)
            if (k + 1) % 200 == 0:
                print("   read {} / {}".format(k + 1, len(files)))
        ref = alle["GCA_000007865.1"]
        var = [i for i in range(len(pos))
               if len(set(a[i] for a in alle.values()) - {"N"}) > 1]
        print(" variable positions in gyrA+gyrB: {}".format(len(var)))
        for i in var:
            c = Counter(a[i] for a in alle.values())
            print("   {} {:>8}  K-10={}  {}".format(pos[i][0], pos[i][1], ref[i], dict(c)))
        for g, a in alle.items():
            gyr[g] = "".join(a[i] for i in var)
        names = {a[0]: a[1] for a in ANCHORS}
        hap = defaultdict(list)
        for g in alle:
            hap[gyr[g]].append(g)
        print("\n haplotype (variable positions only) -> genomes, with tree lineage")
        for h, gs in sorted(hap.items(), key=lambda x: -len(x[1])):
            lin = Counter("external" if g not in cset else
                          ("S/III" if g in III else "S/non-III" if g in S else
                           "C/B" if g in B else "C/non-B") for g in gs)
            anc = [names.get(g, g) for g in gs if g in names or g not in cset]
            print("   {:<12} n={:<5} {}  anchors/external: {}".format(
                h or "(=K-10)", len(gs), dict(lin), ", ".join(anc) or "-"))
        with open(str(out / "gyr_haplotypes.tsv"), "w") as fh:
            fh.write("genome\tin_cohort\tgyr_haplotype\n")
            for g in sorted(alle):
                fh.write("{}\t{}\t{}\n".format(g, g in cset, gyr[g]))

    # ------------------------------------------------------------ write
    with open(str(out / "lineage_assignments.tsv"), "w") as fh:
        fh.write("genome\tsnps_to_K10\tlineage\tsublineage\tgyr_haplotype\tanchor\n")
        names = {a[0]: a[1] + " (" + a[2] + ")" for a in ANCHORS}
        for g in sorted(cohort):
            if g in S:
                lin, sub = "S", ("III" if g in III else "non-III")
            else:
                lin, sub = "C", ("B" if g in B else "non-B")
            fh.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(
                g, snps.get(g, "NA"), lin, sub, gyr.get(g, "NA"), names.get(g, "")))
    print("\n written: {}/anchors.tsv, lineage_assignments.tsv{}".format(
        out, ", gyr_haplotypes.tsv" if gyr else ""))


if __name__ == "__main__":
    main()
