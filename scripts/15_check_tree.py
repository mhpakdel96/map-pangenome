#!/usr/bin/env python3
"""
15_check_tree.py

Sanity checks on the phylogenies produced by 14_gubbins_iqtree.sbatch.
Standard library only; compatible with Python 3.5 (no f-strings).

Checks
  1. tips: every cohort genome appears exactly once, nothing extra
  2. support: distribution of UFBoot values on internal nodes (IQ-TREE)
  3. biology: the S-type genomes (> S_MIN SNPs from K-10 in snippy-core,
     i.e. the sheep lineage seen in core.txt) should form ONE clade.
     A correct tree must recover the classical C/S split; if it does not,
     something upstream is wrong.

Usage
  python3 15_check_tree.py
"""

import sys
from pathlib import Path

ROOT = Path("/cfs/earth/scratch/tabasran/MAP")
COHORT = ROOT / "01_metadata/cohort_final.txt"
CORE_TXT = ROOT / "05_phylogeny/core_masked/per_genome_aligned.tsv"
TREES = [("IQ-TREE", ROOT / "05_phylogeny/iqtree/tree_final.nwk"),
         ("Gubbins", ROOT / "05_phylogeny/gubbins/gubbins.final_tree.tre")]
S_MIN = 1000      # SNPs to K-10; C-type genomes sit at tens to low hundreds


def parse(nwk):
    """Return (tip names, list of (leaf-set, label)) for every internal node."""
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


def main():
    cohort = set(l.strip() for l in open(str(COHORT)) if l.strip())
    s_type = set()
    if CORE_TXT.exists():
        for line in open(str(CORE_TXT)):
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3 and f[2].isdigit() and int(f[2]) > S_MIN and f[0] in cohort:
                s_type.add(f[0])
    print("cohort: {}   S-type candidates (> {} SNPs to K-10): {}".format(
        len(cohort), S_MIN, len(s_type)))

    all_ok = True
    for label, path in TREES:
        print("\n" + "=" * 62)
        print(" " + label + "  (" + str(path) + ")")
        print("=" * 62)
        if not path.exists():
            print(" [!] file missing")
            all_ok = False
            continue
        tips, nodes = parse(open(str(path)).read())
        tipset = set(tips)
        dup = len(tips) - len(tipset)
        missing = cohort - tipset
        extra = tipset - cohort
        print(" tips {}   duplicated {}   missing {}   extra {}".format(
            len(tips), dup, len(missing), len(extra)))
        for g in sorted(missing)[:10]:
            print("   missing: " + g)
        for g in sorted(extra)[:10]:
            print("   extra  : " + g)
        if dup or missing or extra:
            all_ok = False

        sup = []
        for _, lab in nodes:
            try:
                sup.append(float(lab.split("/")[-1]))
            except ValueError:
                pass
        if sup:
            sup.sort()
            k = len(sup)
            print(" support values on {} internal nodes: median {:.0f}, "
                  ">=95: {:.1f}%, <70: {:.1f}%".format(
                      k, sup[k // 2],
                      100.0 * sum(1 for x in sup if x >= 95) / k,
                      100.0 * sum(1 for x in sup if x < 70) / k))
        else:
            print(" no support values (expected for the Gubbins tree)")

        if s_type:
            s = s_type & tipset
            rest = tipset - s
            mono = any(ns == s or ns == rest for ns, _ in nodes)
            if mono:
                lab = [l for ns, l in nodes if ns == s or ns == rest]
                print(" S-type ({} genomes): MONOPHYLETIC  (node support: {})".format(
                    len(s), lab[0] if lab and lab[0] else "n/a"))
            else:
                cont = [ns for ns, _ in nodes if s <= ns]
                best = min(cont, key=len) if cont else tipset
                intr = sorted(best - s)
                # A node's leaf-set, or its complement, is every split of the
                # unrooted tree, so "no match" means no branch separates them.
                print(" S-type ({} genomes): NOT monophyletic — no branch separates "
                      "them; smallest clade containing all has {} other genomes".format(
                          len(s), len(intr)))
                for g in intr[:15]:
                    print("   intruder: " + g)

    print("\n" + ("OK" if all_ok else "NOT OK — see [!] above"))
    sys.exit(0 if all_ok else 1)


if __name__ == "__main__":
    main()
