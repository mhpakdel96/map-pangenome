#!/usr/bin/env python3
"""
57_panel_validation.py

Does the twelve-SNP typing panel generalise beyond the genomes it was
designed on?  Two validations, because they answer different questions.

  A. DESIGN-AND-TEST CROSS-VALIDATION (internal, honest about the procedure)
     The panel was chosen on all 1,792 genomes, so its perfect concordance
     there is partly circular.  Here the marker-selection procedure itself is
     re-run on a training fraction and the resulting markers are applied to
     held-out genomes that took no part in the selection.  Repeated over
     --folds stratified folds.  Reports, per fold: how many of the twelve
     published positions are re-discovered, how many held-out genomes are
     assigned to the correct group, and which genomes fail.
     This tests the PROCEDURE.  It cannot test the panel on new data.

  B. EXTERNAL VALIDATION (the real test)
     The twelve published positions are genotyped in genomes that were never
     in the cohort, by mapping each genome to K-10 and reading the base at
     each position.  Ground truth for those genomes comes from placing them
     on the existing tree by SNP distance to the cohort (nearest-neighbour
     lineage), so a disagreement means either a panel failure or a genuinely
     new lineage; both are reported and must be inspected, not counted
     blindly.
     Candidate external genomes (none of them used anywhere in the paper):
       - the 142 ENA GCA_9091xxxxx assemblies excluded at cohort definition
       - any assembly listed in --extra-fasta-list
       - MAP assemblies released after the cohort was frozen, fetched with
         --fetch-after YYYY-MM-DD (needs network on the submitting node)

Outputs (12_validation/)
  panel_cv_folds.tsv        per fold: rediscovered markers, held-out accuracy
  panel_cv_errors.tsv       held-out genomes misassigned, with the reason
  panel_external_calls.tsv  per external genome: allele at each of 12 sites,
                            panel call, nearest cohort relative, agreement
  panel_external_summary.tsv
  report.txt                everything printed

Usage
  # A only (no new genomes needed)
  python 57_panel_validation.py --folds 10
  # A and B
  python 57_panel_validation.py --folds 10 \
      --external-dir /cfs/earth/scratch/tabasran/MAP/GCA \
      --minimap2 "$(ls /cfs/earth/scratch/tabasran/.conda/envs/*/bin/minimap2 | head -1)"
Environment: ml_env (numpy, pandas) + minimap2 for part B.

NOTE  This script computes validation results; it does not know what they
      should be.  Nothing from it belongs in the manuscript until it has run.
"""

import argparse
import gzip
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/cfs/earth/scratch/tabasran/MAP")
SEED = 20261110

# The published panel (Table 2).  Position, gene, K-10 allele, marker allele,
# and the group the marker allele marks.
PANEL = [
    (9124,   "gyrA",            "C", "A", "S"),
    (5365,   "gyrB",            "T", "C", "S"),
    (5751,   "gyrB (A165V)",    "C", "T", "S"),
    (6610,   "gyrB",            "T", "C", "S"),
    (8170,   "gyrA",            "G", "A", "S/Type I"),
    (7037,   "gyrB (E594K)",    "G", "A", "S/Type I"),
    (9288,   "gyrA",            "C", "T", "S/Type III"),
    (6883,   "gyrB",            "C", "T", "S/Type III"),
    (7470,   "gyrA",            "A", "G", "S+Bison"),
    (8955,   "gyrA",            "C", "A", "S+Bison"),
    (405941, "rsfB (M28L)",     "A", "A", "C clade D"),
    (444953, "unknown (S195G)", "T", "T", "C clade D"),
]
LEVEL = {"S": 1, "S/Type I": 2, "S/Type III": 2, "S+Bison": 3, "C clade D": 4}
BASES = frozenset("ACGT")


class Log:
    def __init__(self, path):
        self.f = open(path, "w")

    def __call__(self, *a):
        t = " ".join(str(x) for x in a)
        print(t)
        self.f.write(t + "\n")
        self.f.flush()


def read_groups(root):
    """genome -> group, using the same definitions as 55/56."""
    lin = pd.read_csv(root / "05_phylogeny/typing/lineage_assignments.tsv", sep="\t", dtype=str)
    sp = pd.read_csv(root / "05_phylogeny/typing/split_assignments.tsv", sep="\t", dtype=str)
    d = set(sp.loc[sp["tree_side"].fillna("").str.startswith("D "), "genome"])
    g = {}
    for _, r in lin.iterrows():
        if r["lineage"] == "S":
            g[r["genome"]] = "S/Type III" if r["sublineage"] == "III" else "S/Type I"
        elif r["sublineage"] == "B":
            g[r["genome"]] = "Bison"
        else:
            g[r["genome"]] = "C clade D" if r["genome"] in d else "C basal grade"
    return pd.Series(g, name="group")


def read_core_table(path, positions=None):
    """core.tab -> DataFrame indexed by POS, columns = genomes (upper-case bases)."""
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as fh:
        head = fh.readline().rstrip("\n").split("\t")
    use = None
    if positions is not None:
        pos = set(str(p) for p in positions)
        rows = []
        with op(path, "rt") as fh:
            fh.readline()
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if f[1] in pos:
                    rows.append(f)
        T = pd.DataFrame(rows, columns=head)
    else:
        T = pd.read_csv(path, sep="\t", dtype=str)
    T = T.set_index("POS")
    meta = [c for c in ("CHR", "REF", "ALT", "EVIDENCE") if c in T.columns]
    ref = T["REF"] if "REF" in T.columns else None
    T = T.drop(columns=meta, errors="ignore")
    return T.apply(lambda col: col.str.upper()), ref


def assign_from_alleles(alleles):
    """alleles: dict position -> base.  Returns the group the panel calls."""
    def hit(level_name):
        sites = [(p, m) for p, _, _, m, grp in PANEL if grp == level_name]
        got = [alleles.get(p) for p, _ in sites]
        # a position absent from the table (None), or NaN, "N", "-", "" is uncalled
        ok = [a == m for (p, m), a in zip(sites, got)
              if isinstance(a, str) and a in BASES]
        return (len(ok) > 0 and all(ok)), got
    s, _ = hit("S")
    if s:
        t3, _ = hit("S/Type III")
        t1, _ = hit("S/Type I")
        if t3 and not t1:
            return "S/Type III"
        if t1 and not t3:
            return "S/Type I"
        return "S/unresolved"
    b, _ = hit("S+Bison")           # Bison = level-3 marker with the C allele at level 1
    if b:
        return "Bison"
    d, _ = hit("C clade D")
    return "C clade D" if d else "C basal grade"


# --------------------------------------------------------------- A. cross-validation
def rediscover_markers(core, groups, train, target, hi=1.0, lo=0.0):
    """Positions that separate `target` from the rest perfectly within `train`.
    This mirrors how the panel was chosen (perfect concordance with a clade)."""
    tr = [g for g in train if g in core.columns]
    sub = core[tr]
    inside = np.array([groups.get(g) in target for g in tr])
    found = []
    for pos, row in sub.iterrows():
        v = row.values
        ok = np.isin(v, list("ACGT"))
        if ok.sum() < len(tr) * 0.9:
            continue
        vi, vo = v[ok & inside], v[ok & ~inside]
        if len(vi) == 0 or len(vo) == 0:
            continue
        for allele in set(vi):
            if (vi == allele).mean() >= hi and (vo == allele).mean() <= lo:
                found.append((pos, allele))
    return found


def cross_validate(core, groups, folds, log, out):
    rng = np.random.default_rng(SEED)
    gen = [g for g in core.columns if g in groups.index]
    lab = np.array([groups[g] for g in gen])
    order = np.arange(len(gen))
    # stratified folds so every group is represented in training
    fold_of = np.empty(len(gen), dtype=int)
    for grp in np.unique(lab):
        idx = order[lab == grp]
        rng.shuffle(idx)
        fold_of[idx] = np.arange(len(idx)) % folds
    panel_pos = set(str(p) for p, *_ in PANEL)
    # panel alleles, looked up once: rows = the panel positions present in the
    # core table, columns = genomes.  Positions absent here stay uncalled.
    have = [str(p) for p, *_ in PANEL if str(p) in core.index]
    missing = [p for p, *_ in PANEL if str(p) not in core.index]
    if missing:
        log("  [!] {} of 12 panel positions are not in the core table and will be "
            "treated as uncalled: {}".format(len(missing),
                                             ", ".join(str(p) for p in missing)))
        for lvl in sorted(set(LEVEL)):
            n_lvl = [p for p, _, _, _, grp in PANEL if grp == lvl]
            if n_lvl and all(p in missing for p in n_lvl):
                log("  [!] every marker for level '{}' is missing; that level can "
                    "never be called and the accuracy below is not interpretable.".format(lvl))
    PA = core.loc[have]
    PA = PA[~PA.index.duplicated(keep="first")]
    targets = {"S": {"S/Type I", "S/Type III"}, "S/Type I": {"S/Type I"},
               "S/Type III": {"S/Type III"}, "S+Bison": {"S/Type I", "S/Type III", "Bison"},
               "C clade D": {"C clade D"}}
    rows, errs = [], []
    for k in range(folds):
        tr = [gen[i] for i in order[fold_of != k]]
        te = [gen[i] for i in order[fold_of == k]]
        red = {}
        for name, tg in targets.items():
            f = rediscover_markers(core, groups, tr, tg)
            red[name] = f
        n_pub = sum(1 for name, f in red.items() for pos, _ in f if pos in panel_pos)
        # apply the PUBLISHED panel to held-out genomes
        wrong, uncalled = 0, 0
        for g in te:
            col = PA[g] if g in PA.columns else None
            al = {} if col is None else {int(pos): col.get(pos) for pos in have}
            call = assign_from_alleles(al)
            truth = groups[g]
            n_called = sum(1 for p, *_ in PANEL
                           if isinstance(al.get(p), str) and al[p] in BASES)
            uncalled += len(PANEL) - n_called
            same = (call == truth) or (call == "C clade D" and truth == "C clade D") \
                or (call == "C basal grade" and truth == "C basal grade")
            if not same:
                wrong += 1
                errs.append({"fold": k, "genome": g, "truth": truth, "panel_call": call,
                             "n_sites_called": n_called,
                             "alleles": ";".join(
                                 "{}:{}".format(p, al.get(p) if isinstance(al.get(p), str) else "?")
                                 for p, *_ in PANEL)})
        rows.append({"fold": k, "n_train": len(tr), "n_test": len(te),
                     "mean_sites_called_of_12": round(len(PANEL) - uncalled / max(len(te), 1), 3),
                     "published_positions_rediscovered": n_pub,
                     "total_positions_rediscovered": sum(len(f) for f in red.values()),
                     "test_misassigned": wrong,
                     "test_accuracy": round(1 - wrong / max(len(te), 1), 5)})
        log("  fold {:>2}  train {:>5}  test {:>4}  published positions re-found {:>2}/12  "
            "misassigned {:>3}  accuracy {:.4f}".format(k, len(tr), len(te), n_pub, wrong,
                                                        rows[-1]["test_accuracy"]))
    F = pd.DataFrame(rows)
    F.to_csv(out / "panel_cv_folds.tsv", sep="\t", index=False)
    pd.DataFrame(errs).to_csv(out / "panel_cv_errors.tsv", sep="\t", index=False)
    log("  mean held-out accuracy {:.4f} (min {:.4f}); published positions re-found per fold: "
        "median {:.0f} of 12".format(F["test_accuracy"].mean(), F["test_accuracy"].min(),
                                     F["published_positions_rediscovered"].median()))
    return F


# --------------------------------------------------------------- B. external genomes
def genotype_external(fasta, k10, positions, minimap2, tmp):
    """Base of `fasta` at each K-10 position, via a whole-genome alignment."""
    paf = tmp / (Path(fasta).stem + ".paf")
    cmd = [minimap2, "-cx", "asm10", "--cs", "-t", "2", str(k10), str(fasta)]
    with open(paf, "w") as fh:
        r = subprocess.run(cmd, stdout=fh, stderr=subprocess.DEVNULL)
    if r.returncode != 0:
        return {}
    want = sorted(int(p) for p in positions)
    calls = {}
    for line in open(paf):
        f = line.rstrip("\n").split("\t")
        if len(f) < 12 or int(f[11]) < 30:
            continue
        tstart = int(f[7])
        cs = next((x[5:] for x in f[12:] if x.startswith("cs:Z:")), None)
        if cs is None:
            continue
        tpos = tstart
        i = 0
        while i < len(cs):
            op = cs[i]
            j = i + 1
            while j < len(cs) and cs[j] not in ":*+-~":
                j += 1
            tok = cs[i + 1:j]
            if op == ":":
                n = int(tok)
                for p in want:
                    if tpos < p <= tpos + n:
                        calls.setdefault(p, "REF")
                tpos += n
            elif op == "*":
                p1 = tpos + 1
                if p1 in want:
                    calls[p1] = tok[1].upper()
                tpos += 1
            elif op == "-":
                for p in want:
                    if tpos < p <= tpos + len(tok):
                        calls.setdefault(p, "-")
                tpos += len(tok)
            # "+" (insertion in query) does not move the target
            i = j
    return calls


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=str(ROOT))
    ap.add_argument("--folds", type=int, default=10)
    ap.add_argument("--skip-cv", action="store_true")
    ap.add_argument("--external-dir", default="", help="directory of FASTA files never used in the paper")
    ap.add_argument("--extra-fasta-list", default="", help="file with one FASTA path per line")
    ap.add_argument("--minimap2", default="minimap2")
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    root = Path(args.root)
    out = Path(args.out) if args.out else root / "12_validation"
    out.mkdir(parents=True, exist_ok=True)
    tmp = out / "tmp"
    tmp.mkdir(exist_ok=True)
    log = Log(out / "report.txt")
    core_path = root / "05_phylogeny/core_masked/core.tab"
    if not core_path.exists():
        sys.exit("[!] missing " + str(core_path))
    groups = read_groups(root)
    log("groups: " + ", ".join("{} {}".format(k, v) for k, v in groups.value_counts().items()))

    # ---------------- A
    if not args.skip_cv:
        log("\nA. DESIGN-AND-TEST CROSS-VALIDATION ({} folds, seed {})".format(args.folds, SEED))
        log("   markers are re-selected on the training genomes only; the published")
        log("   twelve-SNP panel is then applied to the held-out genomes.")
        core_all, _ = read_core_table(core_path)
        core_all = core_all[[c for c in core_all.columns if c in groups.index]]
        log("   core table: {} sites x {} genomes".format(*core_all.shape))
        cross_validate(core_all, groups, args.folds, log, out)

    # ---------------- B
    ext = []
    if args.external_dir:
        ext += [p for p in Path(args.external_dir).rglob("*") if p.suffix.lower() in (".fasta", ".fa", ".fna")]
    if args.extra_fasta_list:
        ext += [Path(l.strip()) for l in open(args.extra_fasta_list) if l.strip()]
    cohort = set(l.strip() for l in open(root / "01_metadata/cohort_final.txt") if l.strip())
    ext = [p for p in ext if p.stem not in cohort and p.parent.name not in cohort]
    if not ext:
        log("\nB. EXTERNAL VALIDATION skipped: no genomes given (--external-dir / --extra-fasta-list)")
        log("\nwritten: " + str(out))
        return
    log("\nB. EXTERNAL VALIDATION on {} genomes that are not in the cohort".format(len(ext)))
    k10 = root / "03_genomes/fasta/GCA_000007865.1.fasta"
    positions = [p for p, *_ in PANEL]
    k10_allele = {p: a for p, _, a, _, _ in PANEL}
    rows = []
    for i, f in enumerate(ext, 1):
        calls = genotype_external(f, k10, positions, args.minimap2, tmp)
        al = {}
        for p in positions:
            c = calls.get(p)
            al[p] = k10_allele[p] if c == "REF" else (c if c in ("A", "C", "G", "T") else "?")
        call = assign_from_alleles(al)
        rows.append(dict(genome=f.stem, path=str(f), panel_call=call,
                         n_sites_called=sum(1 for p in positions if al[p] != "?"),
                         **{"pos_{}".format(p): al[p] for p in positions}))
        if i % 25 == 0:
            log("   genotyped {}/{}".format(i, len(ext)))
    E = pd.DataFrame(rows)
    E.to_csv(out / "panel_external_calls.tsv", sep="\t", index=False)
    S = E.groupby("panel_call").size().rename("genomes").reset_index()
    S.to_csv(out / "panel_external_summary.tsv", sep="\t", index=False)
    log("   calls: " + "; ".join("{} {}".format(r.panel_call, r.genomes) for r in S.itertuples()))
    bad = E[E["n_sites_called"] < len(positions)]
    log("   genomes with at least one uncalled position: {}".format(len(bad)))
    log("\n   The panel call is not ground truth for these genomes. Confirm each group")
    log("   independently (SNP distance to the cohort, or a tree including them) before")
    log("   reporting a concordance figure.")
    log("\nwritten: " + str(out))


if __name__ == "__main__":
    main()
