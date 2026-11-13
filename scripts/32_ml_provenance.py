#!/usr/bin/env python3
"""
32_ml_provenance.py

Provenance (geographic region) prediction from core-genome SNPs, evaluated so
that the number means what a user of the model would get.

What earlier steps established
  * gene presence/absence carries assembly-quality signal -> NOT used here;
    core SNPs (snippy-core, mapping-based) passed the negative control
  * region has h^2 ~ 1: it is written in the phylogeny; significant SNPs are
    mostly markers of clades sampled by single studies
So the honest questions are (1) does any model beat simply taking the label
of the closest relative on the tree, and (2) how well does it generalise to
a study / a lineage it has never seen.

Models
  tree_1NN   label of the nearest labelled genome by patristic distance on
             the ML tree (from 27's similarity matrix: d_ij = C_ii+C_jj-2C_ij)
  snp_1NN    label of the nearest genome by SNP (Hamming) distance
  logreg_L1  multinomial logistic regression, L1, C chosen by inner 3-fold CV
  rf         random forest, 500 trees, balanced class weights
Evaluation schemes
  random_cv        stratified 5-fold (optimistic: relatives in train and test)
  leave_study_out  each study (>= --min-study genomes) held out in turn
  leave_clade_out  the tree cut into --k-clades clades (average linkage on
                   patristic distance), each held out in turn
Negative controls
  permuted   region labels shuffled (random_cv): must give chance level
  neg_lowN50 within-study low/high N50 half as the label (random_cv): must
             give chance level if SNP features carry no assembly-quality signal
Metric: balanced accuracy (mean per-class recall) on pooled out-of-fold
predictions, plus per-class recall; for leave_study_out also per study.

Outputs (07_ml/)
  metrics.tsv, per_study.tsv, predictions.tsv, top_snps_logreg.tsv, report.txt
"""

import argparse
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, recall_score
from sklearn.exceptions import ConvergenceWarning
from sklearn.model_selection import GridSearchCV, StratifiedKFold

# saga may stop at max_iter on a few inner-CV fits; reported once, not per fit
warnings.filterwarnings("ignore", category=ConvergenceWarning)

ROOT = Path("/cfs/earth/scratch/tabasran/MAP")


def nn_predict(D, tr, te, y):
    """1-NN: for each test index, label of the closest training index."""
    sub = D[np.ix_(te, tr)]
    return y[tr][np.argmin(sub, axis=1)]


def fit_predict(model, X, y, tr, te, D_tree, D_snp, seed, cpus):
    if model == "tree_1NN":
        return nn_predict(D_tree, tr, te, y), None
    if model == "snp_1NN":
        return nn_predict(D_snp, tr, te, y), None
    if model == "logreg_L1":
        base = LogisticRegression(penalty="l1", solver="saga", max_iter=3000,
                                  class_weight="balanced", random_state=seed)
        # inner CV only if every class has >= 3 training samples
        _, cnt = np.unique(y[tr], return_counts=True)
        if len(cnt) > 1 and cnt.min() >= 3:
            m = GridSearchCV(base, {"C": [0.01, 0.1, 1.0]}, cv=3, scoring="balanced_accuracy",
                             n_jobs=cpus)
            m.fit(X[tr], y[tr])
            return m.predict(X[te]), m.best_estimator_
        base.set_params(C=0.1)
        base.fit(X[tr], y[tr])
        return base.predict(X[te]), base
    if model == "rf":
        m = RandomForestClassifier(n_estimators=500, class_weight="balanced", n_jobs=cpus,
                                   random_state=seed)
        m.fit(X[tr], y[tr])
        return m.predict(X[te]), m
    raise ValueError(model)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--snps", default=str(ROOT / "06_pangwas/pyseer/variants_snps.Rtab"))
    ap.add_argument("--kinship", default=str(ROOT / "06_pangwas/pyseer/kinship_all.tsv"))
    ap.add_argument("--traits", default=str(ROOT / "06_pangwas/trait.csv"))
    ap.add_argument("--lineages", default=str(ROOT / "05_phylogeny/typing/lineage_assignments.tsv"))
    ap.add_argument("--ena", default=str(ROOT / "01_metadata/ena_raw.tsv"))
    ap.add_argument("--qc", default=str(ROOT / "02_qc/genome_qc_master.tsv"))
    ap.add_argument("--out", default=str(ROOT / "07_ml"))
    ap.add_argument("--subset", default="C_nonB", choices=["C_nonB", "all"])
    ap.add_argument("--min-study", type=int, default=20)
    ap.add_argument("--k-clades", type=int, default=12)
    ap.add_argument("--neg-min-study", type=int, default=60)
    ap.add_argument("--cpus", type=int, default=8)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--models", default="tree_1NN,snp_1NN,logreg_L1,rf")
    args = ap.parse_args()
    out = Path(args.out) / args.subset
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.RandomState(args.seed)
    models = args.models.split(",")
    R = []

    def say(s=""):
        print(s, flush=True)
        R.append(s)

    # ---------------------------------------------------------- data
    V = pd.read_csv(args.snps, sep="\t", index_col=0)
    tr_ = pd.read_csv(args.traits, dtype=str).set_index("Name")
    lin = pd.read_csv(args.lineages, sep="\t", dtype=str).set_index("genome")
    region = pd.Series(np.nan, index=tr_.index, dtype=object)
    for c in ("Region_Americas", "Region_Europe", "Region_Oceania"):
        region[tr_[c] == "1"] = c.replace("Region_", "")
    keep = region.notna()
    if args.subset == "C_nonB":
        keep &= (lin.reindex(tr_.index)["lineage"] + "/" + lin.reindex(tr_.index)["sublineage"]) == "C/non-B"
    ids = list(tr_.index[keep])
    y = region[ids].values.astype(str)
    X = V[ids].T.values.astype(np.float32)
    var = X.std(axis=0) > 0
    X = X[:, var]
    snp_names = np.array(V.index)[var]
    K = pd.read_csv(args.kinship, sep="\t", index_col=0).loc[ids, ids].values
    dg = np.diag(K)
    D_tree = dg[:, None] + dg[None, :] - 2 * K
    np.fill_diagonal(D_tree, np.inf)
    # Hamming distance on 0/1 allele rows: |x| + |y| - 2 x.y
    nx = X.sum(axis=1).astype(np.float64)
    D_snp = nx[:, None] + nx[None, :] - 2.0 * (X.astype(np.float64) @ X.T.astype(np.float64))
    np.fill_diagonal(D_snp, np.inf)
    ena = pd.read_csv(args.ena, sep="\t", dtype=str).drop_duplicates("run_accession").set_index("run_accession")
    study = np.array([ena.at[g, "study_accession"] if g in ena.index else "assembly" for g in ids])
    qc = pd.read_csv(args.qc, sep="\t").set_index("genome").reindex(ids)
    classes = sorted(set(y))
    say("subset {} | genomes {} | SNP features {} | classes {}".format(
        args.subset, len(ids), X.shape[1], dict(zip(*np.unique(y, return_counts=True)))))

    # ---------------------------------------------------------- schemes
    schemes = {}
    skf = StratifiedKFold(5, shuffle=True, random_state=args.seed)
    schemes["random_cv"] = [(tr, te, "fold{}".format(k)) for k, (tr, te) in enumerate(skf.split(X, y))]
    folds = []
    for s in sorted(set(study)):
        te = np.where(study == s)[0]
        if s == "assembly" or len(te) < args.min_study:
            continue
        folds.append((np.where(study != s)[0], te, s))
    schemes["leave_study_out"] = folds
    Z = linkage(squareform(np.where(np.isinf(D_tree), 0, D_tree), checks=False), "average")
    cl = fcluster(Z, args.k_clades, criterion="maxclust")
    schemes["leave_clade_out"] = [(np.where(cl != c)[0], np.where(cl == c)[0], "clade{}".format(c))
                                  for c in sorted(set(cl))]
    say("leave_study_out folds: {} | leave_clade_out clades: {} (sizes {})".format(
        len(folds), len(set(cl)), sorted(np.bincount(cl)[1:].tolist(), reverse=True)))
    say("clade x region:\n" + pd.crosstab(cl, y).to_string())

    # ---------------------------------------------------------- run
    metrics, per_study, preds, coefs = [], [], [], []

    def run(label_name, yy, scheme_name, folds_, mods):
        for mname in mods:
            t0 = time.time()
            pred = np.array([None] * len(yy), dtype=object)
            unseen = 0
            for tr, te, fname in folds_:
                if len(set(yy[tr])) < 2:
                    continue
                unseen += int(np.sum(~np.isin(yy[te], yy[tr])))
                p, model = fit_predict(mname, X, yy, tr, te, D_tree, D_snp, args.seed, args.cpus)
                pred[te] = p
                if scheme_name == "leave_study_out" and label_name == "region":
                    per_study.append({"model": mname, "study": fname, "n": len(te),
                                      "true": ",".join(sorted(set(yy[te]))),
                                      "accuracy": round(float(np.mean(p == yy[te])), 3)})
                if mname == "logreg_L1" and scheme_name == "random_cv" and label_name == "region" \
                        and model is not None and fname == "fold0":
                    for ci, c in enumerate(model.classes_):
                        w = model.coef_[ci] if model.coef_.shape[0] > 1 else model.coef_[0]
                        for j in np.argsort(-np.abs(w))[:15]:
                            if w[j] != 0:
                                coefs.append({"class": c, "snp": snp_names[j], "coef": round(float(w[j]), 4)})
            ok = np.array([q is not None for q in pred])
            labs = sorted(set(yy[ok]))
            ba = balanced_accuracy_score(yy[ok], pred[ok].astype(str)) if ok.any() else float("nan")
            rec = recall_score(yy[ok], pred[ok].astype(str), labels=labs, average=None,
                               zero_division=0) if ok.any() else []
            row = {"label": label_name, "scheme": scheme_name, "model": mname,
                   "n_predicted": int(ok.sum()), "balanced_accuracy": round(ba, 3),
                   "chance": round(1.0 / max(len(labs), 1), 3),
                   "test_with_class_unseen_in_train": unseen,
                   "seconds": round(time.time() - t0, 1)}
            for c, r in zip(labs, rec):
                row["recall_" + c] = round(float(r), 3)
            metrics.append(row)
            say("  {:<11} {:<16} {:<10} BA {:.3f} (chance {:.3f})  ".format(
                label_name, scheme_name, mname, ba, row["chance"]) +
                "  ".join("{} {:.2f}".format(c, r) for c, r in zip(labs, rec)))
            if label_name == "region":
                for i in np.where(ok)[0]:
                    preds.append({"genome": ids[i], "scheme": scheme_name, "model": mname,
                                  "true": yy[i], "pred": pred[i], "study": study[i], "clade": int(cl[i])})

    say("\n=== region ===")
    for sn, f in schemes.items():
        run("region", y, sn, f, models)

    say("\n=== negative controls (random_cv) ===")
    yperm = y.copy()
    rng.shuffle(yperm)
    run("permuted", yperm, "random_cv", schemes["random_cv"], models)
    neg = np.array([None] * len(ids), dtype=object)
    for s in set(study):
        m = study == s
        if s == "assembly" or m.sum() < args.neg_min_study:
            continue
        med = np.median(qc["n50"].values[m])
        neg[m] = np.where(qc["n50"].values[m] < med, "lowN50", "highN50")
    mk = np.array([v is not None for v in neg])
    if mk.sum() > 100:
        idx_keep = np.where(mk)[0]
        yn = neg[mk].astype(str)
        skf2 = StratifiedKFold(5, shuffle=True, random_state=args.seed)
        f2 = [(idx_keep[a], idx_keep[b], "fold{}".format(k)) for k, (a, b) in enumerate(skf2.split(idx_keep, yn))]
        full = np.array(["NA"] * len(ids), dtype=object)
        full[mk] = yn
        # folds contain only labelled genomes, so 'NA' never enters training or scoring
        run("neg_lowN50", full, "random_cv", f2, models)

    M = pd.DataFrame(metrics)
    M.to_csv(out / "metrics.tsv", sep="\t", index=False)
    pd.DataFrame(per_study).to_csv(out / "per_study.tsv", sep="\t", index=False)
    pd.DataFrame(preds).to_csv(out / "predictions.tsv", sep="\t", index=False)
    pd.DataFrame(coefs).to_csv(out / "top_snps_logreg.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    say("\n=== summary (balanced accuracy) ===")
    say(M.pivot_table(index=["label", "scheme"], columns="model", values="balanced_accuracy").round(3).to_string())
    if per_study:
        PS = pd.DataFrame(per_study)
        say("\n=== leave_study_out, accuracy per held-out study ===")
        say(PS.pivot_table(index=["study", "n", "true"], columns="model", values="accuracy").to_string())
    (out / "report.txt").write_text("\n".join(R) + "\n")
    say("\nwritten: {}".format(out))


if __name__ == "__main__":
    main()
