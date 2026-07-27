#!/usr/bin/env python3
"""Cross-validated TabPFN AUC for the metabolomics classifier benchmark.

Reads a feature matrix and per-sample (label, fold) assignments, runs TabPFN on
each CV split (train on out-of-fold samples, predict the held-out fold), pools
the held-out positive-class probabilities, and writes AUC + ROC points as JSON.

Run as an ISOLATED subprocess by R/domain/metabolomics/04c_classifier_benchmark.R
— the main R pipeline never imports TabPFN. The model runs LOCALLY only; nothing
is sent off-machine.

Errors are carried in the JSON payload (never raised to a non-zero exit) so the R
caller can report them and fall back to the RF-only benchmark gracefully.

Usage:
    python tabpfn_cv.py --x X.tsv --labels labels.tsv --out result.json --seed 1234

Inputs (tab-separated, with header):
    --x       first column `sample_id`, remaining columns are features (rows=samples)
    --labels  columns: sample_id, label (0/1), fold (int)
"""
import argparse
import json
import sys

import numpy as np
import pandas as pd


def _make_classifier(TabPFNClassifier, seed):
    """Construct a TabPFNClassifier tolerant of API differences across versions."""
    for kwargs in (
        {"device": "cpu", "random_state": seed, "ignore_pretraining_limits": True},
        {"device": "cpu", "random_state": seed},
        {"device": "cpu"},
        {},
    ):
        try:
            return TabPFNClassifier(**kwargs)
        except TypeError:
            continue
    return TabPFNClassifier()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--x", required=True)
    ap.add_argument("--labels", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=1234)
    args = ap.parse_args()

    result = {
        "auc": None, "fpr": [], "tpr": [],
        "n_samples": 0, "n_features": 0, "error": None,
    }

    try:
        np.random.seed(args.seed)

        X = pd.read_csv(args.x, sep="\t", index_col=0)
        lab = pd.read_csv(args.labels, sep="\t", index_col=0).loc[X.index]
        y = lab["label"].to_numpy().astype(int)
        folds = lab["fold"].to_numpy().astype(int)
        result["n_samples"] = int(X.shape[0])
        result["n_features"] = int(X.shape[1])

        try:
            import torch
            torch.manual_seed(args.seed)
        except Exception:
            pass

        from tabpfn import TabPFNClassifier
        from sklearn.metrics import roc_auc_score, roc_curve

        Xv = X.to_numpy()
        oof = np.full(len(y), np.nan, dtype=float)
        for f in np.unique(folds):
            tr = folds != f
            te = folds == f
            if len(np.unique(y[tr])) < 2:
                continue  # cannot train a 2-class model on a single-class split
            clf = _make_classifier(TabPFNClassifier, args.seed)
            clf.fit(Xv[tr], y[tr])
            proba = clf.predict_proba(Xv[te])
            classes = list(clf.classes_)
            pos = classes.index(1) if 1 in classes else proba.shape[1] - 1
            oof[te] = proba[:, pos]

        mask = ~np.isnan(oof)
        if mask.sum() < 2 or len(np.unique(y[mask])) < 2:
            result["error"] = "insufficient held-out predictions for ROC"
        else:
            result["auc"] = float(roc_auc_score(y[mask], oof[mask]))
            fpr, tpr, _ = roc_curve(y[mask], oof[mask])
            result["fpr"] = [float(v) for v in fpr]
            result["tpr"] = [float(v) for v in tpr]
    except Exception as e:  # noqa: BLE001 — surface any failure to the R caller
        result["error"] = f"{type(e).__name__}: {e}"

    with open(args.out, "w") as fh:
        json.dump(result, fh)
    return 0


if __name__ == "__main__":
    sys.exit(main())
