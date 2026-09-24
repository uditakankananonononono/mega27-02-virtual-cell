"""Benchmark metrics: in-silico essentiality vs experimental ground truth."""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.metrics import (accuracy_score, average_precision_score, f1_score,
                             matthews_corrcoef, precision_score, recall_score,
                             roc_auc_score)


def merge_prediction_truth(pred: pd.DataFrame, labels: pd.DataFrame,
                           pred_gene_col: str = "gene") -> pd.DataFrame:
    """Join in-silico predictions to experimental labels on b-number."""
    m = pred.merge(labels, left_on=pred_gene_col, right_on="bnumber", how="inner")
    return m


def classification_report(df: pd.DataFrame, score_col: str, truth_col: str = "essential",
                          pred_col: str | None = None) -> dict:
    """Full metric set. score_col: continuous essentiality score (higher = more
    essential); pred_col: optional hard 0/1 call for point metrics."""
    y = df[truth_col].to_numpy()
    s = df[score_col].to_numpy(dtype=float)
    out = {"n_genes": int(len(y)), "n_essential": int(y.sum()),
           "prevalence": float(y.mean()),
           "auroc": float(roc_auc_score(y, s)),
           "auprc": float(average_precision_score(y, s))}
    if pred_col is not None:
        p = df[pred_col].astype(int).to_numpy()
        out.update({"accuracy": float(accuracy_score(y, p)),
                    "precision": float(precision_score(y, p, zero_division=0)),
                    "recall": float(recall_score(y, p, zero_division=0)),
                    "f1": float(f1_score(y, p, zero_division=0)),
                    "mcc": float(matthews_corrcoef(y, p))})
    return out


def bootstrap_ci(df: pd.DataFrame, score_col: str, n_boot: int = 2000,
                 seed: int = 7) -> dict:
    """Bootstrap 95% CI for AUROC and AUPRC."""
    rng = np.random.default_rng(seed)
    y = df["essential"].to_numpy(); s = df[score_col].to_numpy(dtype=float)
    aucs, aps = [], []
    n = len(y)
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        if len(np.unique(y[idx])) < 2:
            continue
        aucs.append(roc_auc_score(y[idx], s[idx]))
        aps.append(average_precision_score(y[idx], s[idx]))
    return {"auroc_ci": [float(np.percentile(aucs, 2.5)), float(np.percentile(aucs, 97.5))],
            "auprc_ci": [float(np.percentile(aps, 2.5)), float(np.percentile(aps, 97.5))],
            "n_boot_used": len(aucs)}
