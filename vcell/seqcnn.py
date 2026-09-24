"""Sequence layer: 1D-CNN predicting gene essentiality from coding sequence.

Real data: CDS sequences extracted from the U00096.3 reference GenBank
record, labels from Gerdes et al. 2003. Architecture follows the published
DeepCellEss idea (stacked Conv1D on one-hot DNA) at a size this hardware
can actually train, with a k-mer logistic-regression baseline.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from Bio import SeqIO

_ALPH = "ACGT"
_IDX = {b: i for i, b in enumerate(_ALPH)}


def extract_cds(genbank_path: str, min_len: int = 96) -> pd.DataFrame:
    """Extract one CDS row per gene: b-number, gene name, nucleotide sequence."""
    rec = SeqIO.read(genbank_path, "genbank")
    rows, seen = [], set()
    for feat in rec.features:
        if feat.type != "CDS":
            continue
        locus = feat.qualifiers.get("locus_tag", [""])[0]
        if not locus or locus in seen:
            continue
        seq = str(feat.extract(rec.seq)).upper()
        if len(seq) < min_len or any(c not in _IDX for c in seq[:60]):
            continue
        seen.add(locus)
        rows.append({"bnumber": locus,
                     "gene_name": feat.qualifiers.get("gene", [""])[0],
                     "cds_len": len(seq), "cds": seq})
    return pd.DataFrame(rows)


def one_hot(seq: str, max_len: int) -> np.ndarray:
    """One-hot encode a CDS, right-truncated/zero-padded to max_len."""
    x = np.zeros((4, max_len), dtype=np.float32)
    for j, ch in enumerate(seq[:max_len]):
        i = _IDX.get(ch)
        if i is not None:
            x[i, j] = 1.0
    return x


def kmer_freqs(seq: str, k: int = 3) -> np.ndarray:
    """Normalized k-mer frequency vector (4^k dims) - baseline features."""
    v = np.zeros(4 ** k, dtype=np.float32)
    n = 0
    for i in range(0, len(seq) - k + 1):
        idx = 0
        ok = True
        for ch in seq[i:i + k]:
            b = _IDX.get(ch)
            if b is None:
                ok = False; break
            idx = idx * 4 + b
        if ok:
            v[idx] += 1.0; n += 1
    return v / max(n, 1)


def build_arrays(cds: pd.DataFrame, labels: pd.DataFrame, max_len: int = 1200):
    """Join CDS to labels; return one-hot tensor, k-mer matrix, y, bnumbers."""
    m = cds.merge(labels[["bnumber", "essential"]], on="bnumber", how="inner")
    X_seq = np.stack([one_hot(s, max_len) for s in m["cds"]])
    X_kmer = np.stack([kmer_freqs(s) for s in m["cds"]])
    y = m["essential"].to_numpy(dtype=np.int64)
    return X_seq, X_kmer, y, m[["bnumber", "gene_name", "cds_len", "essential"]].reset_index(drop=True)


def make_cnn(max_len: int = 1200):
    """DeepCellEss-style Conv1D stack (compact)."""
    import torch.nn as nn
    return nn.Sequential(
        nn.Conv1d(4, 32, kernel_size=9, padding=4), nn.ReLU(), nn.MaxPool1d(4),
        nn.Conv1d(32, 64, kernel_size=9, padding=4), nn.ReLU(), nn.MaxPool1d(4),
        nn.Conv1d(64, 128, kernel_size=5, padding=2), nn.ReLU(),
        nn.AdaptiveAvgPool1d(1), nn.Flatten(), nn.Linear(128, 1),
    )


def train_eval_cnn(X: np.ndarray, y: np.ndarray, folds: int = 3, epochs: int = 15,
                   batch: int = 64, seed: int = 7, lr: float = 1e-3,
                   max_len: int = 1200) -> dict:
    """Stratified K-fold CV; returns per-fold metrics and out-of-fold scores."""
    import torch
    import torch.nn as nn
    from sklearn.model_selection import StratifiedKFold
    from sklearn.metrics import average_precision_score, roc_auc_score
    torch.manual_seed(seed)
    skf = StratifiedKFold(n_splits=folds, shuffle=True, random_state=seed)
    oof = np.zeros(len(y), dtype=np.float32)
    fold_metrics = []
    for fold, (tr, te) in enumerate(skf.split(X, y)):
        model = make_cnn(max_len)
        n_pos = max(int(y[tr].sum()), 1)
        pos_weight = torch.tensor([(len(tr) - n_pos) / n_pos], dtype=torch.float32)
        opt = torch.optim.Adam(model.parameters(), lr=lr)
        lossf = nn.BCEWithLogitsLoss(pos_weight=pos_weight)
        Xt = torch.from_numpy(X[tr]); yt = torch.from_numpy(y[tr].astype(np.float32))
        n = len(tr)
        for _ in range(epochs):
            perm = torch.randperm(n)
            for i in range(0, n, batch):
                idx = perm[i:i + batch]
                opt.zero_grad()
                out = model(Xt[idx]).squeeze(-1)
                loss = lossf(out, yt[idx])
                loss.backward()
                opt.step()
        with torch.no_grad():
            logits = model(torch.from_numpy(X[te])).squeeze(-1)
        score = torch.sigmoid(logits).numpy()
        oof[te] = score
        fold_metrics.append({
            "fold": fold,
            "auroc": float(roc_auc_score(y[te], score)),
            "auprc": float(average_precision_score(y[te], score))})
    return {"folds": fold_metrics,
            "oof_auroc": float(roc_auc_score(y, oof)),
            "oof_auprc": float(average_precision_score(y, oof)),
            "oof_scores": oof}


def train_eval_kmer_baseline(X_kmer: np.ndarray, y: np.ndarray, folds: int = 3,
                             seed: int = 7) -> dict:
    """Logistic-regression baseline on 3-mer frequencies."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import StratifiedKFold, cross_val_predict
    skf = StratifiedKFold(n_splits=folds, shuffle=True, random_state=seed)
    clf = LogisticRegression(max_iter=2000, class_weight="balanced")
    oof = cross_val_predict(clf, X_kmer, y, cv=skf, method="predict_proba")[:, 1]
    from sklearn.metrics import average_precision_score, roc_auc_score
    return {"oof_auroc": float(roc_auc_score(y, oof)),
            "oof_auprc": float(average_precision_score(y, oof)),
            "oof_scores": oof}
