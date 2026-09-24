"""Graph layer: GNN over the gene-reaction-metabolite network.

Graph construction from the stoichiometric model itself (no external
labels): two genes share an edge when their associated reactions share a
metabolite (currency metabolites excluded). A 2-layer GCN with symmetric
normalization, implemented in pure PyTorch so the hermetic test suite has
no compiled dependencies beyond torch.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

CURRENCY = {"h", "h2o", "atp", "adp", "amp", "nad", "nadh", "nadp", "nadph",
            "coa", "pi", "ppi", "co2", "o2", "nh4", "h2o2"}


def build_gene_graph(model, genes: list[str]) -> np.ndarray:
    """Adjacency (NxN, 0/1) over `genes`: edge if reactions share a metabolite."""
    gene_set = set(genes)
    gidx = {g: i for i, g in enumerate(genes)}
    A = np.zeros((len(genes), len(genes)), dtype=np.float32)
    # reaction -> genes (via GPR), reaction -> core metabolites
    rxn_genes = []
    rxn_mets = []
    for rxn in model.reactions:
        gs = {g.id for g in rxn.genes} & gene_set
        ms = {m.id.rsplit("_", 1)[0] for m in rxn.metabolites} - CURRENCY
        rxn_genes.append(gs)
        rxn_mets.append(ms)
    # met -> reactions, then gene pairs via shared metabolite
    from collections import defaultdict
    met_rxns = defaultdict(list)
    for i, ms in enumerate(rxn_mets):
        for m in ms:
            met_rxns[m].append(i)
    gene_pairs = set()
    for rxns in met_rxns.values():
        if len(rxns) > 25:   # skip hyper-connected hubs (e.g. proton pumps)
            continue
        for i in rxns:
            for j in rxns:
                if i < j:
                    for gi in rxn_genes[i]:
                        for gj in rxn_genes[j]:
                            if gi != gj:
                                gene_pairs.add((gi, gj))
    for gi, gj in gene_pairs:
        A[gidx[gi], gidx[gj]] = 1.0
        A[gidx[gj], gidx[gi]] = 1.0
    return A


def gcn_norm(A: np.ndarray) -> np.ndarray:
    """Symmetric-normalized adjacency with self-loops: D^-1/2 (A+I) D^-1/2."""
    n = A.shape[0]
    Ah = A + np.eye(n, dtype=np.float32)
    d = Ah.sum(1)
    dinv = 1.0 / np.sqrt(np.maximum(d, 1e-8))
    return (Ah * dinv[:, None]) * dinv[None, :]


def node_features(model, genes: list[str], cds_index: pd.DataFrame | None = None) -> np.ndarray:
    """Structural node features: degree in GPR network, #reactions, #metabolites,
    boundary participation, CDS length (z-scored) when available."""
    n = len(genes)
    gidx = {g: i for i, g in enumerate(genes)}
    nrxn = np.zeros(n); nmet = np.zeros(n); boundary = np.zeros(n)
    for rxn in model.reactions:
        gs = {g.id for g in rxn.genes} & set(genes)
        ms = {m.id.rsplit("_", 1)[0] for m in rxn.metabolites} - CURRENCY
        is_exch = rxn in model.exchanges or rxn.boundary
        for g in gs:
            nrxn[gidx[g]] += 1
            nmet[gidx[g]] += len(ms)
            boundary[gidx[g]] += 1.0 if is_exch else 0.0
    feats = [nrxn, nmet, boundary]
    if cds_index is not None:
        ln = cds_index.set_index("bnumber")["cds_len"].reindex(genes).fillna(0).to_numpy()
        feats.append(np.log1p(ln))
    X = np.stack(feats, 1).astype(np.float32)
    X = (X - X.mean(0)) / (X.std(0) + 1e-8)
    return X


def train_eval_gnn(A: np.ndarray, X: np.ndarray, y: np.ndarray, folds: int = 3,
                   epochs: int = 200, hidden: int = 32, seed: int = 7,
                   lr: float = 5e-3, weight_decay: float = 1e-4) -> dict:
    """Stratified CV for a 2-layer GCN; returns OOF metrics and scores."""
    import torch
    import torch.nn as nn
    from sklearn.model_selection import StratifiedKFold
    from sklearn.metrics import average_precision_score, roc_auc_score
    torch.manual_seed(seed)
    An = torch.from_numpy(gcn_norm(A))
    Xt_all = torch.from_numpy(X)
    skf = StratifiedKFold(n_splits=folds, shuffle=True, random_state=seed)
    oof = np.zeros(len(y), dtype=np.float32)
    fold_metrics = []
    d = X.shape[1]
    for fold, (tr, te) in enumerate(skf.split(X, y)):
        W1 = torch.nn.Parameter(torch.randn(d, hidden) * 0.1)
        W2 = torch.nn.Parameter(torch.randn(hidden, 1) * 0.1)
        params = [W1, W2]
        opt = torch.optim.Adam(params, lr=lr, weight_decay=weight_decay)
        n_pos = max(int(y[tr].sum()), 1)
        pos_weight = torch.tensor([(len(tr) - n_pos) / n_pos], dtype=torch.float32)
        lossf = nn.BCEWithLogitsLoss(pos_weight=pos_weight)
        yt = torch.from_numpy(y.astype(np.float32))
        tr_idx = torch.from_numpy(tr)
        for _ in range(epochs):
            opt.zero_grad()
            H = torch.relu(An @ Xt_all @ W1)
            logits = (An @ H @ W2).squeeze(-1)
            loss = lossf(logits[tr_idx], yt[tr_idx])
            loss.backward()
            opt.step()
        with torch.no_grad():
            H = torch.relu(An @ Xt_all @ W1)
            logits = (An @ H @ W2).squeeze(-1)
        score = torch.sigmoid(logits[te]).numpy()
        oof[te] = score
        fold_metrics.append({"fold": fold,
                             "auroc": float(roc_auc_score(y[te], score)),
                             "auprc": float(average_precision_score(y[te], score))})
    return {"folds": fold_metrics,
            "oof_auroc": float(roc_auc_score(y, oof)),
            "oof_auprc": float(average_precision_score(y, oof)),
            "oof_scores": oof}
