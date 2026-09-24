"""Bernstein et al. 2023 (Mol Syst Biol) RB-TnSeq benchmark as a reusable library.

Model loading, gene/carbon matching, strain/essential/carbon adjustments, simulation and the
precision-recall metric come from the authors' released notebook (MIT licence), executed verbatim.
"""
from __future__ import annotations
import json, os
import numpy as np
from sklearn.metrics import precision_recall_curve, auc, roc_auc_score

NB_CELLS = (9, 11, 13, 15, 17, 19, 21, 23)


def bernstein_prauc(sim: np.ndarray, fit: np.ndarray, thresh: float = 0.001) -> float:
    """Their metric: labels = model growth call (sim>thresh), score = -fitness, positive = no-growth."""
    b = (np.asarray(sim) > thresh).astype(int).ravel()
    f = np.asarray(fit, dtype=float).ravel()
    p, r, _ = precision_recall_curve(b, -f, pos_label=0)
    return float(auc(r, p))


def load_pipeline(base_dir: str) -> dict:
    import cobra, pandas as pd, copy
    nb = json.load(open(os.path.join(base_dir, 'Analysis', 'Analysis_Notebook.ipynb')))
    src = {i: ''.join(c['source']) for i, c in enumerate(nb['cells']) if c['cell_type'] == 'code'}
    g = {'cobra': cobra, 'pd': pd, 'np': np, 'copy': copy}
    for i in NB_CELLS:
        exec(src[i], g)
    return g


def prepare(model_arg: str, base_dir: str):
    import cobra
    base_dir = os.path.abspath(base_dir) + '/'
    g = load_pipeline(base_dir)
    if model_arg.endswith('.xml'):
        model = cobra.io.read_sbml_model(model_arg)
        for ex in model.exchanges:
            ex.lower_bound = 0; ex.upper_bound = 1000
    else:
        model = g['load_model'](model_arg, base_dir)
    med, carb, carb_exp = g['load_environment'](base_dir)
    dexp, dgenes, dfit = g['load_data'](base_dir)
    gm, cem, cmm, fm = g['match_model_data'](model, carb, carb_exp, dexp, dgenes, dfit)
    g['name_genes_matched'] = gm
    m, gma, cema, cmma, fma = g['model_adjustments'](1, 1, 1, model, gm, cem, cmm, fm)
    return g, m, med, list(gma), list(cmma), np.asarray(fma)


def add_supplements(m, ids, lb=-1000.0, prefix='EX_'):
    import cobra
    added = {}
    for mid in ids:
        if mid + '_c' in m.metabolites:
            r = cobra.Reaction(prefix + mid + '_c'); r.lower_bound = lb; r.upper_bound = 1000
            r.add_metabolites({m.metabolites.get_by_id(mid + '_c'): -1.0}); m.add_reactions([r]); added[mid] = r
    return added


def run_benchmark(model_arg: str, base_dir: str, supplements=()) -> dict:
    g, m, med, gma, cmma, fma = prepare(model_arg, base_dir)
    added = list(add_supplements(m, supplements))
    mei, cei = g['check_environment'](m, med, cmma)
    sim = g['simulate_phenotype'](m, gma, cmma, mei, cei)
    b = (sim > 0.001).astype(int)
    return {'model': model_arg, 'supplements_added': added, 'n_genes': len(gma), 'n_carbon': len(cmma),
            'n_pairs': int(b.size), 'n_model_nogrowth': int((b == 0).sum()),
            'pr_auc_bernstein_metric': bernstein_prauc(sim, fma), 'roc_auc': float(roc_auc_score(b.ravel(), fma.ravel())),
            'genes': gma, 'carbon': cmma, 'sim': sim, 'fit': fma}


def greedy_select(sim, fit, pairs, R, cofs, cols):
    """Greedy forward selection of supplements from a single-supplement rescue matrix R (pairs x cofs),
    maximising Bernstein PR-AUC on carbon columns `cols`. Returns (selected, trace)."""
    cols = sorted(cols); b0 = (np.asarray(sim) > 0.001).astype(int)
    def apply(S):
        b = b0.copy()
        if S:
            hit = R[:, [cofs.index(c) for c in S]].max(1) > 0
            for gi, ci in pairs[hit]:
                b[gi, ci] = 1
        return b
    def score(b):
        return _pr(b[:, cols], fit[:, cols])
    S, cur = [], score(b0); trace = [([], cur)]
    while True:
        best = None
        for c in cofs:
            if c in S: continue
            v = score(apply(S + [c]))
            if v > cur + 1e-12 and (best is None or v > best[1]): best = (c, v)
        if best is None: break
        S.append(best[0]); cur = best[1]; trace.append((list(S), cur))
    return S, trace, apply


def _pr(bcalls, fit):
    p, r, _ = precision_recall_curve(np.asarray(bcalls).ravel(), -np.asarray(fit, float).ravel(), pos_label=0)
    return float(auc(r, p))


def paired_gene_bootstrap(b_a, b_b, fit, cols, n=2000, seed=0):
    """Paired bootstrap over genes of PR-AUC(b_a) - PR-AUC(b_b) on carbon columns cols."""
    rng = np.random.default_rng(seed); cols = sorted(cols)
    A, B, F = np.asarray(b_a)[:, cols], np.asarray(b_b)[:, cols], np.asarray(fit)[:, cols]
    d0 = _pr(A, F) - _pr(B, F); ds = []
    for _ in range(n):
        idx = rng.integers(0, A.shape[0], A.shape[0])
        try:
            ds.append(_pr(A[idx], F[idx]) - _pr(B[idx], F[idx]))
        except ValueError:
            continue
    lo, hi = np.percentile(ds, [2.5, 97.5])
    return {'diff': float(d0), 'ci95': [float(lo), float(hi)], 'n_boot': len(ds)}
