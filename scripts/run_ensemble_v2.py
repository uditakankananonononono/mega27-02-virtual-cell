"""Ensemble v2: add label-free external evidence (STRING/UniProt/KEGG) to the v1 stack.
Same folds as v1 (StratifiedKFold(3, shuffle, rs=7)); paired bootstrap vs v1."""
import sys, json; sys.path.insert(0, '.')
import numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, average_precision_score
al = pd.read_csv('results/aligned_predictions_all_models.csv')
ext = pd.read_csv('results/external_features.csv')
df = al.merge(ext, on='bnumber', how='left').fillna(0)
v1 = pd.read_csv('results/ensemble_oof.csv')[['bnumber', 'ensemble']]
df = df.merge(v1, on='bnumber')
y = df.essential.values
base = ['fba_min', 'fba_rich', 'gnn', 'cnn', 'kmer']
extc_all = [c for c in ext.columns if c != 'bnumber']
NEW3 = ['pax_log_ppm', 'cai', 'gc3']  # added 9:17 PM (PaxDb + codon usage) -> v3; v2 set frozen as first reported
extc = [c for c in extc_all if not c.startswith('strctx') and c not in NEW3]
# leakage-controlled set: STRING genomic-context+coexpression channels only, no UniProt annotation score
extc_clean = [c for c in extc_all if not c.startswith('string_') and c != 'up_annot' and c not in NEW3]
skf = StratifiedKFold(3, shuffle=True, random_state=7)
def oof(model_fn, cols, seeds=(7,)):
    X = df[cols].values; out = np.zeros(len(y))
    for s in seeds:
        for tr, te in StratifiedKFold(3, shuffle=True, random_state=s).split(X, y):
            out[te] += model_fn().fit(X[tr], y[tr]).predict_proba(X[te])[:, 1] / len(seeds)
    return out
lr = lambda: make_pipeline(StandardScaler(), LogisticRegression(max_iter=4000, class_weight='balanced', C=0.5))
gb = lambda: HistGradientBoostingClassifier(max_iter=300, learning_rate=0.05, max_leaf_nodes=15, l2_regularization=1.0, class_weight='balanced', random_state=0)
res = {}
S = {'v1_stack_reproduced': oof(lr, base), 'ext_only_lr': oof(lr, extc), 'v2_lr': oof(lr, base + extc),
     'v2_gbm': oof(gb, base + extc),
     'v2clean_lr': oof(lr, base + extc_clean), 'v2clean_gbm': oof(gb, base + extc_clean)}
S['v3_lr'] = oof(lr, base + extc + NEW3); S['v3clean_lr'] = oof(lr, base + extc_clean + NEW3)
S['v2clean_avg'] = (pd.Series(S['v2clean_lr']).rank().values + pd.Series(S['v2clean_gbm']).rank().values) / 2
S['v2_avg'] = (pd.Series(S['v2_lr']).rank().values + pd.Series(S['v2_gbm']).rank().values) / 2
for k, s in S.items():
    res[k] = {'auroc': roc_auc_score(y, s), 'auprc': average_precision_score(y, s)}
res['v1_published_oof'] = {'auroc': roc_auc_score(y, df.ensemble), 'auprc': average_precision_score(y, df.ensemble)}
def paired(a, b, n=2000):
    rng = np.random.default_rng(11); d = []
    for _ in range(n):
        i = rng.integers(0, len(y), len(y))
        if y[i].min() == y[i].max(): continue
        d.append(roc_auc_score(y[i], a[i]) - roc_auc_score(y[i], b[i]))
    d = np.array(d); return {'mean': float(d.mean()), 'ci95': [float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))], 'p_le_0': float((d <= 0).mean())}
best = 'v2_lr'  # pre-specified primary model (simplest; avg/gbm are secondary)
res['paired_v2lr_vs_v1published'] = paired(S['v2_lr'], df.ensemble.values)
res['paired_v2lr_vs_v1reproduced'] = paired(S['v2_lr'], S['v1_stack_reproduced'])
res['paired_v3lr_vs_v2lr'] = paired(S['v3_lr'], S['v2_lr'])
res['paired_v3cleanlr_vs_v1reproduced'] = paired(S['v3clean_lr'], S['v1_stack_reproduced'])
res['paired_v2cleanlr_vs_v1reproduced'] = paired(S['v2clean_lr'], S['v1_stack_reproduced'])
# multi-seed robustness of best
seeds = (7, 11, 23, 42, 99)
fn = gb if 'gbm' in best else lr
if best != 'v2_avg':
    res['best_5seed_auroc'] = roc_auc_score(y, oof(fn, base + extc, seeds))
json.dump(res, open('results/ensemble_v2.json', 'w'), indent=2, default=float)
pd.DataFrame({'bnumber': df.bnumber, 'essential': y, **{k: v for k, v in S.items()}}).to_csv('results/ensemble_v2_oof.csv', index=False)
print(json.dumps(res, indent=1, default=float))
