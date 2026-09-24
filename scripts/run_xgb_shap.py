"""Exploratory (not pre-registered): XGBoost on the v2 feature set, same folds, paired vs the pre-specified v2 LR;
TreeSHAP attribution (xgboost pred_contribs) on the full-data model to show which evidence drives essentiality calls."""
import json, sys; sys.path.insert(0, '.')
import numpy as np, pandas as pd, xgboost as xgb
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, average_precision_score
al = pd.read_csv('results/aligned_predictions_all_models.csv'); ext = pd.read_csv('results/external_features.csv')
df = al.merge(ext, on='bnumber', how='left').fillna(0); y = df.essential.values
base = ['fba_min', 'fba_rich', 'gnn', 'cnn', 'kmer']
NEW3 = ['pax_log_ppm', 'cai', 'gc3']
extc = [c for c in ext.columns if c != 'bnumber' and not c.startswith('strctx') and c not in NEW3]
cols = base + extc; X = df[cols].values
pos = y.mean()
mk = lambda: xgb.XGBClassifier(n_estimators=300, max_depth=3, learning_rate=0.05, subsample=0.8, colsample_bytree=0.8,
                               scale_pos_weight=(1 - pos) / pos, n_jobs=1, random_state=0, eval_metric='logloss')
oof = np.zeros(len(y))
for tr, te in StratifiedKFold(3, shuffle=True, random_state=7).split(X, y):
    oof[te] = mk().fit(X[tr], y[tr]).predict_proba(X[te])[:, 1]
v2 = pd.read_csv('results/ensemble_v2_oof.csv').set_index('bnumber').loc[df.bnumber, 'v2_lr'].values
rng = np.random.default_rng(11); d = []
for _ in range(2000):
    i = rng.integers(0, len(y), len(y)); d.append(roc_auc_score(y[i], oof[i]) - roc_auc_score(y[i], v2[i]))
d = np.array(d)
full = mk().fit(X, y)
# TreeSHAP via xgboost's own pred_contribs (shap 0.49 cannot parse xgboost 3 base_score); last column is the bias term
sv = full.get_booster().predict(xgb.DMatrix(X), pred_contribs=True)[:, :-1]
imp = pd.Series(np.abs(sv).mean(0), index=cols).sort_values(ascending=False)
res = {'xgb_auroc': roc_auc_score(y, oof), 'xgb_auprc': average_precision_score(y, oof), 'v2_lr_auroc': roc_auc_score(y, v2),
       'paired_xgb_minus_v2lr': {'mean': float(d.mean()), 'ci95': [float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))]},
       'shap_mean_abs_top10': {k: float(v) for k, v in imp.head(10).items()}, 'n_features': len(cols),
       'note': 'exploratory; primary model remains v2 LR (pre-specified)'}
json.dump(res, open('results/xgb_shap.json', 'w'), indent=1); print(json.dumps(res, indent=1))
