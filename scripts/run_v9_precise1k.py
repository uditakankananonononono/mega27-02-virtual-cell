"""v9: PRECISE-1K expression + iModulon membership (notes/prereg_precise1k.md)."""
import json, re, collections, numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, average_precision_score
al = pd.read_csv('results/aligned_predictions_all_models.csv'); ext = pd.read_csv('results/external_features.csv')
df = al.merge(ext, on='bnumber', how='left').fillna(0); y = df.essential.values
base = ['fba_min', 'fba_rich', 'gnn', 'cnn', 'kmer']; NEW3 = ['pax_log_ppm', 'cai', 'gc3']
extc = [c for c in ext.columns if c != 'bnumber' and not c.startswith('strctx') and c not in NEW3]
X = pd.read_csv('data/external/precise1k/log_tpm_qc.csv', index_col=0)
M = pd.read_csv('data/external/precise1k/M.csv', index_col=0)
th = json.load(open('data/external/precise1k/thresholds.json')); th['Superoxide'] = th['SoxS']  # renamed iModulon (only name mismatch between M.csv and thresholds)
nmod = (M.abs() > pd.Series({c: th[c] for c in M.columns})).sum(axis=1)
feat = pd.DataFrame({'expr_mean': X.mean(axis=1), 'expr_sd': X.std(axis=1), 'n_imod': nmod}).reset_index().rename(columns={'index': 'bnumber'})
df = df.merge(feat, on='bnumber', how='left'); df['expr_missing'] = df.expr_mean.isna().astype(int)
old = df[base + extc].values.astype(float)
NEWC = ['expr_mean', 'expr_sd', 'n_imod', 'expr_missing']
def fold_feats(tr, te, cols=None):
    V = df[NEWC].values.astype(float); med = np.nanmedian(V[tr], axis=0)
    V = np.where(np.isnan(V), med, V)
    return V[tr], V[te]
lr = lambda: make_pipeline(StandardScaler(), LogisticRegression(max_iter=4000, class_weight='balanced', C=0.5))
v2 = np.zeros(len(y)); v9 = np.zeros(len(y)); new_only = np.zeros(len(y))
for tr, te in StratifiedKFold(3, shuffle=True, random_state=7).split(old, y):
    ntr, nte = fold_feats(tr, te)
    v2[te] = lr().fit(old[tr], y[tr]).predict_proba(old[te])[:, 1]
    v9[te] = lr().fit(np.hstack([old[tr], ntr]), y[tr]).predict_proba(np.hstack([old[te], nte]))[:, 1]
    new_only[te] = lr().fit(ntr, y[tr]).predict_proba(nte)[:, 1]
ref = pd.read_csv('results/ensemble_v2_oof.csv').set_index('bnumber').loc[df.bnumber, 'v2_lr'].values
def paired(yy, a, b, f, n=2000):
    rng = np.random.default_rng(11); d = []
    for _ in range(n):
        i = rng.integers(0, len(yy), len(yy))
        if 0 < yy[i].sum() < len(i): d.append(f(yy[i], a[i]) - f(yy[i], b[i]))
    return {'diff': float(f(yy, a) - f(yy, b)), 'ci95': [float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))]}
res = {'n_genes': int(len(y)), 'expr_coverage': int((df.expr_missing == 0).sum()), 'n_samples': int(X.shape[1]),
       'v2_lr_reproduced_max_abs_dev_vs_committed': float(np.abs(v2 - ref).max()),
       'auroc': {'v2_lr': roc_auc_score(y, v2), 'v9_lr': roc_auc_score(y, v9), 'new_features_only': roc_auc_score(y, new_only)},
       'auprc': {'v2_lr': average_precision_score(y, v2), 'v9_lr': average_precision_score(y, v9)},
       'primary_auroc_v9_minus_v2': paired(y, v9, v2, roc_auc_score),
       'secondary_auprc_v9_minus_v2': paired(y, v9, v2, average_precision_score),
       'univariate_auroc': {c: roc_auc_score(y, df[c].fillna(df[c].median())) for c in ['expr_mean', 'expr_sd', 'n_imod']}}
r = pd.read_csv('data/external/rousset2018/pgen.1007749.s012.csv')
up = pd.read_csv('data/external/uniprot_ecoli.tsv', sep='\t'); n2b = {}
for _, x in up.iterrows():
    b = [t for t in str(x['Gene Names (ordered locus)']).split() if re.fullmatch(r'b\d{4}', t)]
    if b and isinstance(x['Gene Names (primary)'], str): n2b.setdefault(x['Gene Names (primary)'], b[0])
r['bnumber'] = r.gene.map(n2b); m = df[['bnumber']].assign(v2=v2, v9=v9).merge(r[['bnumber', 'median_coding']].dropna(), on='bnumber')
yr = (m.median_coding <= -5).astype(int).values
res['secondary_rousset_auroc_v9_minus_v2'] = paired(yr, m.v9.values, m.v2.values, roc_auc_score)
res['verdict'] = 'IMPROVES (pre-registered)' if res['primary_auroc_v9_minus_v2']['ci95'][0] > 0 else 'NO SIGNIFICANT IMPROVEMENT (pre-registered null)'
json.dump(res, open('results/v9_precise1k.json', 'w'), indent=1); print(json.dumps(res, indent=1))
pd.DataFrame({'bnumber': df.bnumber, 'essential': y, 'v9_lr': v9, 'new_features_only': new_only}).to_csv('results/v9_oof.csv', index=False)
