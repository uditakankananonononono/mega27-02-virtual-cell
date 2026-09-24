"""v7: OMA HOG depth (amendment 1) (notes/prereg_oma_conservation.md)."""
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
om = pd.read_csv('data/external/oma/oma_hog_levels.tsv', sep='\t').drop_duplicates('bnumber')
df = df.merge(om[['bnumber', 'n_hog_levels']], on='bnumber', how='left')
df['log_orth'] = np.log1p(pd.to_numeric(df.n_hog_levels, errors='coerce').fillna(0))
old = df[base + extc].values.astype(float)
def fold_feats(tr, te):
    return df[['log_orth']].values[tr], df[['log_orth']].values[te]
lr = lambda: make_pipeline(StandardScaler(), LogisticRegression(max_iter=4000, class_weight='balanced', C=0.5))
v2 = np.zeros(len(y)); v7 = np.zeros(len(y)); new_only = np.zeros(len(y))
for tr, te in StratifiedKFold(3, shuffle=True, random_state=7).split(old, y):
    ntr, nte = fold_feats(tr, te)
    v2[te] = lr().fit(old[tr], y[tr]).predict_proba(old[te])[:, 1]
    v7[te] = lr().fit(np.hstack([old[tr], ntr]), y[tr]).predict_proba(np.hstack([old[te], nte]))[:, 1]
    new_only[te] = lr().fit(ntr, y[tr]).predict_proba(nte)[:, 1]
ref = pd.read_csv('results/ensemble_v2_oof.csv').set_index('bnumber').loc[df.bnumber, 'v2_lr'].values
def paired(yy, a, b, f, n=2000):
    rng = np.random.default_rng(11); d = []
    for _ in range(n):
        i = rng.integers(0, len(yy), len(yy))
        if 0 < yy[i].sum() < len(i): d.append(f(yy[i], a[i]) - f(yy[i], b[i]))
    return {'diff': float(f(yy, a) - f(yy, b)), 'ci95': [float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))]}
res = {'n_genes': int(len(y)), 'oma_found': int(df.n_hog_levels.notna().sum()), 'oma_nonzero': int((df.log_orth > 0).sum()),
       'v2_lr_reproduced_max_abs_dev_vs_committed': float(np.abs(v2 - ref).max()),
       'auroc': {'v2_lr': roc_auc_score(y, v2), 'v7_lr': roc_auc_score(y, v7), 'new_features_only': roc_auc_score(y, new_only)},
       'auprc': {'v2_lr': average_precision_score(y, v2), 'v7_lr': average_precision_score(y, v7)},
       'primary_auroc_v7_minus_v2': paired(y, v7, v2, roc_auc_score),
       'secondary_auprc_v7_minus_v2': paired(y, v7, v2, average_precision_score),
       'univariate_auroc': {c: roc_auc_score(y, df[c].fillna(df[c].median())) for c in ['log_orth']}}
r = pd.read_csv('data/external/rousset2018/pgen.1007749.s012.csv')
up = pd.read_csv('data/external/uniprot_ecoli.tsv', sep='\t'); n2b = {}
for _, x in up.iterrows():
    b = [t for t in str(x['Gene Names (ordered locus)']).split() if re.fullmatch(r'b\d{4}', t)]
    if b and isinstance(x['Gene Names (primary)'], str): n2b.setdefault(x['Gene Names (primary)'], b[0])
r['bnumber'] = r.gene.map(n2b); m = df[['bnumber']].assign(v2=v2, v7=v7).merge(r[['bnumber', 'median_coding']].dropna(), on='bnumber')
yr = (m.median_coding <= -5).astype(int).values
res['secondary_rousset_auroc_v7_minus_v2'] = paired(yr, m.v7.values, m.v2.values, roc_auc_score)
res['verdict'] = 'IMPROVES (pre-registered)' if res['primary_auroc_v7_minus_v2']['ci95'][0] > 0 else 'NO SIGNIFICANT IMPROVEMENT (pre-registered null)'
json.dump(res, open('results/v7_oma.json', 'w'), indent=1); print(json.dumps(res, indent=1))
pd.DataFrame({'bnumber': df.bnumber, 'essential': y, 'v7_lr': v7, 'new_features_only': new_only}).to_csv('results/v7_oof.csv', index=False)
