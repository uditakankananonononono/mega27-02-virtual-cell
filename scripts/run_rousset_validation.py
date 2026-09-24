"""Pre-registered external validation on Rousset 2018 CRISPRi (notes/prereg_rousset_validation.md)."""
import json
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score
r = pd.read_csv('data/external/rousset2018/pgen.1007749.s012.csv').rename(columns={'essential': 'rousset_call'})
up = pd.read_csv('data/external/uniprot_ecoli.tsv', sep='\t')
n2b = {}
for _, x in up.iterrows():
    b = [t for t in str(x['Gene Names (ordered locus)']).split() if t.startswith('b') and t[1:].isdigit()]
    if b and isinstance(x['Gene Names (primary)'], str): n2b.setdefault(x['Gene Names (primary)'], b[0])
r['bnumber'] = r.gene.map(n2b)
v2 = pd.read_csv('results/ensemble_v2_oof.csv')
al = pd.read_csv('results/aligned_predictions_all_models.csv')[['bnumber', 'fba_min']]
d = r.dropna(subset=['bnumber']).merge(v2[['bnumber', 'essential', 'v2_lr']], on='bnumber').merge(al, on='bnumber')
y = (d.median_coding <= -5).astype(int).values
a, b = d.v2_lr.values, d.fba_min.values
rng = np.random.default_rng(11); diff = []
for _ in range(2000):
    i = rng.integers(0, len(y), len(y)); diff.append(roc_auc_score(y[i], a[i]) - roc_auc_score(y[i], b[i]))
diff = np.array(diff); ne = d.essential.values == 0
res = {'n_rousset_genes': int(len(r)), 'n_mapped': int(r.bnumber.notna().sum()), 'n_evaluated': int(len(d)), 'n_crispri_essential': int(y.sum()),
       'auroc_v2_lr': float(roc_auc_score(y, a)), 'auroc_fba_min': float(roc_auc_score(y, b)),
       'paired_v2_minus_fba': {'mean': float(diff.mean()), 'ci95': [float(np.percentile(diff, 2.5)), float(np.percentile(diff, 97.5))]},
       'spearman_v2_vs_depletion': float(spearmanr(a, -d.median_coding, nan_policy='omit').correlation),
       'gerdes_nonessential_subset': {'n': int(ne.sum()), 'n_crispri_essential': int(y[ne].sum()),
                                      'auroc_v2_lr': float(roc_auc_score(y[ne], a[ne])) if 0 < y[ne].sum() < ne.sum() else None,
                                      'spearman': float(spearmanr(a[ne], -d.median_coding[ne], nan_policy='omit').correlation)}}
res['exploratory_agreement_with_authors_call'] = float(((d.rousset_call.astype(str).str.upper()=='TRUE').astype(int).values == y).mean())
res['verdict'] = 'GENERALISES' if (res['auroc_v2_lr'] > 0.70 and res['paired_v2_minus_fba']['ci95'][0] > 0) else 'DOES NOT MEET PRE-REGISTERED BAR'
json.dump(res, open('results/rousset_validation.json', 'w'), indent=1); print(json.dumps(res, indent=1))
