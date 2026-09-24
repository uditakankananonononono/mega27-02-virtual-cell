"""Per-dataset PaxDb abundance robustness (notes/prereg_paxdb_datasets.md)."""
import glob, json, numpy as np, pandas as pd
from scipy.stats import binomtest, spearmanr, norm
from sklearn.metrics import roc_auc_score
al = pd.read_csv('results/aligned_predictions_all_models.csv')[['bnumber', 'essential']]
per = {}
for f in sorted(glob.glob('data/external/paxdb511145/511145-*.txt')):
    meta = {l[1:].split(':', 1)[0]: l.split(':', 1)[1].strip() for l in open(f) if l.startswith('#') and ':' in l}
    d = pd.read_csv(f, sep='\t', comment='#', header=None, usecols=[0, 1, 2], names=['gene_name', 'sid', 'abundance'])  # some files carry a 4th raw-count column
    d['bnumber'] = d.sid.astype(str).str.replace('511145.', '', regex=False)
    d = d[d.abundance > 0].groupby('bnumber').abundance.max().reset_index()
    m = al.merge(d, on='bnumber')
    y, x = m.essential.values, np.log10(m.abundance.values)
    rec = {'paxdb_id': meta.get('id'), 'name': meta.get('name'), 'n_covered': int(len(m)), 'n_essential_covered': int(y.sum())}
    if 0 < y.sum() < len(y):
        a = roc_auc_score(y, x); n1, n0 = y.sum(), len(y) - y.sum()
        q1, q2 = a / (2 - a), 2 * a * a / (1 + a)
        se = np.sqrt((a * (1 - a) + (n1 - 1) * (q1 - a * a) + (n0 - 1) * (q2 - a * a)) / (n1 * n0))
        rec.update(auroc=float(a), se=float(se))
    per[f.split('/')[-1]] = rec
sc = [v for v in per.values() if 'auroc' in v]
k = sum(v['auroc'] > 0.5 for v in sc); n = sum(v['auroc'] != 0.5 for v in sc)
w = np.array([1 / v['se'] ** 2 for v in sc]); a = np.array([v['auroc'] for v in sc])
pooled = float((w * a).sum() / w.sum()); pse = float(np.sqrt(1 / w.sum()))
res = {'n_datasets': len(per), 'n_scored': len(sc), 'n_auroc_gt_0.5': int(k),
       'sign_test_p': float(binomtest(k, n, 0.5).pvalue), 'pooled_auroc_ivw': pooled, 'pooled_ci95': [pooled - 1.96 * pse, pooled + 1.96 * pse],
       'spearman_auroc_vs_coverage': float(spearmanr([v['auroc'] for v in sc], [v['n_covered'] for v in sc]).correlation),
       'auroc_range': [float(a.min()), float(a.max())], 'per_dataset': per}
res['verdict'] = 'ROBUST (pre-registered)' if (res['sign_test_p'] < 0.05 and k > n / 2) else 'NOT ROBUST (pre-registered)'
json.dump(res, open('results/paxdb_datasets.json', 'w'), indent=1)
print(json.dumps({k_: v for k_, v in res.items() if k_ != 'per_dataset'}, indent=1))
for kf, v in per.items(): print(kf[7:50], v['n_covered'], round(v.get('auroc', float('nan')), 3))
