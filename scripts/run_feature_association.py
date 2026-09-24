"""Univariate association of each label-free external feature with Gerdes essentiality:
Mann-Whitney U (two-sided), AUROC effect size, Benjamini-Hochberg FDR (statsmodels)."""
import json, pandas as pd, numpy as np
from scipy.stats import mannwhitneyu
from sklearn.metrics import roc_auc_score
from statsmodels.stats.multitest import multipletests
al = pd.read_csv('results/aligned_predictions_all_models.csv')[['bnumber', 'essential']]
ext = pd.read_csv('results/external_features.csv')
df = al.merge(ext, on='bnumber', how='left').fillna(0)
rows = []
for c in [c for c in ext.columns if c != 'bnumber']:
    a, b = df[c][df.essential == 1], df[c][df.essential == 0]
    rows.append({'feature': c, 'auroc': roc_auc_score(df.essential, df[c]), 'median_ess': a.median(),
                 'median_non': b.median(), 'p': mannwhitneyu(a, b).pvalue})
r = pd.DataFrame(rows)
r['q_bh'] = multipletests(r.p, method='fdr_bh')[1]
r = r.sort_values('p'); r.to_csv('results/feature_association.csv', index=False)
print(r.to_string(index=False, float_format=lambda x: f'{x:.3g}'))
