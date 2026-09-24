"""Extended real-data figures for the 50-page paper (all from results/)."""
import json, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from sklearn.calibration import calibration_curve
plt.rcParams['font.family'] = 'serif'; plt.rcParams['font.serif'] = ['Times New Roman', 'DejaVu Serif']
al = pd.read_csv('results/aligned_predictions_all_models.csv')
ens = pd.read_csv('results/ensemble_oof.csv')
al = al.merge(ens[['bnumber', 'ensemble']], on='bnumber')
y = al.essential.values
# fig6: score distributions by class
fig, ax = plt.subplots(1, 4, figsize=(13, 3.2))
for a, c in zip(ax, ['gnn', 'cnn', 'kmer', 'ensemble']):
    a.hist(al[c][y == 0], bins=30, alpha=.6, density=True, label='non-essential')
    a.hist(al[c][y == 1], bins=30, alpha=.6, density=True, label='essential')
    a.set_title(c); a.set_xlabel('score')
ax[0].legend(fontsize=7); plt.tight_layout(); plt.savefig('figures/fig6_score_distributions.png', dpi=160); plt.close()
# fig7: calibration
fig, a = plt.subplots(figsize=(5, 4.2))
for c in ['gnn', 'cnn', 'kmer', 'ensemble']:
    s = al[c].values; s = (s - s.min()) / (s.max() - s.min() + 1e-12)
    pt, pp = calibration_curve(y, s, n_bins=10, strategy='quantile')
    a.plot(pp, pt, 'o-', label=c)
a.plot([0, 1], [0, 1], 'k--', lw=.8); a.set_xlabel('mean predicted (min-max scaled)'); a.set_ylabel('observed essential fraction')
a.legend(); plt.tight_layout(); plt.savefig('figures/fig7_calibration.png', dpi=160); plt.close()
# fig8: per-fold AUROC
g = json.load(open('results/gnn_results.json')); sc = json.load(open('results/seqcnn_results.json'))
fig, a = plt.subplots(figsize=(5.5, 3.4))
gf = [f['auroc'] for f in g['folds']]; cf = [f['auroc'] for f in sc['cnn']['folds']]
x = np.arange(3); a.bar(x - .18, gf, .36, label='GCN'); a.bar(x + .18, cf, .36, label='CNN')
a.set_xticks(x); a.set_xticklabels(['fold 0', 'fold 1', 'fold 2']); a.set_ylim(.5, .75); a.set_ylabel('AUROC'); a.legend()
plt.tight_layout(); plt.savefig('figures/fig8_fold_auroc.png', dpi=160); plt.close()
# fig9: threshold sensitivity
t = json.load(open('results/threshold_sensitivity.json'))
th = sorted(float(k) for k in t); f1 = [t[str(k) if str(k) in t else repr(k)]['f1'] for k in th]
n = [t[str(k) if str(k) in t else repr(k)]['n_predicted_essential'] for k in th]
fig, a = plt.subplots(figsize=(5.5, 3.4)); a.semilogx(th, f1, 'o-'); a.set_xlabel('growth-fraction threshold tau'); a.set_ylabel('F1')
b = a.twinx(); b.semilogx(th, n, 's--', color='C1'); b.set_ylabel('# predicted essential', color='C1')
plt.tight_layout(); plt.savefig('figures/fig9_threshold.png', dpi=160); plt.close()
# fig10: precision at top-k
fig, a = plt.subplots(figsize=(5.5, 3.6))
for c in ['fba_min', 'gnn', 'cnn', 'ensemble']:
    o = np.argsort(-al[c].values, kind='stable'); ks = np.arange(10, 400, 10)
    a.plot(ks, [y[o[:k]].mean() for k in ks], label=c)
a.axhline(y.mean(), color='k', ls=':', label='prevalence'); a.set_xlabel('top-k genes'); a.set_ylabel('precision@k'); a.legend(fontsize=8)
plt.tight_layout(); plt.savefig('figures/fig10_precision_at_k.png', dpi=160); plt.close()
print('ok')
