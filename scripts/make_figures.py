"""Generate all paper figures from saved results (hermetic)."""
import sys, json
sys.path.insert(0, '.')
import warnings; warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, precision_recall_curve

plt.rcParams.update({'font.size': 9, 'axes.titlesize': 10, 'figure.dpi': 200})

# --- Fig 1: AUROC comparison ---
rows = []
for name, path, key in [
    ('core FBA', 'results/benchmark_core_vs_gerdes.json', None),
    ('iJO1366 FBA (min)', 'results/benchmark_ijo1366_vs_gerdes.json', None),
    ('iJO1366 FBA (rich)', 'results/benchmark_ijo1366_rich_vs_gerdes.json', None)]:
    d = json.load(open(path))
    rows.append((name, d['auroc'], d['auroc_ci']))
ens = json.load(open('results/ensemble_results.json'))
fig, ax = plt.subplots(figsize=(6.2, 3.2))
names = [r[0] for r in rows] + ['GNN (graph)', 'CNN (sequence)', 'k-mer LR', 'ENSEMBLE']
aucs = [r[1] for r in rows] + [0.6596, 0.6430, 0.6281, ens['ensemble']['oof_auroc']]
errs = [(r[1]-r[2][0], r[2][1]-r[1]) for r in rows] + [(0,0)]*4
colors = ['#4C72B0']*3 + ['#55A868']*3 + ['#C44E52']
ax.bar(names, aucs, yerr=np.array(errs).T, capsize=3, color=colors, alpha=0.85)
ax.axhline(0.5, ls='--', c='gray', lw=0.8, label='chance')
ax.set_ylabel('AUROC vs Gerdes 2003 truth'); ax.set_ylim(0.4, 0.8)
ax.set_title('Essentiality prediction: individual layers vs stacked ensemble')
plt.setp(ax.get_xticklabels(), rotation=20, ha='right')
ax.legend(); fig.tight_layout(); fig.savefig('figures/fig1_auroc_comparison.png'); plt.close(fig)

# --- Fig 2: ROC + PR on aligned set ---
df = pd.read_csv('results/aligned_predictions_all_models.csv')
ens_oof = pd.read_csv('results/ensemble_oof.csv')
df = df.merge(ens_oof[['bnumber','ensemble']], on='bnumber')
y = df['essential'].to_numpy()
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.2, 3.0))
for col, lab in [('fba_min','FBA minimal'), ('gnn','GNN'), ('cnn','CNN'), ('ensemble','Ensemble')]:
    fpr, tpr, _ = roc_curve(y, df[col]); ax1.plot(fpr, tpr, label=lab, lw=1.2)
    prec, rec, _ = precision_recall_curve(y, df[col]); ax2.plot(rec, prec, label=lab, lw=1.2)
ax1.plot([0,1],[0,1], ls='--', c='gray', lw=0.8); ax1.set_xlabel('FPR'); ax1.set_ylabel('TPR'); ax1.set_title('ROC'); ax1.legend(fontsize=7)
ax2.axhline(y.mean(), ls='--', c='gray', lw=0.8); ax2.set_xlabel('Recall'); ax2.set_ylabel('Precision'); ax2.set_title('Precision-recall'); ax2.legend(fontsize=7)
fig.tight_layout(); fig.savefig('figures/fig2_roc_pr.png'); plt.close(fig)

# --- Fig 3: diauxie, unregulated vs regulated ---
fig, axes = plt.subplots(1, 2, figsize=(6.2, 3.0), sharey=False)
for ax, path, title in [(axes[0], 'results/dfba_diauxie_trajectory.csv', 'Unregulated dFBA (co-utilization)'),
                        (axes[1], 'results/dfba_diauxie_regulated_trajectory.csv', 'With catabolite-repression rule')]:
    t = pd.read_csv(path)
    ax.plot(t.time, t.conc_EX_glc__D_e, label='glucose', c='#4C72B0')
    ax.plot(t.time, t.conc_EX_ac_e, label='acetate', c='#C44E52')
    ax2 = ax.twinx(); ax2.plot(t.time, t.biomass, label='biomass', c='k', ls=':'); ax2.set_ylabel('biomass g/L', fontsize=8)
    ax.set_xlabel('time (h)'); ax.set_ylabel('substrate mmol/L'); ax.set_title(title, fontsize=8.5)
    ax.legend(fontsize=7, loc='center right')
fig.suptitle('Dynamic FBA: diauxic shift on glucose + acetate', y=1.02)
fig.tight_layout(); fig.savefig('figures/fig3_diauxie.png', bbox_inches='tight'); plt.close(fig)

# --- Fig 4: discovery scatter ---
feats = ['fba_min','fba_rich','gnn','cnn','kmer']
norm = df[feats].apply(lambda c: (c - c.min())/(c.max()-c.min()+1e-12))
df['consensus'] = norm.mean(axis=1)
A = pd.read_csv('results/discovery_classA.csv'); B = pd.read_csv('results/discovery_classB.csv')
fig, ax = plt.subplots(figsize=(6.2, 3.2))
jitter = np.random.default_rng(1).normal(0, 0.02, len(df))
ax.scatter(df.consensus, df.essential + jitter, s=3, alpha=0.25, c='gray')
ax.scatter(A.consensus, np.zeros(len(A)), s=25, c='#C44E52', marker='v', label='Class A: model-essential, viable (MoCo cluster)')
ax.scatter(B.consensus, np.ones(len(B)), s=25, c='#4C72B0', marker='^', label='Class B: essential, model-missed')
for _, r in A.head(5).iterrows():
    ax.annotate(r.gene_name, (r.consensus, 0), textcoords='offset points', xytext=(0,-14), ha='center', fontsize=7)
ax.set_xlabel('ensemble consensus essentiality score'); ax.set_ylabel('experimental (Gerdes 2003)')
ax.set_yticks([0,1]); ax.set_yticklabels(['nonessential','essential'])
ax.set_title('Consensus vs experiment: disagreement classes'); ax.legend(fontsize=7)
fig.tight_layout(); fig.savefig('figures/fig4_discovery.png'); plt.close(fig)

# --- Fig 5: condition matrix heatmap ---
cm = pd.read_csv('results/core_condition_matrix.csv')
piv = cm.pivot(index='carbon', columns='aerobic', values='growth_h')
fig, ax = plt.subplots(figsize=(4.5, 3.4))
im = ax.imshow(piv.values, cmap='viridis', aspect='auto')
ax.set_xticks([0,1]); ax.set_xticklabels(['anaerobic','aerobic'])
ax.set_yticks(range(len(piv))); ax.set_yticklabels([c.replace('EX_','').replace('_e','') for c in piv.index], fontsize=7)
for i in range(len(piv)):
    for j in range(2):
        ax.text(j, i, f'{piv.values[i,j]:.2f}', ha='center', va='center', fontsize=6.5,
                color='white' if piv.values[i,j] < piv.values.max()/2 else 'black')
ax.set_title('Core-model growth (h$^{-1}$) across carbon sources')
fig.colorbar(im, label='growth h$^{-1}$')
fig.tight_layout(); fig.savefig('figures/fig5_condition_matrix.png'); plt.close(fig)
print('figures done')
