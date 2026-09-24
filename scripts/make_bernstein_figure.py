import json, glob, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'serif'
rows = [json.load(open(f)) for f in sorted(glob.glob('results/bernstein/*.json')) if not f.endswith('_ids.json')]
rows.sort(key=lambda r: r['pr_auc_bernstein_metric'])
fig, a = plt.subplots(figsize=(6.4, 0.5 + 0.45 * len(rows)))
a.barh([r['variant'] for r in rows], [r['pr_auc_bernstein_metric'] for r in rows], color='C0')
for i, r in enumerate(rows): a.text(r['pr_auc_bernstein_metric'] + 0.005, i, f"{r['pr_auc_bernstein_metric']:.3f}", va='center', fontsize=8)
a.set_xlabel('precision-recall AUC (Bernstein et al. 2023 metric)'); a.set_xlim(0.5, 1.0)
plt.tight_layout(); plt.savefig('figures/fig11_bernstein_benchmark.png', dpi=160); print(len(rows))
