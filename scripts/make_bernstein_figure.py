import json, glob, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'serif'
rows = [json.load(open(f)) for f in sorted(glob.glob('results/bernstein/*.json')) if not f.endswith('_ids.json')]
rows.sort(key=lambda r: r['pr_auc_bernstein_metric'])
BAD = {'iML1515_bernstein_allcorr', 'iML1515_allcorr_selected'}  # loader-bug runs (Appendix H.4)
lab = lambda r: r['variant'] + (' (INVALID)' if r['variant'] in BAD else '')
col = ['0.7' if r['variant'] in BAD else 'C0' for r in rows]
fig, a = plt.subplots(figsize=(6.4, 0.5 + 0.45 * len(rows)))
a.barh([lab(r) for r in rows], [r['pr_auc_bernstein_metric'] for r in rows], color=col)
a.axvline(0.843, ls='--', color='k', lw=0.8); a.text(0.845, -0.6, 'published 0.843', fontsize=7)
for i, r in enumerate(rows): a.text(r['pr_auc_bernstein_metric'] + 0.005, i, f"{r['pr_auc_bernstein_metric']:.3f}", va='center', fontsize=8)
a.set_xlabel('precision-recall AUC (Bernstein et al. 2023 metric)'); a.set_xlim(0.5, 1.0)
plt.tight_layout(); plt.savefig('figures/fig11_bernstein_benchmark.png', dpi=160); print(len(rows))
