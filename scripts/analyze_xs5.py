"""Amendment 5 analysis (notes/prereg_cross_species.md): pooled within-organism fitness percentiles, on > off,
model-internal demand-flux label, all eligible Fitness Browser bacteria. Reads results/xs5/<org>.json."""
import glob, json, sys, numpy as np
from scipy.stats import mannwhitneyu
per, pooled, excluded = {}, [], {}
def verdict(d, p): return 'REPLICATED' if (p < 0.05 and d > 0) else ('FALSIFIED' if d <= 0 else 'NOT SIGNIFICANT (direction consistent)')
for f in sorted(glob.glob('results/xs5/*.json')):
    r = json.load(open(f)); org = r['org']
    if 'excluded' in r: excluded[org] = r['excluded']; continue
    on = [v['median_fitness'] for v in r['genes'].values() if v['label'] == 'on_pathway' and v['median_fitness'] is not None]
    off = [v['median_fitness'] for v in r['genes'].values() if v['label'] == 'off_pathway' and v['median_fitness'] is not None]
    rec = {'wt': r['wt'], 'n_rescue_pairs': r['n_rescue_pairs'], 'with_fitness_on': len(on), 'with_fitness_off': len(off)}
    if len(on) < 5 or len(off) < 5: rec['verdict'] = 'UNDERPOWERED'
    else:
        d = float(np.median(on) - np.median(off)); p = float(mannwhitneyu(on, off, alternative='greater').pvalue)
        rec.update(median_on=float(np.median(on)), median_off=float(np.median(off)), median_diff=d, mwu_p=p, verdict=verdict(d, p))
    per[org] = rec
    for g, v in r['genes'].items():
        if v['percentile'] is not None: pooled.append((org, g, v['label'], v['percentile'], v['supplements']))
def ptest(items, min_n=10):
    a = [x[3] for x in items if x[2] == 'on_pathway']; b = [x[3] for x in items if x[2] == 'off_pathway']
    out = {'n_on': len(a), 'n_off': len(b), 'n_organisms': len({x[0] for x in items})}
    if len(a) < min_n or len(b) < min_n: out['verdict'] = 'UNDERPOWERED'; return out
    d = float(np.median(a) - np.median(b)); p = float(mannwhitneyu(a, b, alternative='greater').pvalue)
    out.update(median_pct_on=float(np.median(a)), median_pct_off=float(np.median(b)), diff=d, mwu_p=p, verdict=verdict(d, p))
    return out
res = {'n_organisms_analysed': len(per), 'excluded': excluded, 'primary_pooled_percentile': ptest(pooled),
       'secondary_pooled_excluding_SAM': ptest([x for x in pooled if x[4] != ['amet_c']]),
       'per_organism_verdict_counts': {k: sum(v['verdict'] == k for v in per.values()) for k in {v['verdict'] for v in per.values()}},
       'per_organism': per, 'complete': len(sys.argv) > 1 and sys.argv[1] == '--final'}
out = 'results/cross_species_xs5.json' if res['complete'] else '/tmp/xs5_interim.json'
json.dump(res, open(out, 'w'), indent=1)
print(json.dumps({k: res[k] for k in ('n_organisms_analysed', 'excluded', 'primary_pooled_percentile', 'secondary_pooled_excluding_SAM', 'per_organism_verdict_counts')}, indent=1))
