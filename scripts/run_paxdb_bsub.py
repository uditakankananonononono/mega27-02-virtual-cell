"""PaxDb (224308) abundance vs SubtiWiki essential genes in B. subtilis (prereg notes/prereg_paxdb_bsub.md)."""
import glob, json, math, os
import numpy as np
from scipy.stats import mannwhitneyu, binomtest
D = 'data/external/bsub'
ess = {g.lower() for g in json.load(open(f'{D}/subtiwiki_essential_names.json'))}
def auroc_se(a, n1, n0):
    q1, q2 = a / (2 - a), 2 * a * a / (1 + a)
    return math.sqrt((a * (1 - a) + (n1 - 1) * (q1 - a * a) + (n0 - 1) * (q2 - a * a)) / (n1 * n0))
res = {}
for f in sorted(glob.glob(f'{D}/224308-*.txt')):
    ab = {}
    for line in open(f):
        if line.startswith('#') or not line.strip(): continue
        p = line.rstrip('\n').split('\t')
        try: v = float(p[2])
        except (IndexError, ValueError): continue
        name = p[0].lower()
        if v > 0 and not name.startswith('bsu'): ab[name] = max(v, ab.get(name, 0))
    y = np.array([1 if g in ess else 0 for g in ab]); x = np.log10(list(ab.values()))
    n1, n0 = int(y.sum()), int(len(y) - y.sum())
    a = float(mannwhitneyu(x[y == 1], x[y == 0]).statistic / (n1 * n0))
    res[os.path.basename(f)[7:-4]] = dict(n_named=len(ab), n_essential=n1, n_nonessential=n0, auroc=a, se=auroc_se(a, n1, n0))
k = sum(d['auroc'] > 0.5 for d in res.values()); N = len(res)
p = float(binomtest(k, N, 0.5, alternative='greater').pvalue)
w = np.array([1 / d['se'] ** 2 for d in res.values()]); m = float((w * [d['auroc'] for d in res.values()]).sum() / w.sum()); se = 1 / math.sqrt(w.sum())
out = dict(datasets=res, n_essential_list=len(ess), primary=dict(n=N, n_auroc_above_half=k, sign_test_p=p, pooled_auroc=m,
           ci95=[m - 1.96 * se, m + 1.96 * se], verdict='HOLDS IN B. SUBTILIS' if p < 0.05 else 'NOT SHOWN'))
json.dump(out, open('results/paxdb_bsub.json', 'w'), indent=1)
for n, d in res.items(): print(n, d['n_named'], d['n_essential'], round(d['auroc'], 3))
print(out['primary'])
