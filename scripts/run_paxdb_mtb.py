"""PaxDb (83332) abundance vs DeJesus 2017 essentiality in M. tuberculosis (prereg notes/prereg_paxdb_mtb.md)."""
import glob, json, math, os
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu, spearmanr, binomtest

D = 'data/external/mtb'
t = pd.read_excel(f'{D}/DeJesus_mbio.xlsx', header=None, skiprows=2)
call = dict(zip(t[0].astype(str), t[12].astype(str)))
lab = {g: 1 if c in ('ES', 'ESD') else 0 for g, c in call.items() if c in ('ES', 'ESD', 'NE')}

def auroc_se(a, n1, n0):
    q1, q2 = a / (2 - a), 2 * a * a / (1 + a)
    return math.sqrt((a * (1 - a) + (n1 - 1) * (q1 - a * a) + (n0 - 1) * (q2 - a * a)) / (n1 * n0))

def load(f):
    ab = {}
    for line in open(f):
        if line.startswith('#') or not line.strip(): continue
        p = line.rstrip('\n').split('\t')
        try: v = float(p[2])
        except (IndexError, ValueError): continue
        if v > 0: ab[p[1].split('.', 1)[1]] = v
    return ab

files = sorted(f for f in glob.glob(f'{D}/83332-*.txt') if 'integrated' not in f)
abs_ = {os.path.basename(f)[6:-4]: load(f) for f in files}
names = list(abs_)
# duplicate rule
dup_of = {}
for i, a in enumerate(names):
    for b in names[:i]:
        if b in dup_of: continue
        sh = set(abs_[a]) & set(abs_[b])
        if len(sh) > 0.9 * min(len(abs_[a]), len(abs_[b])) and len(sh) > 50:
            r = spearmanr([abs_[a][g] for g in sh], [abs_[b][g] for g in sh]).statistic
            if r > 0.99: dup_of[a] = b; break
res = {}
for n, ab in abs_.items():
    gs = [g for g in ab if g in lab]
    y = np.array([lab[g] for g in gs]); x = np.log10([ab[g] for g in gs])
    n1, n0 = int(y.sum()), int(len(y) - y.sum())
    a = float(mannwhitneyu(x[y == 1], x[y == 0]).statistic / (n1 * n0))
    gd = [math.log10(ab[g]) for g in ab if call.get(g) == 'GD']
    res[n] = dict(n_genes=len(gs), n_essential=n1, n_nonessential=n0, auroc=a, se=auroc_se(a, n1, n0),
                  duplicate_of=dup_of.get(n), median_log_ab_ES=float(np.median(x[y == 1])),
                  median_log_ab_GD=float(np.median(gd)) if gd else None, median_log_ab_NE=float(np.median(x[y == 0])))
uniq = {n: d for n, d in res.items() if d['duplicate_of'] is None}
k = sum(d['auroc'] > 0.5 for d in uniq.values()); N = len(uniq)
p = binomtest(k, N, 0.5, alternative='greater').pvalue
w = np.array([1 / d['se'] ** 2 for d in uniq.values()]); m = float((w * [d['auroc'] for d in uniq.values()]).sum() / w.sum()); se = 1 / math.sqrt(w.sum())
out = dict(datasets=res, label_counts=dict(essential=sum(lab.values()), nonessential=len(lab) - sum(lab.values())),
           primary=dict(n_unique=N, n_auroc_above_half=k, sign_test_p=float(p), pooled_auroc=m, ci95=[m - 1.96 * se, m + 1.96 * se],
                        auroc_range=[min(d['auroc'] for d in uniq.values()), max(d['auroc'] for d in uniq.values())],
                        verdict='HOLDS IN MTB' if p < 0.05 else 'NOT SHOWN'),
           secondary_gd_between=sum(1 for d in uniq.values() if d['median_log_ab_GD'] is not None and d['median_log_ab_NE'] <= d['median_log_ab_GD'] <= d['median_log_ab_ES']))
json.dump(out, open('results/paxdb_mtb.json', 'w'), indent=1)
for n, d in res.items(): print(n[:50], d['n_genes'], d['n_essential'], round(d['auroc'], 3), d['duplicate_of'])
print(out['primary'], out['secondary_gd_between'], out['label_counts'])
