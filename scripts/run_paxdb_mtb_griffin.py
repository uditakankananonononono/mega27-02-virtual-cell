"""H.21 label replication with Griffin 2011 (prereg notes/prereg_paxdb_mtb_griffin.md)."""
import json, math
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu, binomtest
from sklearn.metrics import cohen_kappa_score
import importlib.util
D = 'data/external/mtb'
g = pd.read_excel(f'{D}/griffin2011_table2.xlsx')
lab = {r: int(p < 0.05) for r, p in zip(g['Rv_ID'].astype(str), g['p_val']) if pd.notna(p)}
prev = json.load(open('results/paxdb_mtb.json'))
uniq = [n for n, d in prev['datasets'].items() if d['duplicate_of'] is None]
def load(f):
    ab = {}
    for line in open(f):
        if line.startswith('#') or not line.strip(): continue
        p = line.rstrip('\n').split('\t')
        try: v = float(p[2])
        except (IndexError, ValueError): continue
        if v > 0: ab[p[1].split('.', 1)[1]] = v
    return ab
res = {}
for n in uniq:
    ab = load(f'{D}/83332-{n}.txt')
    gs = [x for x in ab if x in lab]
    y = np.array([lab[x] for x in gs]); x = np.log10([ab[z] for z in gs])
    n1, n0 = int(y.sum()), int(len(y) - y.sum())
    res[n] = dict(n_genes=len(gs), n_essential=n1, auroc=float(mannwhitneyu(x[y == 1], x[y == 0]).statistic / (n1 * n0)))
k = sum(d['auroc'] > 0.5 for d in res.values()); N = len(res)
p = float(binomtest(k, N, 0.5, alternative='greater').pvalue)
dj = pd.read_excel(f'{D}/DeJesus_mbio.xlsx', header=None, skiprows=2)
djl = {r: int(c in ('ES', 'ESD')) for r, c in zip(dj[0].astype(str), dj[12].astype(str)) if c in ('ES', 'ESD', 'NE')}
both = [r for r in lab if r in djl]
kap = float(cohen_kappa_score([lab[r] for r in both], [djl[r] for r in both]))
out = dict(datasets=res, n_griffin_essential=sum(lab.values()), n_griffin_labelled=len(lab),
           primary=dict(n=N, n_auroc_above_half=k, sign_test_p=p, median_auroc=float(np.median([d['auroc'] for d in res.values()])),
                        verdict='REPLICATES' if p < 0.05 else 'NOT SHOWN'),
           label_agreement=dict(n_both=len(both), cohen_kappa=kap))
json.dump(out, open('results/paxdb_mtb_griffin.json', 'w'), indent=1)
for n, d in res.items(): print(n[:45], d['n_genes'], d['n_essential'], round(d['auroc'], 3))
print(out['primary'], out['label_agreement'], out['n_griffin_essential'])
