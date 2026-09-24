"""PaxDb (208964) abundance vs Poulsen 2019 core essential genes (prereg notes/prereg_paxdb_pao1.md)."""
import glob, json, math, os, re
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu, binomtest
D = 'data/external/pao1'
t = pd.read_excel(f'{D}/pnas.1900570116.sd05.xlsx', sheet_name='Essential Genes')
t = t[t['Gene Name'].notna()]
core = {str(n).lower() for n, c in zip(t['Gene Name'], t['Essential Category']) if str(c).strip() == 'Core'}
listed = {str(n).lower() for n in t['Gene Name']}
res = {}
for f in sorted(glob.glob(f'{D}/208964-*.txt')):
    ab = {}
    for line in open(f):
        if line.startswith('#') or not line.strip(): continue
        p = line.rstrip('\n').split('\t')
        try: v = float(p[2])
        except (IndexError, ValueError): continue
        n = p[0].lower()
        if v > 0 and not re.fullmatch(r'pa\d{4}.*', n): ab[n] = max(v, ab.get(n, 0))
    gs = [g for g in ab if g in core or g not in listed]
    y = np.array([1 if g in core else 0 for g in gs]); x = np.log10([ab[g] for g in gs])
    n1, n0 = int(y.sum()), int(len(y) - y.sum())
    a = float(mannwhitneyu(x[y == 1], x[y == 0]).statistic / (n1 * n0))
    name = os.path.basename(f)[7:-4]
    study = 'PXD009705' if 'PXD009705' in name else name
    res[name] = dict(study=study, n_named=len(gs), n_essential=n1, n_nonessential=n0, auroc=a)
studies = {}
for d in res.values(): studies.setdefault(d['study'], []).append(d['auroc'])
sa = {s: float(np.mean(v)) for s, v in studies.items()}
k = sum(v > 0.5 for v in sa.values()); N = len(sa)
p = float(binomtest(k, N, 0.5, alternative='greater').pvalue)
out = dict(datasets=res, study_auroc=sa, n_core_named=len(core), primary=dict(n_studies=N, n_auroc_above_half=k, sign_test_p=p,
           median_study_auroc=float(np.median(list(sa.values()))), auroc_range=[min(sa.values()), max(sa.values())],
           verdict='HOLDS IN P. AERUGINOSA' if p < 0.05 else 'NOT SHOWN'))
json.dump(out, open('results/paxdb_pao1.json', 'w'), indent=1)
for n, d in res.items(): print(n[:45], d['n_named'], d['n_essential'], round(d['auroc'], 3))
print(len(core), out['primary'])
