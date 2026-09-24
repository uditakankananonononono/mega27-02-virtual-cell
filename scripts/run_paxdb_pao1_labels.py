"""H.23 label replication with four earlier P. aeruginosa screens (prereg notes/prereg_paxdb_pao1_labels.md)."""
import glob, json, os, re
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu, binomtest
D = 'data/external/pao1'
idn = pd.read_csv(f'{D}/poulsen2019_PA14_id_name.csv')
id2name = {i: str(n).lower() for i, n in zip(idn['pa14_id'], idn['name']) if isinstance(n, str) and re.fullmatch(r'[A-Za-z][A-Za-z0-9]{2,6}', n) and not n.lower().startswith('pa14')}
s6 = pd.read_excel(f'{D}/pnas.1900570116.sd06.xlsx')
screens = ['Turner K.H. et al 2015', 'Lee S.A. et al 2015', 'Skurnik D. et al 2013', 'Liberati N.T. et al 2006']
s6_names = {id2name[i] for i in s6['PA14_ID (for reference)'] if i in id2name}
ab_by_file = {}
for f in sorted(glob.glob(f'{D}/208964-*.txt')):
    ab = {}
    for line in open(f):
        if line.startswith('#') or not line.strip(): continue
        p = line.rstrip('\n').split('\t')
        try: v = float(p[2])
        except (IndexError, ValueError): continue
        n = p[0].lower()
        if v > 0 and not re.fullmatch(r'pa\d{4}.*', n): ab[n] = max(v, ab.get(n, 0))
    ab_by_file[os.path.basename(f)[7:-4]] = ab
out = {'screens': {}}
for sc in screens:
    ess = {id2name[i] for i, c in zip(s6['PA14_ID (for reference)'], s6[sc]) if c == 'E' and i in id2name}
    studies = {}
    for name, ab in ab_by_file.items():
        gs = [g for g in ab if g in ess or g not in s6_names]
        y = np.array([1 if g in ess else 0 for g in gs]); x = np.log10([ab[g] for g in gs])
        n1, n0 = int(y.sum()), int(len(y) - y.sum())
        a = float(mannwhitneyu(x[y == 1], x[y == 0]).statistic / (n1 * n0))
        studies.setdefault('PXD009705' if 'PXD009705' in name else name, []).append(a)
    sa = {s: float(np.mean(v)) for s, v in studies.items()}
    k = sum(v > 0.5 for v in sa.values())
    p = float(binomtest(k, len(sa), 0.5, alternative='greater').pvalue)
    out['screens'][sc] = dict(n_essential_named=len(ess), n_studies=len(sa), n_auroc_above_half=k, sign_test_p=p,
                              median_auroc=float(np.median(list(sa.values()))), verdict='REPLICATES' if p < 0.05 else 'NOT SHOWN')
out['n_replicating'] = sum(v['verdict'] == 'REPLICATES' for v in out['screens'].values())
json.dump(out, open('results/paxdb_pao1_labels.json', 'w'), indent=1)
print(json.dumps(out, indent=1))
