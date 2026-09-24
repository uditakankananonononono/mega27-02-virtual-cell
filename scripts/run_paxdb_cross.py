"""PaxDb x Fitness Browser: abundance vs putative essentiality outside E. coli (prereg notes/prereg_paxdb_cross.md)."""
import gzip, glob, json, math, os
import numpy as np
from scipy.stats import spearmanr, binomtest, mannwhitneyu
from sklearn.linear_model import LogisticRegression

D = 'data/external/paxdb_cross'
ORGS = {'160488': 'Putida', '1140': 'SynE', '211586': 'MR1', '226186': 'Btheta', '882': 'DvH'}

def read_fasta(fh):
    name, seq = None, []
    for line in fh:
        line = line.strip()
        if line.startswith('>'):
            if name: yield name, ''.join(seq)
            name, seq = line[1:].split()[0], []
        else: seq.append(line)
    if name: yield name, ''.join(seq)

fb_seq = {}
with gzip.open('data/cross/FB_aaseqs.gz', 'rt') as fh:
    for n, s in read_fasta(fh):
        org, loc = n.split(':', 1)
        fb_seq.setdefault(org, {})[loc] = s.rstrip('*').upper()

def auroc_se(a, n1, n0):
    q1, q2 = a / (2 - a), 2 * a * a / (1 + a)
    return math.sqrt((a * (1 - a) + (n1 - 1) * (q1 - a * a) + (n0 - 1) * (q2 - a * a)) / (n1 * n0))

out = {'datasets': {}, 'per_organism': {}}
for taxid, org in ORGS.items():
    seq2loc = {}
    for loc, s in fb_seq[org].items(): seq2loc.setdefault(s, []).append(loc)
    with gzip.open(f'{D}/{taxid}.seq.fa.gz', 'rt') as fh:
        str2loc = {}
        for n, s in read_fasta(fh):
            locs = seq2loc.get(s.rstrip('*').upper(), [])
            if len(locs) == 1: str2loc[n] = locs[0]
    fit = {}
    for line in open(f'data/cross/{org}_gene_median_fitness.tsv').readlines()[1:]:
        loc, _, mf = line.rstrip('\n').split('\t'); fit[loc] = float(mf)
    aucs = []
    for f in sorted(glob.glob(f'{D}/{taxid}-*.txt')):
        if 'integrated' in f: continue
        ab = {}
        for line in open(f):
            if line.startswith('#') or not line.strip(): continue
            p = line.rstrip('\n').split('\t')
            try: v = float(p[2])
            except (IndexError, ValueError): continue
            if v > 0 and p[1] in str2loc: ab[str2loc[p[1]]] = v
        locs = list(ab)
        y = np.array([0 if l in fit else 1 for l in locs]); x = np.log10([ab[l] for l in locs])
        n1, n0 = int(y.sum()), int(len(y) - y.sum())
        a = mannwhitneyu(x[y == 1], x[y == 0]).statistic / (n1 * n0) if n1 and n0 else float('nan')
        ne = [l for l in locs if l in fit]
        rho = spearmanr([ab[l] for l in ne], [fit[l] for l in ne]).statistic
        L = np.log10([len(fb_seq[org][l]) for l in locs])
        X = np.column_stack([(x - x.mean()) / x.std(), (L - L.mean()) / L.std()])
        coef = LogisticRegression(C=1e6, max_iter=1000).fit(X, y).coef_[0]
        name = os.path.basename(f)[:-4]
        out['datasets'][name] = dict(organism=org, n_mapped_detected=len(locs), n_essential=n1, n_nonessential=n0,
            auroc=a, se=auroc_se(a, n1, n0), spearman_abund_vs_fitness=rho,
            logit_coef_abundance=float(coef[0]), logit_coef_length=float(coef[1]))
        aucs.append((a, auroc_se(a, n1, n0)))
    w = np.array([1 / s ** 2 for _, s in aucs]); m = float(np.sum(w * [a for a, _ in aucs]) / w.sum())
    out['per_organism'][org] = dict(n_datasets=len(aucs), n_string_mapped=len(str2loc), pooled_auroc=m, pooled_se=float(1 / math.sqrt(w.sum())))
ds = out['datasets'].values()
k = sum(d['auroc'] > 0.5 for d in ds); n = len(out['datasets'])
p = binomtest(k, n, 0.5, alternative='greater').pvalue
out['primary'] = dict(n_datasets=n, n_auroc_above_half=k, sign_test_p=p, verdict='GENERALISES' if p < 0.05 else 'NOT SHOWN')
out['secondary'] = dict(n_spearman_negative=sum(d['spearman_abund_vs_fitness'] < 0 for d in ds),
    n_abundance_coef_positive_after_length=sum(d['logit_coef_abundance'] > 0 for d in ds), n=n)
json.dump(out, open('results/paxdb_cross.json', 'w'), indent=1, default=float)
for k2, d in out['datasets'].items(): print(k2[:45], d['organism'], d['n_mapped_detected'], d['n_essential'], round(d['auroc'], 3), round(d['spearman_abund_vs_fitness'], 3), round(d['logit_coef_abundance'], 2))
print(out['per_organism']); print(out['primary'], out['secondary'])
