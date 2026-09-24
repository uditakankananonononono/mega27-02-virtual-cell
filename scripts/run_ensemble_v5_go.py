"""v5: nested fold-internal GO error-term feature (notes/prereg_go_condition_feature.md)."""
import gzip, json, collections, numpy as np, pandas as pd
from Bio import SeqIO
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
from goatools.obo_parser import GODag
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, average_precision_score
rec = SeqIO.read('data/U00096.3.gb', 'genbank'); sym2b = {}
for f in rec.features:
    if f.type == 'CDS' and 'locus_tag' in f.qualifiers:
        b = f.qualifiers['locus_tag'][0]
        for key in ('gene', 'gene_synonym'):
            for v in f.qualifiers.get(key, []):
                for s in v.replace(';', ' ').split(): sym2b.setdefault(s.lower(), b)
dag = GODag('data/external/go-basic.obo', load_obsolete=False, prt=None)
assoc = collections.defaultdict(set)
for line in gzip.open('data/external/ecocyc.gaf.gz', 'rt'):
    if line.startswith('!'): continue
    c = line.rstrip('\n').split('\t')
    if 'NOT' in c[3] or c[4] not in dag: continue
    b = sym2b.get(c[2].lower()) or next((sym2b[s.lower()] for s in c[10].split('|') if s.lower() in sym2b), None)
    if b: assoc[b] |= {c[4]} | dag[c[4]].get_all_parents()
al = pd.read_csv('results/aligned_predictions_all_models.csv'); ext = pd.read_csv('results/external_features.csv')
df = al.merge(ext, on='bnumber', how='left').fillna(0); y = df.essential.values
base = ['fba_min', 'fba_rich', 'gnn', 'cnn', 'kmer']; NEW3 = ['pax_log_ppm', 'cai', 'gc3']
extc = [c for c in ext.columns if c != 'bnumber' and not c.startswith('strctx') and c not in NEW3]
cols = base + extc; X = df[cols].values; genes = list(df.bnumber)
lr = lambda: make_pipeline(StandardScaler(), LogisticRegression(max_iter=4000, class_weight='balanced', C=0.5))
def go_terms_from_errors(idx):
    inner = np.zeros(len(idx))
    for a, b_ in StratifiedKFold(3, shuffle=True, random_state=7).split(X[idx], y[idx]):
        inner[b_] = lr().fit(X[idx][a], y[idx][a]).predict_proba(X[idx][b_])[:, 1]
    r = pd.Series(inner).rank(pct=True).values
    fp = set(np.array(genes)[idx][(y[idx] == 0) & (r > 0.9)])
    ann = [g for g in np.array(genes)[idx] if g in assoc]
    cnt = collections.Counter(t for g in ann for t in assoc[g]); cfp = collections.Counter(t for g in ann if g in fp for t in assoc[g])
    nfp = sum(1 for g in ann if g in fp); N = len(ann); terms, ps = [], []
    for t, k in cnt.items():
        if k < 3: continue
        a = cfp.get(t, 0); p = fisher_exact([[a, nfp - a], [k - a, N - nfp - k + a]], alternative='greater')[1]
        terms.append(t); ps.append(p)
    q = multipletests(ps, method='fdr_bh')[1]
    keep = [terms[i] for i in np.argsort(q) if q[i] < 0.05][:10]
    return keep
v2 = np.zeros(len(y)); v5 = np.zeros(len(y)); kept = []
for tr, te in StratifiedKFold(3, shuffle=True, random_state=7).split(X, y):
    v2[te] = lr().fit(X[tr], y[tr]).predict_proba(X[te])[:, 1]
    keep = go_terms_from_errors(tr); kept.append([(t, dag[t].name) for t in keep])
    f = np.array([1.0 if (g in assoc and assoc[g] & set(keep)) else 0.0 for g in genes])[:, None]
    X5 = np.hstack([X, f]); v5[te] = lr().fit(X5[tr], y[tr]).predict_proba(X5[te])[:, 1]
def paired(a, b, n=2000):
    rng = np.random.default_rng(11); d = []
    for _ in range(n):
        i = rng.integers(0, len(y), len(y))
        if y[i].min() == y[i].max(): continue
        d.append(roc_auc_score(y[i], a[i]) - roc_auc_score(y[i], b[i]))
    d = np.array(d); return {'mean': float(d.mean()), 'ci95': [float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))], 'p_le_0': float((d <= 0).mean())}
res = {'v2_lr': {'auroc': roc_auc_score(y, v2), 'auprc': average_precision_score(y, v2)},
       'v5_lr': {'auroc': roc_auc_score(y, v5), 'auprc': average_precision_score(y, v5)},
       'paired_v5_vs_v2': paired(v5, v2), 'terms_kept_per_fold': kept}
res['verdict'] = 'gain' if res['paired_v5_vs_v2']['ci95'][0] > 0 else 'no significant gain (negative)'
json.dump(res, open('results/ensemble_v5_go.json', 'w'), indent=1, default=float); print(json.dumps(res, indent=1, default=float))
