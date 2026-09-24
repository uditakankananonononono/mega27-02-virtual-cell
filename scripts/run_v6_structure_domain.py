"""v6: AlphaFold pLDDT + Pfam features (notes/prereg_structure_domain_features.md)."""
import json, re, collections, numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, average_precision_score
al = pd.read_csv('results/aligned_predictions_all_models.csv'); ext = pd.read_csv('results/external_features.csv')
df = al.merge(ext, on='bnumber', how='left').fillna(0); y = df.essential.values
base = ['fba_min', 'fba_rich', 'gnn', 'cnn', 'kmer']; NEW3 = ['pax_log_ppm', 'cai', 'gc3']
extc = [c for c in ext.columns if c != 'bnumber' and not c.startswith('strctx') and c not in NEW3]
af = pd.read_csv('data/external/alphafold/af_metrics.tsv', sep='\t').drop_duplicates('bnumber')
df = df.merge(af[['bnumber', 'uniprot', 'plddt_mean', 'plddt_frac_vlow']], on='bnumber', how='left')
pf = pd.read_csv('data/external/uniprot_ecoli_pfam.tsv', sep='\t').fillna('')
u2p = {r.Entry: [p for p in r.Pfam.split(';') if p] for r in pf.itertuples()}
fams = [u2p.get(u, []) if isinstance(u, str) else [] for u in df.uniprot]
df['n_pfam'] = [len(f) for f in fams]
old = df[base + extc].values.astype(float)
def fold_feats(tr, te):
    cnt = collections.defaultdict(lambda: [0, 0])
    for i in tr:
        for p in fams[i]: cnt[p][0] += y[i]; cnt[p][1] += 1
    prior = y[tr].mean()
    def te_val(i, own):
        v = []
        for p in fams[i]:
            e, n = cnt[p]
            if own: e -= y[i]; n -= 1
            if n > 0: v.append((e + 1) / (n + 2))
        return max(v) if v else 0.5 * prior
    med = np.nanmedian(df.plddt_mean.values[tr]), np.nanmedian(df.plddt_frac_vlow.values[tr])
    def mk(idx, own):
        pl = np.where(np.isnan(df.plddt_mean.values[idx]), med[0], df.plddt_mean.values[idx])
        fv = np.where(np.isnan(df.plddt_frac_vlow.values[idx]), med[1], df.plddt_frac_vlow.values[idx])
        return np.column_stack([pl, fv, df.n_pfam.values[idx], [te_val(i, own) for i in idx]])
    return mk(tr, True), mk(te, False)
lr = lambda: make_pipeline(StandardScaler(), LogisticRegression(max_iter=4000, class_weight='balanced', C=0.5))
v2 = np.zeros(len(y)); v6 = np.zeros(len(y)); new_only = np.zeros(len(y))
for tr, te in StratifiedKFold(3, shuffle=True, random_state=7).split(old, y):
    ntr, nte = fold_feats(tr, te)
    v2[te] = lr().fit(old[tr], y[tr]).predict_proba(old[te])[:, 1]
    v6[te] = lr().fit(np.hstack([old[tr], ntr]), y[tr]).predict_proba(np.hstack([old[te], nte]))[:, 1]
    new_only[te] = lr().fit(ntr, y[tr]).predict_proba(nte)[:, 1]
ref = pd.read_csv('results/ensemble_v2_oof.csv').set_index('bnumber').loc[df.bnumber, 'v2_lr'].values
def paired(yy, a, b, f, n=2000):
    rng = np.random.default_rng(11); d = []
    for _ in range(n):
        i = rng.integers(0, len(yy), len(yy))
        if 0 < yy[i].sum() < len(i): d.append(f(yy[i], a[i]) - f(yy[i], b[i]))
    return {'diff': float(f(yy, a) - f(yy, b)), 'ci95': [float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))]}
res = {'n_genes': int(len(y)), 'af_coverage': int(df.plddt_mean.notna().sum()), 'pfam_coverage': int((df.n_pfam > 0).sum()),
       'v2_lr_reproduced_max_abs_dev_vs_committed': float(np.abs(v2 - ref).max()),
       'auroc': {'v2_lr': roc_auc_score(y, v2), 'v6_lr': roc_auc_score(y, v6), 'new_features_only': roc_auc_score(y, new_only)},
       'auprc': {'v2_lr': average_precision_score(y, v2), 'v6_lr': average_precision_score(y, v6)},
       'primary_auroc_v6_minus_v2': paired(y, v6, v2, roc_auc_score),
       'secondary_auprc_v6_minus_v2': paired(y, v6, v2, average_precision_score),
       'univariate_auroc': {c: roc_auc_score(y, df[c].fillna(df[c].median())) for c in ['plddt_mean', 'plddt_frac_vlow', 'n_pfam']}}
r = pd.read_csv('data/external/rousset2018/pgen.1007749.s012.csv')
up = pd.read_csv('data/external/uniprot_ecoli.tsv', sep='\t'); n2b = {}
for _, x in up.iterrows():
    b = [t for t in str(x['Gene Names (ordered locus)']).split() if re.fullmatch(r'b\d{4}', t)]
    if b and isinstance(x['Gene Names (primary)'], str): n2b.setdefault(x['Gene Names (primary)'], b[0])
r['bnumber'] = r.gene.map(n2b); m = df[['bnumber']].assign(v2=v2, v6=v6).merge(r[['bnumber', 'median_coding']].dropna(), on='bnumber')
yr = (m.median_coding <= -5).astype(int).values
res['secondary_rousset_auroc_v6_minus_v2'] = paired(yr, m.v6.values, m.v2.values, roc_auc_score)
res['verdict'] = 'IMPROVES (pre-registered)' if res['primary_auroc_v6_minus_v2']['ci95'][0] > 0 else 'NO SIGNIFICANT IMPROVEMENT (pre-registered null)'
json.dump(res, open('results/v6_structure_domain.json', 'w'), indent=1); print(json.dumps(res, indent=1))
pd.DataFrame({'bnumber': df.bnumber, 'essential': y, 'v6_lr': v6, 'new_features_only': new_only}).to_csv('results/v6_oof.csv', index=False)
