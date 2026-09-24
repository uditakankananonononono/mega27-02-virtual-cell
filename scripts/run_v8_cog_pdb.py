"""v8: COG category/spread + PDB coverage (notes/prereg_cog_pdb.md)."""
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
cog = pd.read_csv('data/external/cog/ecoli_k12_cog.csv', header=None)
cdef = pd.read_csv('data/external/cog/cog-20.def.tab', sep='\t', header=None, encoding='latin-1')
spread = dict(pd.read_csv('data/external/cog/cog_genome_spread.tsv', sep='\t', header=None).values)
c2f = dict(zip(cdef[0], cdef[1])); LET = list('JKLCEFGHIPQMDOTVSR')
g2c = cog.groupby(0)[6].apply(list).to_dict()
F = []
for b in df.bnumber:
    cs = g2c.get(b, []); lets = set(''.join(c2f.get(c, '') for c in cs))
    F.append([int(l in lets) for l in LET] + [int(bool(lets - set(LET))), int(not cs), np.log1p(max([spread.get(c, 0) for c in cs], default=0))])
cogX = np.array(F, float); COLS = [f'cog_{l}' for l in LET] + ['cog_other', 'no_cog', 'log_cog_spread']
up = pd.read_csv('data/external/uniprot_ecoli.tsv', sep='\t'); b2u = {}
for _, x in up.iterrows():
    for t in str(x['Gene Names (ordered locus)']).split():
        if re.fullmatch(r'b\d{4}', t): b2u.setdefault(t, x['Entry'])
pdb = pd.read_csv('data/external/uniprot_ecoli_pdb.tsv', sep='\t').fillna('')
u2n = {r.Entry: len([p for p in r.PDB.split(';') if p]) for r in pdb.itertuples()}
df['log_pdb'] = np.log1p([u2n.get(b2u.get(b), 0) for b in df.bnumber])
for i, c in enumerate(COLS): df[c] = cogX[:, i]
old = df[base + extc].values.astype(float)
NEWC = COLS + ['log_pdb']; NOPDB = COLS
def fold_feats(tr, te, cols=None):
    cols = cols or NEWC
    return df[cols].values[tr], df[cols].values[te]
lr = lambda: make_pipeline(StandardScaler(), LogisticRegression(max_iter=4000, class_weight='balanced', C=0.5))
v2 = np.zeros(len(y)); v8 = np.zeros(len(y)); new_only = np.zeros(len(y)); v8nopdb = np.zeros(len(y))
for tr, te in StratifiedKFold(3, shuffle=True, random_state=7).split(old, y):
    ntr, nte = fold_feats(tr, te)
    v2[te] = lr().fit(old[tr], y[tr]).predict_proba(old[te])[:, 1]
    v8[te] = lr().fit(np.hstack([old[tr], ntr]), y[tr]).predict_proba(np.hstack([old[te], nte]))[:, 1]
    a_, b_ = fold_feats(tr, te, NOPDB); v8nopdb[te] = lr().fit(np.hstack([old[tr], a_]), y[tr]).predict_proba(np.hstack([old[te], b_]))[:, 1]
    new_only[te] = lr().fit(ntr, y[tr]).predict_proba(nte)[:, 1]
ref = pd.read_csv('results/ensemble_v2_oof.csv').set_index('bnumber').loc[df.bnumber, 'v2_lr'].values
def paired(yy, a, b, f, n=2000):
    rng = np.random.default_rng(11); d = []
    for _ in range(n):
        i = rng.integers(0, len(yy), len(yy))
        if 0 < yy[i].sum() < len(i): d.append(f(yy[i], a[i]) - f(yy[i], b[i]))
    return {'diff': float(f(yy, a) - f(yy, b)), 'ci95': [float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))]}
res = {'n_genes': int(len(y)), 'cog_coverage': int((df.no_cog == 0).sum()), 'pdb_nonzero': int((df.log_pdb > 0).sum()),
       'v2_lr_reproduced_max_abs_dev_vs_committed': float(np.abs(v2 - ref).max()),
       'auroc': {'v2_lr': roc_auc_score(y, v2), 'v8_lr': roc_auc_score(y, v8), 'new_features_only': roc_auc_score(y, new_only)},
       'auprc': {'v2_lr': average_precision_score(y, v2), 'v8_lr': average_precision_score(y, v8)},
       'primary_auroc_v8_minus_v2': paired(y, v8, v2, roc_auc_score),
       'secondary_auprc_v8_minus_v2': paired(y, v8, v2, average_precision_score),
       'secondary_v8_noPDB_minus_v2': paired(y, v8nopdb, v2, roc_auc_score),
       'univariate_auroc': {c: roc_auc_score(y, df[c].fillna(df[c].median())) for c in ['log_cog_spread', 'log_pdb']}}
r = pd.read_csv('data/external/rousset2018/pgen.1007749.s012.csv')
up = pd.read_csv('data/external/uniprot_ecoli.tsv', sep='\t'); n2b = {}
for _, x in up.iterrows():
    b = [t for t in str(x['Gene Names (ordered locus)']).split() if re.fullmatch(r'b\d{4}', t)]
    if b and isinstance(x['Gene Names (primary)'], str): n2b.setdefault(x['Gene Names (primary)'], b[0])
r['bnumber'] = r.gene.map(n2b); m = df[['bnumber']].assign(v2=v2, v8=v8).merge(r[['bnumber', 'median_coding']].dropna(), on='bnumber')
yr = (m.median_coding <= -5).astype(int).values
res['secondary_rousset_auroc_v8_minus_v2'] = paired(yr, m.v8.values, m.v2.values, roc_auc_score)
res['verdict'] = 'IMPROVES (pre-registered)' if res['primary_auroc_v8_minus_v2']['ci95'][0] > 0 else 'NO SIGNIFICANT IMPROVEMENT (pre-registered null)'
json.dump(res, open('results/v8_cog_pdb.json', 'w'), indent=1); print(json.dumps(res, indent=1))
pd.DataFrame({'bnumber': df.bnumber, 'essential': y, 'v8_lr': v8, 'new_features_only': new_only}).to_csv('results/v8_oof.csv', index=False)
