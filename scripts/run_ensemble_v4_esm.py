"""v4: v2_lr + fold-internal PCA of ESM-2 embeddings (pre-registered in notes/prereg_esm_v4.md)."""
import json, numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.pipeline import make_pipeline
from sklearn.compose import ColumnTransformer
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, average_precision_score
al = pd.read_csv('results/aligned_predictions_all_models.csv')
ext = pd.read_csv('results/external_features.csv')
df = al.merge(ext, on='bnumber', how='left').fillna(0)
E = np.load('results/esm2_t6_embeddings.npy'); eid = json.load(open('results/esm2_t6_ids.json'))
emb = pd.DataFrame(E, columns=[f'esm{i}' for i in range(E.shape[1])]); emb['bnumber'] = eid
emb = emb.drop_duplicates('bnumber')
df = df.merge(emb, on='bnumber', how='left')
cover = float(df['esm0'].notna().mean()); df = df.fillna(0)
y = df.essential.values
base = ['fba_min', 'fba_rich', 'gnn', 'cnn', 'kmer']
NEW3 = ['pax_log_ppm', 'cai', 'gc3']
extc = [c for c in ext.columns if c != 'bnumber' and not c.startswith('strctx') and c not in NEW3]
esmc = [f'esm{i}' for i in range(E.shape[1])]
def oof(make, cols, seed=7):
    X = df[cols]; out = np.zeros(len(y))
    for tr, te in StratifiedKFold(3, shuffle=True, random_state=seed).split(X, y):
        out[te] = make().fit(X.iloc[tr], y[tr]).predict_proba(X.iloc[te])[:, 1]
    return out
LR = lambda: LogisticRegression(max_iter=4000, class_weight='balanced', C=0.5)
lr = lambda: make_pipeline(StandardScaler(), LR())
def v4(k):
    return lambda: make_pipeline(ColumnTransformer([('b', StandardScaler(), base + extc),
                                                    ('e', make_pipeline(StandardScaler(), PCA(k, random_state=0)), esmc)]),
                                 StandardScaler(), LR())
S = {'v2_lr': oof(lr, base + extc), 'v4_lr': oof(v4(16), base + extc + esmc),
     'v4_pca32_lr': oof(v4(32), base + extc + esmc),
     'esm_only_lr': oof(lambda: make_pipeline(StandardScaler(), PCA(16, random_state=0), LR()), esmc)}
res = {'esm_coverage': cover, 'n_genes': int(len(y))}
for k, s in S.items():
    res[k] = {'auroc': roc_auc_score(y, s), 'auprc': average_precision_score(y, s)}
def paired(a, b, n=2000):
    rng = np.random.default_rng(11); d = []
    for _ in range(n):
        i = rng.integers(0, len(y), len(y))
        if y[i].min() == y[i].max(): continue
        d.append(roc_auc_score(y[i], a[i]) - roc_auc_score(y[i], b[i]))
    d = np.array(d); return {'mean': float(d.mean()), 'ci95': [float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))], 'p_le_0': float((d <= 0).mean())}
res['paired_v4lr_vs_v2lr'] = paired(S['v4_lr'], S['v2_lr'])
res['verdict'] = 'gain' if res['paired_v4lr_vs_v2lr']['ci95'][0] > 0 else 'no significant gain (negative)'
json.dump(res, open('results/ensemble_v4_esm.json', 'w'), indent=2, default=float)
pd.DataFrame({'bnumber': df.bnumber, 'essential': y, **S}).to_csv('results/ensemble_v4_oof.csv', index=False)
print(json.dumps(res, indent=1, default=float))
