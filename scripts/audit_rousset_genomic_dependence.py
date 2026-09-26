"""Post-result sensitivity of Rousset AUROC difference to local gene dependence.
B-number genomic windows are proxies, NOT verified operons or independent assays.
"""
import json
import numpy as np,pandas as pd
from sklearn.metrics import roc_auc_score
r=pd.read_csv('data/external/rousset2018/pgen.1007749.s012.csv')
u=pd.read_csv('data/external/uniprot_ecoli.tsv',sep='\t')
n2b={}
for _,x in u.iterrows():
 b=[t for t in str(x['Gene Names (ordered locus)']).split() if t.startswith('b') and t[1:].isdigit()]
 if b and isinstance(x['Gene Names (primary)'],str):n2b.setdefault(x['Gene Names (primary)'],b[0])
r['bnumber']=r.gene.map(n2b)
v=pd.read_csv('results/ensemble_v2_oof.csv')
a=pd.read_csv('results/aligned_predictions_all_models.csv')[['bnumber','fba_min']]
d=r.dropna(subset=['bnumber']).merge(v[['bnumber','v2_lr']],on='bnumber').merge(a,on='bnumber')
assert len(d)==1214 and d.bnumber.nunique()==1214
idx=d.bnumber.str.extract(r'^b(\d+)$')[0].astype(int).to_numpy()
y=(d.median_coding<=-5).astype(int).to_numpy(); s=d.v2_lr.to_numpy(); f=d.fba_min.to_numpy()
res={'status':'descriptive post-result genomic-neighborhood dependence sensitivity; no new preregistered gate',
     'n_genes':len(d),'positives':int(y.sum()),'score_difference':float(roc_auc_score(y,s)-roc_auc_score(y,f)),
     'window_definition':'contiguous fixed-size b-number ID bins; surrogate for local dependence, NOT operon annotation',
     'window_bootstrap':[]}
for width in [5,10,20,50]:
 groups=idx//width; unique=np.unique(groups); members={z:np.flatnonzero(groups==z) for z in unique}
 rng=np.random.default_rng(4100+width); draws=[]
 for _ in range(3000):
  chosen=rng.choice(unique,size=len(unique),replace=True)
  inds=np.concatenate([members[z] for z in chosen]); yy=y[inds]
  if yy.min()==yy.max():continue
  draws.append(float(roc_auc_score(yy,s[inds])-roc_auc_score(yy,f[inds])))
 res['window_bootstrap'].append({'width_bnumber_ids':width,'n_occupied_windows':len(unique),'n_draws_valid':len(draws),
                                'diff_ci95':[float(z) for z in np.percentile(draws,[2.5,97.5])],
                                'frac_diff_le_zero':float(np.mean(np.array(draws)<=0))})
res['limits']='No validated operon annotations, no donor/assay batch groups, no refit on external labels. The test remains ranking of a Gerdes-trained score in another assay, not independent model training or Gerdes G2 repair.'
with open('results/rousset_genomic_dependence_audit.json','w') as fp:json.dump(res,fp,indent=2);fp.write('\n')
print(json.dumps(res,indent=2))
