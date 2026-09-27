"""Source-operon exploratory dependence sensitivity, per locked notes protocol."""
import hashlib,json
from pathlib import Path
import numpy as np,pandas as pd
from sklearn.metrics import roc_auc_score
F=['data/external/rousset2018/pgen.1007749.s012.csv','data/external/uniprot_ecoli.tsv','results/ensemble_v2_oof.csv','results/aligned_predictions_all_models.csv']
hashes={p:hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in F}
r=pd.read_csv(F[0]);u=pd.read_csv(F[1],sep='\t');v=pd.read_csv(F[2]);a=pd.read_csv(F[3])[['bnumber','fba_min']]
n2b={}
for _,x in u.iterrows():
 b=[z for z in str(x['Gene Names (ordered locus)']).split() if z.startswith('b') and z[1:].isdigit()]
 if b and isinstance(x['Gene Names (primary)'],str):n2b.setdefault(x['Gene Names (primary)'],b[0])
r['bnumber']=r.gene.map(n2b)
d=r.dropna(subset=['bnumber']).merge(v[['bnumber','v2_lr']],on='bnumber',validate='one_to_one').merge(a,on='bnumber',validate='one_to_one').copy()
assert len(d)==1214 and d.bnumber.nunique()==1214
s=d.operon.fillna('').astype(str).str.strip(); missing=(s=='')
d['source_operon_group']=np.where(missing,'missing:'+d.bnumber,s)
y=(d.median_coding<=-5).astype(int).to_numpy();assert y.sum()==70
model=d.v2_lr.to_numpy();base=d.fba_min.to_numpy();groups=d.groupby('source_operon_group',sort=True).indices
keys=sorted(groups); mem=[np.asarray(groups[k]) for k in keys];n=len(mem)
assert sum(len(z) for z in mem)==len(d)
def diff(ids):return float(roc_auc_score(y[ids],model[ids])-roc_auc_score(y[ids],base[ids]))
def run(seed,one_per_group=False):
 rng=np.random.default_rng(seed);vals=[]
 for _ in range(3000):
  if one_per_group:ids=np.array([rng.choice(z) for z in mem])
  else:ids=np.concatenate([mem[j] for j in rng.integers(0,n,n)])
  if y[ids].min()==y[ids].max():continue
  vals.append(diff(ids))
 return {'n_draws_valid':len(vals),'diff_ci95':[float(t) for t in np.percentile(vals,[2.5,97.5])],
         'frac_diff_le_zero':float(np.mean(np.array(vals)<=0)),'median_diff':float(np.median(vals))}
res={'status':'post-result exploratory source-operon dependence sensitivity, no new gate',
 'protocol':'notes/prereg_rousset_operon_sensitivity.md','source_sha256':hashes,
 'n_genes':len(d),'n_positives':int(y.sum()),'n_missing_operon':int(missing.sum()),
 'n_source_operon_groups':n,'n_multigene_groups':int(sum(len(z)>1 for z in mem)),
 'n_groups_with_positive':int(sum(y[z].sum()>0 for z in mem)),
 'gene_level_difference':diff(np.arange(len(d))),
 'group_bootstrap':run(271027),'one_gene_per_operon':run(271028,True),
 'limits':'Operon strings are source annotations, not verified expressed/polar units in this assay; same Gerdes-trained scores and Rousset accession, no independent training or wet-lab validation.'}
Path('results/rousset_source_operon_sensitivity.json').write_text(json.dumps(res,indent=2)+'\n')
d[['gene','bnumber','operon','source_operon_group','median_coding','v2_lr','fba_min']].to_csv('results/rousset_source_operon_mapped.csv',index=False)
print(json.dumps({k:v for k,v in res.items() if k!='source_sha256'},indent=2))
