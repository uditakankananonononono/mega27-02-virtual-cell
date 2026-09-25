"""Post-selection sensitivity: degree-stratified null of false-priority modules."""
import gzip,hashlib,json
from pathlib import Path
import numpy as np,pandas as pd,numba
from scipy.stats import rankdata
ROOT=Path(__file__).resolve().parents[1]
p=ROOT/'data/external/string_links.txt.gz';a=pd.read_csv(ROOT/'results/network_modules_gene_assignments.csv')
assert len(a)==1249 and a.bnumber.is_unique and int(a.priority.sum())==100
ids={x:i for i,x in enumerate(a.bnumber)};edges=set()
with gzip.open(p,'rt') as f:
 next(f)
 for line in f:
  x,y,score=line.split()
  if int(score)<700:continue
  x=x.split('.')[-1];y=y.split('.')[-1]
  if x!=y and x in ids and y in ids:edges.add(tuple(sorted((x,y))))
assert len(edges)==8714
degree=np.zeros(len(a),dtype=np.int32)
for x,y in edges:degree[ids[x]]+=1;degree[ids[y]]+=1
eligible=np.flatnonzero(a.essential.to_numpy()==0)
# Equal-frequency ranks prevent a many-zero-degree tie from collapsing bins.
order=np.lexsort((a.bnumber.to_numpy()[eligible],degree[eligible]))
strata=np.empty(len(eligible),np.int32);strata[order]=np.minimum(np.arange(len(eligible))*10//len(eligible),9)
priority=a.priority.to_numpy().astype(np.int32); assert priority[~(a.essential.to_numpy()==0)].sum()==0
comm=a.community.to_numpy().astype(np.int32)
mods=np.array(sorted(x for x,n in a.community.value_counts().items() if n>=20),dtype=np.int32);assert len(mods)==15 and 3 in mods
# Preserve each stratum's observed priority count; only shuffle non-essential genes.
groups=[eligible[np.flatnonzero(strata==s)].astype(np.int32) for s in range(10)]
counts=[int(priority[g].sum()) for g in groups]
assert sum(counts)==100 and all(k>0 for k in counts)
B=20000; rng=np.random.default_rng(173);draws=np.zeros((B,len(mods)),np.int16)
@numba.njit(cache=True)
def count_members(picks,communities,modules):
 out=np.zeros(len(modules),np.int16)
 for i in picks:
  for j in range(len(modules)):
   if communities[i]==modules[j]:out[j]+=1;break
 return out
# Warm JIT, then verify a pure-NumPy reference on the same draw.
probe=np.concatenate([rng.choice(g,k,replace=False) for g,k in zip(groups,counts)]).astype(np.int32)
first=count_members(probe,comm,mods)
ref=np.array([np.sum(comm[probe]==x) for x in mods]);assert np.array_equal(first,ref)
for b in range(B):
 chosen=np.concatenate([rng.choice(g,k,replace=False) for g,k in zip(groups,counts)]).astype(np.int32)
 draws[b]=count_members(chosen,comm,mods)
obs=np.array([int(priority[comm==x].sum()) for x in mods])
mu=draws.mean(0);sd=draws.std(0,ddof=1)
z=(obs-mu)/sd;zs=(draws-mu)/sd
mx=zs.max(1)
rows=[]
for j,x in enumerate(mods):
 rows.append({'community':int(x),'size':int((comm==x).sum()),'observed_priority':int(obs[j]),
              'degree_stratified_expected':float(mu[j]),'degree_stratified_sd':float(sd[j]),
              'empirical_one_sided_p':float((1+np.sum(draws[:,j]>=obs[j]))/(B+1)),
              'max_statistic_fwer_p':float((1+np.sum(mx>=z[j]))/(B+1))})
out={'design':'post-selection degree-stratified sensitivity on 15 pre-existing STRING/Leiden modules',
     'n_genes':len(a),'n_nonessential':len(eligible),'n_priorities':int(priority.sum()),'n_edges':len(edges),
     'n_modules_tested':len(mods),'n_degree_strata':len(groups),'stratum_sizes':[len(g) for g in groups],
     'stratum_priority_counts':counts,'permutations':B,'random_seed':173,'numba_version':numba.__version__,
     'numba_kernel_executed':True,'reference_match_on_probe':bool(np.array_equal(first,ref)),
     'rows':rows,'module3':next(x for x in rows if x['community']==3),
     'source_sha256':{str(x.relative_to(ROOT)):hashlib.sha256(x.read_bytes()).hexdigest() for x in [p,ROOT/'results/network_modules_gene_assignments.csv']},
     'caveat':'Exploratory post-selection sensitivity, not independent validation. Degree matching controls one network confound only; modules, labels and graph are reused.'}
q=ROOT/'results/network_degree_null.json';q.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'module3':out['module3'],'n_modules_tested':len(mods),'stratum_priority_counts':counts},indent=2))
