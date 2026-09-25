"""Post-selection Louvain sensitivity using scikit-network on saved STRING graph."""
import gzip,hashlib,json
from pathlib import Path
import numpy as np,pandas as pd,sknetwork
from scipy.sparse import coo_matrix
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
from sklearn.metrics import adjusted_rand_score,normalized_mutual_info_score
from sknetwork.clustering import Louvain
ROOT=Path(__file__).resolve().parents[1];edge=ROOT/'data/external/string_links.txt.gz'
df=pd.read_csv(ROOT/'results/network_modules_gene_assignments.csv');assert len(df)==1249 and df.bnumber.is_unique
ids={g:i for i,g in enumerate(df.bnumber)};pairs=set()
with gzip.open(edge,'rt') as f:
 next(f)
 for l in f:
  x,y,score=l.split()
  if int(score)<700:continue
  x=x.split('.')[-1];y=y.split('.')[-1]
  if x!=y and x in ids and y in ids:pairs.add(tuple(sorted((ids[x],ids[y]))))
assert len(pairs)==8714
row=np.array([i for i,j in pairs]+[j for i,j in pairs]);col=np.array([j for i,j in pairs]+[i for i,j in pairs]);adj=coo_matrix((np.ones(len(row)),(row,col)),shape=(len(df),len(df))).tocsr()
new=Louvain(resolution=1,random_state=17,shuffle_nodes=False,return_probs=False,return_aggregate=False).fit_predict(adj)
old=df.community.to_numpy();assert len(new)==len(old)
priority=df.priority.to_numpy(bool);essential=df.essential.to_numpy(bool);assert priority.sum()==100
mods=sorted(set(new));tests=[]
for mod in mods:
 mask=(new==mod);n=int(mask.sum())
 if n<20:continue
 k=int((mask&priority).sum());a=k;b=n-k;c=int(priority.sum())-k;d=len(df)-n-c
 _,p=fisher_exact([[a,b],[c,d]],alternative='greater')
 tests.append({'community':int(mod),'size':n,'false_priority':k,'fisher_p':float(p),
               'overlap_with_leiden3':int((mask&(old==3)).sum())})
assert tests
q=multipletests([r['fisher_p'] for r in tests],method='fdr_bh')[1]
for r,v in zip(tests,q):r['fdr_bh']=float(v)
closest=max(tests,key=lambda x:x['overlap_with_leiden3'])
# Deterministic degree-stratified draws for the top-overlap Louvain module.
eligible=np.flatnonzero(~essential);degree=np.asarray(adj.sum(axis=1)).ravel()
order=np.lexsort((df.bnumber.to_numpy()[eligible],degree[eligible]))
strata=np.empty(len(eligible),np.int32);strata[order]=np.minimum(np.arange(len(eligible))*10//len(eligible),9)
groups=[eligible[np.flatnonzero(strata==i)] for i in range(10)];counts=[int(priority[g].sum()) for g in groups]
mask=(new==closest['community']);B=20000;rng=np.random.default_rng(293)
draws=np.zeros(B,dtype=np.int16)
for i in range(B):
 pick=np.concatenate([rng.choice(g,k,replace=False) for g,k in zip(groups,counts)])
 draws[i]=mask[pick].sum()
closest['degree_stratified_expected']=float(draws.mean())
closest['degree_stratified_one_sided_p']=float((1+np.sum(draws>=closest['false_priority']))/(B+1))
res={'design':'post-selection cross-algorithm sensitivity on same frozen STRING network and label set',
     'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [edge,ROOT/'results/network_modules_gene_assignments.csv']},
     'package':'scikit-network '+sknetwork.__version__,'algorithm':'Louvain, default resolution, seed 17, no shuffle',
     'n_edges':len(pairs),'n_genes':len(df),'n_louvain_modules_all_sizes':len(mods),'n_tested_size_ge20':len(tests),
     'adjusted_rand_index_vs_leiden':float(adjusted_rand_score(old,new)),
     'normalized_mutual_information_vs_leiden':float(normalized_mutual_info_score(old,new)),
     'closest_to_leiden3':closest,'tested_modules':tests,
     'degree_stratified_permutations':B,'degree_stratified_seed':293,
     'caveat':'Same graph, genes and phenotypes; post-selection sensitivity. A second clustering algorithm is not biological replication.'}
(ROOT/'results/graph_algorithm_sensitivity.json').write_text(json.dumps(res,indent=2)+'\n')
pd.DataFrame({'bnumber':df.bnumber,'leiden':old,'louvain':new,'priority':priority}).to_csv(ROOT/'results/graph_algorithm_assignments.csv',index=False)
print(json.dumps({k:res[k] for k in ('n_louvain_modules_all_sizes','n_tested_size_ge20','adjusted_rand_index_vs_leiden','closest_to_leiden3')},indent=2))
