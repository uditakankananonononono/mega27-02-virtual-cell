"""Offline STRING high-confidence network module analysis with igraph/Leiden."""
import gzip,hashlib,json
from pathlib import Path
import numpy as np,pandas as pd,igraph as ig,leidenalg as la
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
al=pd.read_csv('results/aligned_predictions_all_models.csv')[['bnumber','essential']]
v=pd.read_csv('results/ensemble_v2_oof.csv')[['bnumber','v2_lr']]
m=al.merge(v,on='bnumber',validate='one_to_one')
priorities=set(m[m.essential==0].nlargest(100,'v2_lr').bnumber)
allowed=set(m.bnumber)
p=Path('data/external/string_links.txt.gz')
edges=set()
with gzip.open(p,'rt') as f:
    next(f)
    for line in f:
        a,b,score=line.split()
        if int(score)<700:continue
        a=a.split('.')[-1];b=b.split('.')[-1]
        if a!=b and a in allowed and b in allowed:edges.add(tuple(sorted((a,b))))
vertices=sorted(allowed);vid={x:i for i,x in enumerate(vertices)}
g=ig.Graph(n=len(vertices),edges=[(vid[a],vid[b]) for a,b in sorted(edges)],directed=False)
part=la.find_partition(g,la.RBConfigurationVertexPartition,resolution_parameter=1,seed=17,n_iterations=-1)
lab=dict(zip(vertices,part.membership))
m['community']=[lab[x] for x in m.bnumber];m['priority']=m.bnumber.isin(priorities)
rows=[]
for k,gp in m.groupby('community'):
    if len(gp)<20:continue
    rec={'community':int(k),'n_genes':int(len(gp)),'n_essential':int(gp.essential.sum()),
         'n_top100_false_priority':int(gp.priority.sum())}
    for key,flag in [('essential','essential'),('false_priority','priority')]:
        inside=int(gp[flag].sum());outside=int(m[flag].sum()-inside)
        odds,pv=fisher_exact([[inside,len(gp)-inside],[outside,len(m)-len(gp)-outside]],alternative='greater')
        rec[key+'_odds']=float(odds);rec[key+'_p']=float(pv)
    rows.append(rec)
for key in ('essential','false_priority'):
    q=multipletests([r[key+'_p'] for r in rows],method='fdr_bh')[1]
    for r,vv in zip(rows,q):r[key+'_bh_q']=float(vv)
result={'design':'exploratory label-blind Leiden modules on saved STRING >=700 graph',
        'n_genes':len(vertices),'n_high_confidence_edges':len(edges),
        'n_communities_all_sizes':len(part),'n_communities_tested_size_ge20':len(rows),
        'communities':rows,'false_priority_enriched_modules':[r['community'] for r in rows if r['false_priority_bh_q']<0.05],
        'essential_enriched_modules':[r['community'] for r in rows if r['essential_bh_q']<0.05],
        'source_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),
        'caveat':'Same 1,249 genes/labels and one STRING network; network modules are not causal cell states or independent datasets.'}
Path('results/network_modules.json').write_text(json.dumps(result,indent=2)+'\n')
m.to_csv('results/network_modules_gene_assignments.csv',index=False)
print(json.dumps(result,indent=2))
