"""Exploratory agreement of OOF ranking scores across tested learners."""
import hashlib,itertools,json
from pathlib import Path
import numpy as np,pandas as pd,pingouin as pg
from scipy.stats import spearmanr
v=pd.read_csv('results/ensemble_v2_oof.csv')[['bnumber','v2_lr']]
a=pd.read_csv('results/learner_sensitivity_oof.csv')
d=v.merge(a,on='bnumber',validate='one_to_one');n=len(d);assert n==1249
cols=['v2_lr','lightgbm','catboost','imblearn_oversampled_lr']
rng=np.random.default_rng(11);out=[]
for u,w in itertools.combinations(cols,2):
    r=pg.corr(d[u],d[w],method='spearman')
    bs=np.array([spearmanr(d[u].values[ix],d[w].values[ix]).statistic for ix in (rng.integers(0,n,n) for _ in range(1000))])
    au=set(d.nlargest(50,u).bnumber);aw=set(d.nlargest(50,w).bnumber)
    out.append({'a':u,'b':w,'spearman':float(r['r'].iloc[0]),'pingouin_p':float(r['p_val'].iloc[0]),
                'bootstrap_ci95':[float(np.quantile(bs,.025)),float(np.quantile(bs,.975))],
                'top50_intersection':len(au&aw),'top50_jaccard':len(au&aw)/len(au|aw)})
res={'question':'Agreement of OOF rankings across four learner variants','n_genes':n,'comparisons':out,
     'input_sha256':{p:hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in ['results/ensemble_v2_oof.csv','results/learner_sensitivity_oof.csv']},
     'caveat':'Exploratory; models share genes, training labels and feature construction. Correlation significance is not biological replication.'}
Path('results/score_agreement.json').write_text(json.dumps(res,indent=2)+'\n')
print(json.dumps(res,indent=2))
