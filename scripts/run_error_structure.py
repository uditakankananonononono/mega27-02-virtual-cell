"""Unsupervised feature-space structure of out-of-fold false priorities."""
import hashlib,json, warnings
from pathlib import Path
import numpy as np,pandas as pd
from sklearn.preprocessing import StandardScaler
from scipy.stats import fisher_exact,mannwhitneyu
from statsmodels.stats.multitest import multipletests
import umap,hdbscan
from pyod.models.knn import KNN
warnings.filterwarnings('ignore', category=UserWarning)
al=pd.read_csv('results/aligned_predictions_all_models.csv');ext=pd.read_csv('results/external_features.csv')
df=al.merge(ext,on='bnumber',how='left').fillna(0)
v=pd.read_csv('results/ensemble_v2_oof.csv')[['bnumber','v2_lr']]
df=df.merge(v,on='bnumber',validate='one_to_one')
base=['fba_min','fba_rich','gnn','cnn','kmer'];NEW3=['pax_log_ppm','cai','gc3']
cols=base+[c for c in ext.columns if c!='bnumber' and not c.startswith('strctx') and c not in NEW3]
X=StandardScaler().fit_transform(df[cols].to_numpy(float));n=len(X);assert n==1249
error=set(df[df.essential==0].nlargest(100,'v2_lr').bnumber)
yerr=df.bnumber.isin(error).to_numpy()
emb=umap.UMAP(n_neighbors=30,min_dist=0.1,random_state=17,n_jobs=1).fit_transform(X)
labels=hdbscan.HDBSCAN(min_cluster_size=25,min_samples=10).fit_predict(emb)
outlier=KNN(n_neighbors=20,contamination=0.1);outlier.fit(X)
score=outlier.decision_scores_.astype(float)
rows=[]
for k in sorted(np.unique(labels)):
    mask=labels==k;a=int((mask & yerr).sum());b=int((mask & ~yerr).sum());c=int((~mask & yerr).sum());d=int((~mask & ~yerr).sum())
    odds,p=fisher_exact([[a,b],[c,d]],alternative='greater')
    rows.append({'cluster':int(k),'size':int(mask.sum()),'top100_false_priorities':a,'fraction':a/max(int(mask.sum()),1),
                 'odds_ratio':float(odds),'p_one_sided':float(p)})
q=multipletests([x['p_one_sided'] for x in rows],method='fdr_bh')[1]
for x,vv in zip(rows,q):x['fdr_bh']=float(vv)
test=mannwhitneyu(score[yerr],score[~yerr],alternative='two-sided')
res={'design':'exploratory unsupervised feature-space test; labels inspected after embedding/clustering',
     'n_genes':n,'n_features':len(cols),'n_clusters_excluding_noise':sum(x['cluster']!=-1 for x in rows),
     'clusters':rows,'pyod_outlier_score':{'median_false_priority':float(np.median(score[yerr])),
                       'median_other':float(np.median(score[~yerr])),'two_sided_mwu_p':float(test.pvalue)},
     'packages_executed':['umap-learn','hdbscan','PyOD'],
     'input_sha256':{p:hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in ['results/aligned_predictions_all_models.csv','results/external_features.csv','results/ensemble_v2_oof.csv']},
     'caveat':'Embedding distances and HDBSCAN clusters are algorithmic summaries, not cell types; same genes/labels as prior model; no independent validation.'}
Path('results/error_structure.json').write_text(json.dumps(res,indent=2)+'\n')
pd.DataFrame({'bnumber':df.bnumber,'essential':df.essential,'v2_lr':df.v2_lr,
              'false_priority_top100':yerr,'umap_x':emb[:,0],'umap_y':emb[:,1],
              'cluster':labels,'pyod_score':score}).to_csv('results/error_structure_gene_scores.csv',index=False)
print(json.dumps(res,indent=2))
