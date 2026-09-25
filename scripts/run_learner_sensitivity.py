"""Exploratory robustness of v2 features across four genuinely used science packages."""
import hashlib, json, warnings
import numpy as np, pandas as pd
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from imblearn.over_sampling import RandomOverSampler
import lightgbm as lgb
from catboost import CatBoostClassifier
import shap
warnings.filterwarnings('ignore', category=UserWarning)
al=pd.read_csv('results/aligned_predictions_all_models.csv')
ext=pd.read_csv('results/external_features.csv')
df=al.merge(ext,on='bnumber',how='left').fillna(0);y=df.essential.to_numpy()
base=['fba_min','fba_rich','gnn','cnn','kmer']; NEW3=['pax_log_ppm','cai','gc3']
extc=[c for c in ext.columns if c!='bnumber' and not c.startswith('strctx') and c not in NEW3]
cols=base+extc;X=df[cols].to_numpy(dtype=float)
prior=pd.read_csv('results/ensemble_v2_oof.csv').set_index('bnumber').loc[df.bnumber,'v2_lr'].to_numpy()
models={'lightgbm':np.zeros(len(y)), 'catboost':np.zeros(len(y)), 'imblearn_oversampled_lr':np.zeros(len(y))}
for fold,(tr,te) in enumerate(StratifiedKFold(3,shuffle=True,random_state=7).split(X,y)):
    a=lgb.LGBMClassifier(n_estimators=120,max_depth=3,num_leaves=7,learning_rate=0.03,
                         min_child_samples=20,verbose=-1,n_jobs=1,random_state=7)
    a.fit(X[tr],y[tr]);models['lightgbm'][te]=a.predict_proba(X[te])[:,1]
    b=CatBoostClassifier(iterations=120,depth=3,learning_rate=0.03,loss_function='Logloss',
                         verbose=False,thread_count=1,random_seed=7,allow_writing_files=False)
    b.fit(X[tr],y[tr]);models['catboost'][te]=b.predict_proba(X[te])[:,1]
    # Oversample ONLY the training fold; holdout observations never enter resampling.
    xtr,ytr=RandomOverSampler(random_state=7).fit_resample(X[tr],y[tr])
    c=make_pipeline(StandardScaler(),LogisticRegression(max_iter=4000,C=0.5))
    c.fit(xtr,ytr);models['imblearn_oversampled_lr'][te]=c.predict_proba(X[te])[:,1]
metrics={'v2_lr':{'auroc':float(roc_auc_score(y,prior)), 'auprc':float(average_precision_score(y,prior))}}
rng=np.random.default_rng(11)
ids=[rng.integers(0,len(y),len(y)) for _ in range(1000)]
ids=[a for a in ids if 0<y[a].sum()<len(a)]
for name,scores in models.items():
    dif=np.array([roc_auc_score(y[i],scores[i])-roc_auc_score(y[i],prior[i]) for i in ids])
    metrics[name]={'auroc':float(roc_auc_score(y,scores)),'auprc':float(average_precision_score(y,scores)),
                   'delta_auroc_vs_v2_lr':float(roc_auc_score(y,scores)-roc_auc_score(y,prior)),
                   'paired_ci95':list(map(float,np.quantile(dif,[0.025,0.975])))}
full=CatBoostClassifier(iterations=120,depth=3,learning_rate=0.03,verbose=False,thread_count=1,random_seed=7,allow_writing_files=False)
full.fit(X,y)
explainer=shap.TreeExplainer(full)
sv=explainer.shap_values(X[:64]);sv=np.array(sv)
if sv.ndim==3:sv=sv[:,:,1]
imp=pd.Series(np.mean(np.abs(sv),axis=0),index=cols).sort_values(ascending=False)
res={'design':'exploratory, frozen folds and hyperparameters, no model selection',
     'n_genes':len(y),'n_features':len(cols),'science_packages_executed':['lightgbm','catboost','imbalanced-learn','shap'],
     'metrics':metrics,'shap_catboost_top10':{k:float(v) for k,v in imp.head(10).items()},
     'input_sha256':{p:hashlib.sha256(open(p,'rb').read()).hexdigest() for p in ['results/aligned_predictions_all_models.csv','results/external_features.csv','results/ensemble_v2_oof.csv']},
     'caveat':'SHAP values describe in-sample full-data model only; comparative AUROC/AUPRC use OOF predictions. Not an independent biological validation.'}
with open('results/learner_sensitivity.json','w') as f:json.dump(res,f,indent=2)
pd.DataFrame({'bnumber':df.bnumber,'essential':y,**models}).to_csv('results/learner_sensitivity_oof.csv',index=False)
print(json.dumps(res,indent=2))
