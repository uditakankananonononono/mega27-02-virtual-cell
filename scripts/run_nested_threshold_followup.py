"""Exploratory nested threshold evaluation; frozen protocol in notes/prereg_nested_threshold_followup.md."""
import hashlib,json
from pathlib import Path
import numpy as np,pandas as pd
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.metrics import roc_auc_score,average_precision_score,f1_score,matthews_corrcoef
from scipy.stats import binomtest
A=Path('results/aligned_predictions_all_models.csv');E=Path('results/external_features.csv')
hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (A,E)}
a=pd.read_csv(A);e=pd.read_csv(E)
assert a.bnumber.is_unique and e.bnumber.is_unique
assert len(a)==1249
D=a.merge(e,on='bnumber',how='left',validate='one_to_one').fillna(0)
y=D.essential.to_numpy(dtype=int)
NEW3={'pax_log_ppm','cai','gc3'}
clean=[c for c in e if c!='bnumber' and not c.startswith('string_') and c!='up_annot' and c not in NEW3]
cols=['fba_min','fba_rich']+clean
X=D[cols].to_numpy(dtype=float)
assert np.isfinite(X).all() and len(y)==len(X)

def model():return make_pipeline(StandardScaler(),LogisticRegression(class_weight='balanced',C=0.5,max_iter=4000,random_state=7))
def threshold(score,truth):
    grid=np.quantile(score,np.linspace(.005,.995,199))
    best=None
    for t in grid:
        m=matthews_corrcoef(truth,score>=t)
        key=(-m,abs(t-.5),t)
        if best is None or key<best[0]:best=(key,t)
    return float(best[1])

scores=np.full(len(y),np.nan);hard=np.full(len(y),-1);folds=[]
for k,(tr,te) in enumerate(StratifiedKFold(3,shuffle=True,random_state=7).split(X,y)):
    inner=np.full(len(tr),np.nan)
    for itr,ite in StratifiedKFold(3,shuffle=True,random_state=17).split(X[tr],y[tr]):
        inner[ite]=model().fit(X[tr][itr],y[tr][itr]).predict_proba(X[tr][ite])[:,1]
    assert np.isfinite(inner).all()
    t=threshold(inner,y[tr])
    scores[te]=model().fit(X[tr],y[tr]).predict_proba(X[te])[:,1]
    hard[te]=(scores[te]>=t).astype(int)
    folds.append({'fold':k,'n_outer_train':len(tr),'n_outer_eval':len(te),'inner_oof_n':len(inner),'threshold':t,
                  'inner_mcc':float(matthews_corrcoef(y[tr],inner>=t))})
assert np.isfinite(scores).all() and (hard>=0).all()
base=(D.fba_min.to_numpy()>.5).astype(int)
ok,baseok=hard==y,base==y
b=int((ok & ~baseok).sum());c=int((~ok & baseok).sum())
out={'status':'new exploratory same-accession nested follow-up; original G2 FAIL retained; no benchmark gate repair',
     'protocol':'notes/prereg_nested_threshold_followup.md','input_sha256':hashes,'n_genes':len(y),'n_features':len(cols),'feature_names':cols,'folds':folds,
     'model':{'auroc':roc_auc_score(y,scores),'auprc':average_precision_score(y,scores),'f1':f1_score(y,hard),'mcc':matthews_corrcoef(y,hard)},
     'fba_fixed':{'auroc':roc_auc_score(y,D.fba_min),'f1':f1_score(y,base),'mcc':matthews_corrcoef(y,base)},
     'discordance':{'b_model_only_correct':b,'c_fba_only_correct':c,'p_exact_two_sided':float(binomtest(b,b+c).pvalue)},
     'original_gate_verdict':json.load(open('results/calibration_gate.json'))['verdict']}
Path('results/nested_threshold_followup.json').write_text(json.dumps(out,indent=2)+'\n')
pd.DataFrame({'bnumber':D.bnumber,'essential':y,'score':scores,'hard':hard,'fba_fixed':base}).to_csv('results/nested_threshold_followup_oof.csv',index=False)
print(json.dumps({k:v for k,v in out.items() if k not in ('feature_names','input_sha256')},indent=2))
