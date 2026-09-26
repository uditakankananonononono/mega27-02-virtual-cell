"""Read-only diagnostic of VC2 calibration gate score provenance and comparator calls.
No new preregistered gate, no overwrite of original calibration result.
"""
import hashlib,json,numpy as np,pandas as pd
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import matthews_corrcoef
al=pd.read_csv('results/aligned_predictions_all_models.csv')
oof=pd.read_csv('results/ensemble_v2_oof.csv')
assert al.bnumber.is_unique and oof.bnumber.is_unique and set(al.bnumber)==set(oof.bnumber)
d=al.merge(oof[['bnumber','essential','v2_lr']],on='bnumber',suffixes=('_al','_oof'),validate='one_to_one')
assert (d.essential_al==d.essential_oof).all()
y=d.essential_al.to_numpy().astype(int); folds=np.zeros(len(d),int)
for k,(_,te) in enumerate(StratifiedKFold(3,shuffle=True,random_state=7).split(np.zeros(len(y)),y)):folds[te]=k
r=json.load(open('results/calibration_gate.json'));scores=d.fba_min.to_numpy();fixed=(scores>.5)
rows=[];pred=np.zeros(len(y),bool)
for k in range(3):
 tr,te=folds!=k,folds==k;t=r['comparator']['thresholds'][f'fold{k}'];p=(scores[te]>=t);pred[te]=p
 rows.append({'fold':k,'n_train':int(tr.sum()),'n_eval':int(te.sum()),'train_label_influences_oof_scores':True,
   'threshold':t,'threshold_vs_fixed_calls_differ':int(np.sum(p!=fixed[te])),
   'near_threshold_score_count_eval_0.11_to_0.5':int(np.sum((scores[te]>=.110729)&(scores[te]<.5))),
   'near_threshold_score_count_eval_0.5_to_0.5903':int(np.sum((scores[te]>=.5)&(scores[te]<.590301))),
   'MCC_eval_fold_threshold':float(matthews_corrcoef(y[te],p)),
   'MCC_eval_fixed_0.5':float(matthews_corrcoef(y[te],fixed[te]))})
assert abs(matthews_corrcoef(y,pred)-r['comparator']['metrics']['mcc'])<1e-12
out={'status':'method-lineage audit; no new gate score or correction',
 'underlying_fit':'scripts/run_ensemble_v2.py oof() trains each held-out gene score from its complementary folds',
 'leak_path':'A later threshold for outer eval fold k uses the other two folds scores and labels; each of those OOF scores was trained on complement including fold k labels, so threshold selection is not completely independent of k.',
 'additional_leak_path':'The v2 LR stack includes CNN/GNN/kmer OOF base predictions trained on labels in different fold complements; final outer-fold independence requires all supervised base layers nested or a genuinely independent evaluation dataset.',
 'not_affected':'fba_min fixed >0.5 predictions are label-free model outputs, but the fitted comparator threshold uses labels; this audit does not invalidate independent LP biochemical mass-balance checks or MoCo in-silico rescue.',
 'first_result_retained':'results/calibration_gate.json gate G2 FAIL, no retroactive amendment; 93:67 R1 win is descriptive and not independently isolated.',
 'fold_rows':rows,'pooled_comparator_calls_differ_from_fixed':int(np.sum(pred!=fixed)),
 'sha256':{p:hashlib.sha256(open(p,'rb').read()).hexdigest() for p in ['results/aligned_predictions_all_models.csv','results/ensemble_v2_oof.csv','results/calibration_gate.json']}}
open('results/calibration_lineage_audit.json','w').write(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
