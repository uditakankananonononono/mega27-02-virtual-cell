"""Frozen iJN1463 -> iYO844 transfer, source-verified and checkpointed."""
import csv,hashlib,json,sys
from pathlib import Path
import cobra
import numpy as np
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score,average_precision_score,brier_score_loss

source={r['accession']:r for r in csv.DictReader(open('results/bigg_model_survey.csv'))}
p=Path(source['iYO844']['path']);sha=hashlib.sha256(p.read_bytes()).hexdigest();assert sha==source['iYO844']['sha256']
training=json.load(open('results/biomass_topology_manifest.json'))
old=json.load(open('results/biomass_single_term_sensitivity.json'))
assert training['sha256']==old['model_sha256']==source['iJN1463']['sha256']
lookup={z['metabolite_id']:z['n_baseline_essential_rescued'] for z in old['rows']}
assert len(lookup)==len(training['terms'])==102
keys=['abs_biomass_coefficient','non_biomass_reaction_count','gene_associated_reaction_count','active_producer_reaction_count_wt','active_producer_abs_flux_sum_wt']
def vec(z):return [np.log1p(float(z[k])) if k=='abs_biomass_coefficient' else float(z[k]) for k in keys]
X=np.array([vec(z) for z in training['terms']]);y=np.array([lookup[z['metabolite_id']]>0 for z in training['terms']],int)
assert int(y.sum())==23
model=make_pipeline(StandardScaler(),LogisticRegression(C=1,class_weight='balanced',solver='liblinear',random_state=0,max_iter=1000))
model.fit(X,y)
m=cobra.io.load_json_model(str(p));ob=next(z for z in m.reactions if z.objective_coefficient>0)
assert ob.id=='BIOMASS_BS_10';wt_sol=m.optimize();assert wt_sol.status=='optimal' and wt_sol.objective_value>0
wt=float(wt_sol.objective_value);terms=[]
for met,c in sorted(ob.metabolites.items(),key=lambda x:x[0].id):
 if c>=0:continue
 links=[r for r in met.reactions if r.id!=ob.id]
 producing=[r for r in links if r.metabolites[met]*float(wt_sol.fluxes[r.id])>1e-9]
 z={'metabolite_id':met.id,'name':met.name,'abs_biomass_coefficient':abs(float(c)),
    'non_biomass_reaction_count':len(links),'gene_associated_reaction_count':sum(bool(r.genes) for r in links),
    'active_producer_reaction_count_wt':len(producing),'active_producer_abs_flux_sum_wt':round(sum(abs(float(wt_sol.fluxes[r.id])) for r in producing),8)}
 z['predicted_probability']=float(model.predict_proba(np.array([vec(z)]))[0,1]);terms.append(z)
assert len(terms)==60
pred=Path('results/biomass_transfer_iyo844_frozen_predictions.json')
manifest={'source_url':source['iYO844']['source_url'],'source_sha256':sha,'training_sha256':training['sha256'],
          'protocol':'notes/prereg_biomass_transfer_iyo844.md','objective':ob.id,'wt_growth':wt,
          'feature_keys':keys,'training_n':len(y),'training_positive':int(y.sum()),'terms':terms}
if pred.exists():assert json.loads(pred.read_text())==manifest
else:pred.write_text(json.dumps(manifest,indent=2)+'\n')
if '--predictions-only' in sys.argv:print('Frozen predictions',len(terms),'in',pred);sys.exit(0)
base=cobra.flux_analysis.single_gene_deletion(m,processes=1)
essential=sorted(next(iter(ids)) for ids,gr in zip(base['ids'],base['growth']) if gr is not None and gr<1e-6)
assert essential
progress=Path('results/biomass_transfer_iyo844_progress.json')
done=json.loads(progress.read_text()) if progress.exists() else {'source_sha256':sha,'wt_growth':wt,'baseline_essential_count':len(essential),'rows':{}}
assert done['source_sha256']==sha and abs(done['wt_growth']-wt)<1e-10 and done['baseline_essential_count']==len(essential)
for z in terms:
 key=z['metabolite_id']
 if key in done['rows']:continue
 with m:
  o=m.reactions.get_by_id(ob.id);o.subtract_metabolites({m.metabolites.get_by_id(key):-z['abs_biomass_coefficient']})
  edited_wt=float(m.slim_optimize() or 0);rescued={}
  for gene in essential:
   with m:
    m.genes.get_by_id(gene).knock_out();gr=float(m.slim_optimize() or 0)
    if gr>=.95*wt:rescued[gene]=gr
 done['rows'][key]={'edited_wt_growth':edited_wt,'rescued_ko_growth':rescued}
 progress.write_text(json.dumps(done,indent=2)+'\n');print(len(done['rows']),'/',len(terms),key,len(rescued),flush=True)
ytest=np.array([bool(done['rows'][z['metabolite_id']]['rescued_ko_growth']) for z in terms],int)
proba=np.array([z['predicted_probability'] for z in terms]);rng=np.random.default_rng(20260926)
null_auc=[];null_ap=[]
for _ in range(1000):
 yp=rng.permutation(ytest)
 if len(set(yp))==2:null_auc.append(float(roc_auc_score(yp,proba)))
 null_ap.append(float(average_precision_score(yp,proba)))
auc=float(roc_auc_score(ytest,proba)) if len(set(ytest))==2 else None
ap=float(average_precision_score(ytest,proba));brier=float(brier_score_loss(ytest,proba))
res={'source_url':source['iYO844']['source_url'],'source_sha256':sha,'prediction_manifest':str(pred),
     'baseline_essential_count':len(essential),'n_terms':len(terms),'n_positive':int(ytest.sum()),
     'target_prevalence':float(ytest.mean()),'auroc':auc,'average_precision':ap,
     'brier':brier,'brier_constant_train_prevalence':float(brier_score_loss(ytest,np.full(len(ytest),float(y.mean())))),
     'null_feature_shuffle_auc_median':float(np.median(null_auc)) if null_auc else None,
     'null_feature_shuffle_ap_median':float(np.median(null_ap)),
     'null_feature_shuffle_auc_p_one_sided':float((1+sum(x>=auc for x in null_auc))/(1+len(null_auc))) if auc is not None else None,
     'null_feature_shuffle_ap_p_one_sided':float((1+sum(x>=ap for x in null_ap))/(1+len(null_ap))),
     'rows':[dict(z,**done['rows'][z['metabolite_id']]) for z in terms],
     'scope':'exploratory computational transfer across model families only; shared BiGG conventions; no independent biological validation or G2 repair'}
Path('results/biomass_transfer_iyo844.json').write_text(json.dumps(res,indent=2)+'\n')
print('COMPLETE', {k:res[k] for k in ['n_terms','n_positive','auroc','average_precision','target_prevalence','brier','brier_constant_train_prevalence','null_feature_shuffle_auc_p_one_sided','null_feature_shuffle_ap_p_one_sided']})
