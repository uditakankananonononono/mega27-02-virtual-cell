"""Post-result simple baseline and ATP/water influence audit, no new gate."""
import json
from pathlib import Path
import numpy as np
from sklearn.metrics import roc_auc_score,average_precision_score
out={'status':'post-result comparison of two fixed simple scores and ATP/water influence; no new held-out result',
     'protocol':'notes/postresult_biomass_baseline_audit.md','targets':{}}
for name in ['iyo844','imm904']:
 data=json.load(open('results/biomass_transfer_'+name+'.json'));rows=data['rows'];y=np.array([int(bool(r['rescued_ko_growth'])) for r in rows]);assert len(set(y))==2
 def metrics(ids,score):
  yy=y[ids];v=np.array(score)[ids];return {'n':int(len(yy)),'positives':int(yy.sum()),'prevalence':float(yy.mean()),
       'auroc':float(roc_auc_score(yy,v)),'average_precision':float(average_precision_score(yy,v))}
 n=np.arange(len(rows));freeze=[r['predicted_probability'] for r in rows]
 coef=[r['abs_biomass_coefficient'] for r in rows];degree=[-r['gene_associated_reaction_count'] for r in rows]
 trimmed=np.array([i for i,r in enumerate(rows) if r['metabolite_id'] not in {'atp_c','h2o_c'}])
 assert len(trimmed)==len(rows)-2
 full={k:metrics(n,v) for k,v in [('frozen_five_feature',freeze),('coefficient_magnitude',coef),('inverse_gene_associated_degree',degree)]}
 out['targets'][name]={'source_url':data['source_url'],'source_sha256':data['source_sha256'],'full':full,
    'frozen_minus_simple':{k:{metric:full['frozen_five_feature'][metric]-full[k][metric] for metric in ['auroc','average_precision']} for k in ['coefficient_magnitude','inverse_gene_associated_degree']},
    'without_atp_and_water':metrics(trimmed,freeze),
    'atp_and_water':[{key:r[key] for key in ['metabolite_id','predicted_probability','abs_biomass_coefficient','gene_associated_reaction_count']}|{'rescued_count':len(r['rescued_ko_growth'])} for r in rows if r['metabolite_id'] in {'atp_c','h2o_c'}],
    'per_term':[{key:r[key] for key in ['metabolite_id','predicted_probability','abs_biomass_coefficient','gene_associated_reaction_count']}|{'rescue_binary':int(bool(r['rescued_ko_growth']))} for r in rows]}
Path('results/biomass_baseline_dominance_audit.json').write_text(json.dumps(out,indent=2)+'\n')
for name,d in out['targets'].items():
 print(name,'full', {k:(round(v['auroc'],4),round(v['average_precision'],4)) for k,v in d['full'].items()},'trimmed',d['without_atp_and_water'])
