"""Fixed target-panel reconstruction-version transfer; no phenotype evidence."""
import csv,hashlib,json
from pathlib import Path
from lxml import etree
import cobra
src=Path('/tmp/nogales-supp/EMI-22-255-s003.xml');bigg=Path('data/bigg_survey/iJN1463.json')
assert hashlib.sha256(src.read_bytes()).hexdigest()=='0c8113ee6412c21df1afbd0682d1e89e73a515b558a2cab29ee8213bb5afb5f0'
expected=next(z['sha256'] for z in csv.DictReader(open('results/bigg_model_survey.csv')) if z['accession']=='iJN1463')
assert hashlib.sha256(bigg.read_bytes()).hexdigest()==expected
xml=etree.parse(str(src));m=xml.getroot()[0];cs=m.find('{*}listOfCompartments');blank=[z for z in cs if z.get('id')==''];assert len(blank)==1 and not any(z.get('compartment')=='' for z in m.find('{*}listOfSpecies'));cs.remove(blank[0]);tmp=Path('/tmp/nogales-supp/s3-fixed-temp.xml');xml.write(str(tmp),encoding='utf-8',xml_declaration=True)
mo=json.load(open('results/moco_subset_interaction_ijn1463.json'));mo_ids=sorted(next(z for z in mo['subsets'] if z['removed']==['bmocogdp_c','mocogdp_c'])['rescued_ko_growth'])
ex=json.load(open('results/postresult_pair_exception_reactions.json'));pair_ids=sorted(next(z for z in ex['exceptions'] if z['pair']==['2fe2s_c','btamp_c'])['lost_union_genes'])
assert len(mo_ids)==6 and len(pair_ids)==8
panels={'MoCo':{'genes':mo_ids,'remove':['bmocogdp_c','mocogdp_c']},'suppressive_pair':{'genes':pair_ids,'remove':['2fe2s_c','btamp_c']}}
out={'status':'post-result reconstruction version transfer; no independent phenotype','protocol':'notes/postresult_published_model_counterfactual_protocol.md','original_published_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'bigg_sha256':expected,'panels':panels,'models':{}}
for label,path in [('published_iJN1462',tmp),('BiGG_iJN1463',bigg)]:
 model=cobra.io.read_sbml_model(str(path)) if label.startswith('published') else cobra.io.load_json_model(str(path))
 bio=next(r for r in model.reactions if r.objective_coefficient)
 def solve():
  x=model.optimize();return {'status':x.status,'growth':float(x.objective_value) if x.objective_value is not None else None}
 baseline=solve();assert baseline['status']=='optimal' and baseline['growth']>.1
 result={'model_id':model.id,'objective':bio.id,'baseline':baseline,'panels':{}}
 for panel,p in panels.items():
  raw={g:None for g in p['genes']}
  for g in p['genes']:
   with model:
    model.genes.get_by_id(g).knock_out();raw[g]=solve()
  coeff={t:float(bio.metabolites[model.metabolites.get_by_id(t)]) for t in p['remove']}
  with model:
   bio.subtract_metabolites({model.metabolites.get_by_id(t):coeff[t] for t in p['remove']})
   edited_wt=solve();edited={}
   for g in p['genes']:
    with model:
     model.genes.get_by_id(g).knock_out();edited[g]=solve()
  assert all(z['status']=='optimal' for z in [edited_wt,*raw.values(),*edited.values()])
  result['panels'][panel]={'removed_coefficients':coeff,'unedited_genes':raw,'edited_wt':edited_wt,'edited_genes':edited,
   'rescued_over_95pct_original_wt':[g for g,x in edited.items() if x['growth']>=.95*baseline['growth']]}
 out['models'][label]=result
Path('results/postresult_published_model_counterfactual.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
for label,x in out['models'].items():
 print(label,'WT',x['baseline']['growth'])
 for p,z in x['panels'].items():print(p,'edited WT',z['edited_wt']['growth'],'rescued',z['rescued_over_95pct_original_wt'],'KO values',sorted({round(q['growth'],8) for q in z['edited_genes'].values()}))
