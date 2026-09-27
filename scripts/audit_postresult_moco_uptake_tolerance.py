"""Post-result within-model glucose bound sensitivity, not phenotype validation."""
import hashlib,json
from pathlib import Path
import cobra
src=Path('data/bigg_survey/iJN1463.json');expected='8232a1ea39af9b9b0313601a3cefca1028e50e53b2ca2c06002c1a9adf7568bf'
assert hashlib.sha256(src.read_bytes()).hexdigest()==expected
m=cobra.io.load_json_model(str(src));g=sorted(json.load(open('results/postresult_published_model_counterfactual.json'))['models']['BiGG_iJN1463']['panels']['MoCo']['edited_genes'])
assert len(g)==6 and m.reactions.get_by_id('EX_glc__D_e').lower_bound==-6
bm=m.reactions.get_by_id('BIOMASS_KT2440_WT3')
for id in ['bmocogdp_c','mocogdp_c']:
 assert abs(bm.get_coefficient(id)+.000223)<1e-10
out={'status':'post-result same-model medium-bound sensitivity; not phenotype or new gate','protocol':'notes/postresult_moco_uptake_tolerance_protocol.md','source_url':'https://bigg.ucsd.edu/static/models/iJN1463.json','source_sha256':expected,'genes':g,'uptake_variants':{},'limits':'Retains all other native model bounds, oxygen and nutrients. Increasing uptake is not a measured growth medium. Same model and preselected genes, not independent validation. G2 remains FAIL.'}
def opt(model):
 sol=model.optimize();return {'status':sol.status,'growth':float(sol.objective_value) if sol.status=='optimal' else None}
for mult in [.5,1.,1.5,2.]:
 with m:
  bound=-6*mult;m.reactions.get_by_id('EX_glc__D_e').lower_bound=bound
  wt=opt(m);assert wt['status']=='optimal'
  baseline={}
  for id in g:
   with m:
    m.genes.get_by_id(id).knock_out();baseline[id]=opt(m)
  bm.subtract_metabolites({m.metabolites.get_by_id(id):-0.000223 for id in ['bmocogdp_c','mocogdp_c']})
  edited=opt(m);assert edited['status']=='optimal'
  patched={}
  for id in g:
   with m:
    m.genes.get_by_id(id).knock_out();patched[id]=opt(m)
  denom=wt['growth'];rescued=[id for id in g if denom>0 and patched[id]['growth'] is not None and patched[id]['growth']/denom>=.95]
  out['uptake_variants'][str(mult)]={'glucose_exchange_lower_bound':bound,'unedited_wt':wt,'unedited_gene_knockouts':baseline,'patched_wt':edited,'patched_gene_knockouts':patched,'patched_over_unedited_wt':{id:patched[id]['growth']/denom if denom>0 and patched[id]['growth'] is not None else None for id in g},'rescued_95pct':rescued,'ratio_defined':denom>0}
Path('results/postresult_moco_uptake_tolerance.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
for k,z in out['uptake_variants'].items():print(k,z['glucose_exchange_lower_bound'],z['unedited_wt']['growth'],z['patched_wt']['growth'],len(z['rescued_95pct']),z['patched_over_unedited_wt'])
