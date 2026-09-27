"""Artificial intracellular precursor supply in native iJN1463 only."""
import hashlib,json,math
from pathlib import Path
import cobra
from cobra.util.array import create_stoichiometric_matrix
import numpy as np
p=Path('data/bigg_survey/iJN1463.json');sha='8232a1ea39af9b9b0313601a3cefca1028e50e53b2ca2c06002c1a9adf7568bf';assert hashlib.sha256(p.read_bytes()).hexdigest()==sha
m=cobra.io.load_json_model(str(p));genes=json.load(open('results/postresult_moco_uptake_tolerance.json'))['genes'];assert genes==['PP_0735','PP_1292','PP_1293','PP_1294','PP_2123','PP_3457']
coef={z:m.reactions.get_by_id('BIOMASS_KT2440_WT3').get_coefficient(z) for z in ('bmocogdp_c','mocogdp_c')};assert set(coef.values())=={-.000223}
variants=[('native',[]),('virtual_moco_c',['moco_c']),('virtual_bmoco1gdp_c',['bmoco1gdp_c']),('virtual_both_terminal_GDP',['mocogdp_c','bmocogdp_c'])]
source={}
for z in {x for _,lst in variants for x in lst}:
 met=m.metabolites.get_by_id(z);source[z]={'reactions':[{'id':r.id,'stoichiometric_coefficient':r.metabolites[met],'gene_rule':r.gene_reaction_rule,'bounds':list(r.bounds)} for r in sorted(met.reactions,key=lambda x:x.id)],'existing_supply_reactions':[r.id for r in met.reactions if r.id.startswith(('DM_','SK_','EX_'))]};assert not source[z]['existing_supply_reactions']
wt=m.optimize();assert wt.status=='optimal' and wt.objective_value>0;denom=float(wt.objective_value)
out={'status':'post-result artificial-supply LP diagnostic, not biological phenotype','protocol':'notes/postresult_moco_substrate_bypass_protocol.md','source_url':'https://bigg.ucsd.edu/static/models/iJN1463.json','source_sha256':sha,'genes':genes,'native_biomass_coefficients':coef,'native_wt_growth':denom,'pathway_source_reactions':source,'variants':[],'limits':'Positive-coefficient cytosolic supply is an artificial internal model reaction not a nutrient or demonstrated biological intervention. Only one reconstruction and preselected targets. Original G2 FAIL.'}
for name,lst in variants:
 with m:
  virtual=[]
  for z in lst:
   r=cobra.Reaction('virtual_supply_'+z);r.lower_bound=0;r.upper_bound=1000;r.add_metabolites({m.metabolites.get_by_id(z):1});m.add_reactions([r]);virtual.append(r)
  def solve():
   sol=m.optimize();assert sol.status=='optimal'
   flux=sol.fluxes.loc[[r.id for r in m.reactions]].to_numpy();res=float(np.max(np.abs(create_stoichiometric_matrix(m,array_type='lil').tocsr()@flux)));assert res<1e-6
   v=float(sol.objective_value);assert math.isfinite(v)
   return {'status':sol.status,'growth':v,'ratio_to_native_wt':v/denom,'rescued_95pct':v/denom>=.95,'virtual_supply_fluxes':{r.id:float(sol.fluxes[r.id]) for r in virtual},'mass_balance_max_abs':res}
  wild=solve();ko={}
  for g in genes:
   with m:m.genes.get_by_id(g).knock_out();ko[g]=solve()
  out['variants'].append({'name':name,'virtual_supplies':lst,'wild_type':wild,'knockouts':ko,'n_rescued_95pct':sum(z['rescued_95pct'] for z in ko.values())})
prior=json.load(open('results/postresult_moco_uptake_tolerance.json'))['uptake_variants']['1.0'];assert abs(out['variants'][0]['wild_type']['growth']-prior['unedited_wt']['growth'])<1e-7
for g in genes:assert abs(out['variants'][0]['knockouts'][g]['growth']-prior['unedited_gene_knockouts'][g]['growth'])<1e-7
Path('results/postresult_moco_substrate_bypass.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
for z in out['variants']:print(z['name'],'WT',z['wild_type']['growth'],'rescued',z['n_rescued_95pct'],{g:(round(v['growth'],7),v['virtual_supply_fluxes']) for g,v in z['knockouts'].items()})
