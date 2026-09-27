"""Fixed post-result suppressive-pair flux diagnostic; see protocol note."""
import csv,hashlib,json
from pathlib import Path
import cobra
import numpy as np
from cobra.util.array import create_stoichiometric_matrix
p=Path('data/bigg_survey/iJN1463.json')
expected=next(z['sha256'] for z in csv.DictReader(open('results/bigg_model_survey.csv')) if z['accession']=='iJN1463')
assert hashlib.sha256(p.read_bytes()).hexdigest()==expected
prior=json.load(open('results/postresult_pair_interactions_fresh_check.json'))
row=next(z for z in prior['pairs'] if z['pair']==['2fe2s_c','btamp_c'])
source=json.load(open('results/postresult_pair_exception_reactions.json'))
ex=next(z for z in source['exceptions'] if z['pair']==row['pair'])
genes=ex['lost_union_genes'];assert genes==row['candidate_genes'] and len(genes)==8
out={'source_sha256':expected,'source_url':'https://bigg.ucsd.edu/static/models/iJN1463.json',
     'protocol':'notes/postresult_suppressive_pair_flux_audit.md','status':'post-result model-internal diagnosis; no gate update',
     'pair':row['pair'],'candidate_genes':genes,'variants':[]}
for saved in row['variants']:
 removed=saved['removed'];m=cobra.io.load_json_model(str(p));bio=m.reactions.get_by_id('BIOMASS_KT2440_WT3')
 native={t:float(bio.metabolites[m.metabolites.get_by_id(t)]) for t in row['pair']}
 bio.subtract_metabolites({m.metabolites.get_by_id(t):native[t] for t in removed})
 def solve():
  sol=m.optimize();assert sol.status=='optimal',sol.status
  assert abs(sol.fluxes[bio.id]-sol.objective_value)<1e-8
  stoich=create_stoichiometric_matrix(m,array_type='lil')
  flux=sol.fluxes.loc[[r.id for r in m.reactions]].to_numpy()
  residual=float(np.max(np.abs(stoich.tocsr()@flux)))
  assert residual<1e-6,residual
  # The cobra solution carries primal metabolite shadow prices, not a validated dual certificate.
  exchanges={r.id:float(sol.fluxes[r.id]) for r in m.exchanges if abs(sol.fluxes[r.id])>1e-7}
  return {'status':sol.status,'growth':float(sol.objective_value),'biomass_flux':float(sol.fluxes[bio.id]),'max_mass_balance_residual':residual,'nonzero_exchanges':exchanges}
 wt=solve();assert abs(wt['growth']-saved['wt_growth'])<1e-8
 ko={}
 for gene in genes:
  with m:
   m.genes.get_by_id(gene).knock_out();ko[gene]=solve()
  assert abs(ko[gene]['growth']-saved['gene_knockouts'][gene]['growth'])<1e-7
  ko[gene]['growth_over_variant_wt']=ko[gene]['growth']/wt['growth']
  ko[gene]['growth_over_original_wt']=ko[gene]['growth']/row['variants'][0]['wt_growth']
 out['variants'].append({'removed':removed,'original_coefficients':native,'wild_type':wt,'knockouts':ko})
out['specific_reactions']={rid:{'reaction':m.reactions.get_by_id(rid).reaction,'bounds':m.reactions.get_by_id(rid).bounds,'gpr':m.reactions.get_by_id(rid).gene_reaction_rule} for rid in ['BTS5','BACCL','I2FE2ST','LIPOS']}
out['conclusion']='The double edit has lower WT optimum and all eight previously rescued KO objectives return to zero, not just below a relative rescue threshold. These are model LP outcomes; not unique flux mechanisms or cell phenotypes.'
Path('results/postresult_suppressive_pair_flux.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
for v in out['variants']:
 print(v['removed'],'WT',v['wild_type']['growth'],'KO',sorted(set(round(z['growth'],9) for z in v['knockouts'].values())), 'active exchanges',len(v['wild_type']['nonzero_exchanges']))
