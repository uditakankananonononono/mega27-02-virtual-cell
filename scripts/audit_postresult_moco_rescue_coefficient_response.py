"""Post-result paired MoCo objective coefficient response in one source model."""
import hashlib,json,math
from pathlib import Path
import cobra
p=Path('data/bigg_survey/iJN1463.json');expected='8232a1ea39af9b9b0313601a3cefca1028e50e53b2ca2c06002c1a9adf7568bf'
assert hashlib.sha256(p.read_bytes()).hexdigest()==expected
m=cobra.io.load_json_model(str(p));b=m.reactions.get_by_id('BIOMASS_KT2440_WT3');met=['bmocogdp_c','mocogdp_c'];coef={z:float(b.get_coefficient(z)) for z in met};assert all(abs(c+0.000223)<1e-11 for c in coef.values()) and abs(b.get_coefficient('mobd_c')+0.002817)<1e-11
prior=json.load(open('results/postresult_moco_uptake_tolerance.json'));genes=prior['genes'];assert len(genes)==6 and m.reactions.get_by_id('EX_glc__D_e').lower_bound==-6
native=prior['uptake_variants']['1.0'];native_wt=native['unedited_wt']['growth'];assert native_wt>0
out={'status':'post-result same-model LP coefficient grid, no phenotype validation','protocol':'notes/postresult_moco_rescue_coefficient_response_protocol.md','source_url':'https://bigg.ucsd.edu/static/models/iJN1463.json','source_sha256':expected,'native_coefficients':coef,'native_wt_growth':native_wt,'genes':genes,'grid':[],'limits':'Changing paired biomass coefficients in one model is not a feasible biological intervention; fractions chosen after original rescue seen. No continuous threshold or generality; original G2 FAIL.'}
for f in (0,.1,.25,.5,.75,.9,.99,1):
 with m:
  for z,c in coef.items():b.add_metabolites({m.metabolites.get_by_id(z):c*(f-1)},combine=True)
  assert all(abs((b.get_coefficient(z) if m.metabolites.get_by_id(z) in b.metabolites else 0)-coef[z]*f)<1e-10 for z in met)
  sol=m.optimize();assert sol.status=='optimal';wt=float(sol.objective_value)
  gene={}
  for g in genes:
   with m:
    m.genes.get_by_id(g).knock_out();s=m.optimize();assert s.status=='optimal';v=float(s.objective_value);assert math.isfinite(v)
    gene[g]={'status':s.status,'growth':v,'ratio_to_native_wt':v/native_wt,'rescued_at_95pct':v/native_wt>=.95}
  out['grid'].append({'remaining_fraction':f,'wt_growth':wt,'wt_ratio_to_native_wt':wt/native_wt,'knockouts':gene,'n_rescued_95pct':sum(v['rescued_at_95pct'] for v in gene.values())})
assert abs(out['grid'][-1]['wt_growth']-native_wt)<1e-7
for g in genes:
 assert abs(out['grid'][0]['knockouts'][g]['growth']-native['patched_gene_knockouts'][g]['growth'])<1e-7
 assert abs(out['grid'][-1]['knockouts'][g]['growth']-native['unedited_gene_knockouts'][g]['growth'])<1e-7
out['first_tested_positive_fraction_losing_rescue']={g:next((v['remaining_fraction'] for v in out['grid'] if v['remaining_fraction']>0 and not v['knockouts'][g]['rescued_at_95pct']),None) for g in genes}
Path('results/postresult_moco_rescue_coefficient_response.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print('native',native_wt,'grid',[(a['remaining_fraction'],a['wt_growth'],a['n_rescued_95pct'],sorted(set(round(z['growth'],9) for z in a['knockouts'].values()))) for a in out['grid']]);print(out['first_tested_positive_fraction_losing_rescue'])
