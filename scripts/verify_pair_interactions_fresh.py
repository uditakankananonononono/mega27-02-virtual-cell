"""Independent model reload for three interaction surprises, with LP statuses and WT."""
import csv,hashlib,json
from pathlib import Path
import cobra
p=Path('data/bigg_survey/iJN1463.json')
expected=next(z['sha256'] for z in csv.DictReader(open('results/bigg_model_survey.csv')) if z['accession']=='iJN1463')
assert hashlib.sha256(p.read_bytes()).hexdigest()==expected
prior=json.load(open('results/postresult_pair_interactions.json'))
singles=json.load(open('results/biomass_single_term_sensitivity.json'))
source={z['metabolite_id']:z for z in singles['rows']}
pairs=json.load(open('results/postresult_equal_coefficient_pair_null.json'))
pp={tuple(z['removed']):z for z in pairs['pairs']}
checks=[('2fe2s_c','btamp_c'),('hemeO_c','pheme_c'),('adocbl_c','fad_c')]
rows=[]
for terms in checks:
 pair=pp[terms]; gset=set(pair['rescued_ko_growth'])|set(source[terms[0]]['rescued_genes'])|set(source[terms[1]]['rescued_genes'])
 values=[]
 for removed in [(),(terms[0],),(terms[1],),terms]:
  m=cobra.io.load_json_model(str(p)); obj=m.reactions.get_by_id('BIOMASS_KT2440_WT3'); coeff={z.id:float(c) for z,c in obj.metabolites.items()}
  obj.subtract_metabolites({m.metabolites.get_by_id(z):coeff[z] for z in removed})
  wt=m.optimize(); assert wt.status=='optimal'
  outcomes={}
  for gene in sorted(gset):
   with m:
    m.genes.get_by_id(gene).knock_out();res=m.optimize();outcomes[gene]={'status':res.status,'growth':float(res.objective_value) if res.objective_value is not None else None}
  values.append({'removed':list(removed),'wt_growth':float(wt.objective_value),'gene_knockouts':outcomes})
 assert abs(values[-1]['wt_growth']-pair['edited_wt_growth'])<1e-8
 restored={g for g,r in values[-1]['gene_knockouts'].items() if r['growth'] is not None and r['growth']>=.95*.5861175448479794}
 assert restored==set(pair['rescued_ko_growth'])
 rows.append({'pair':list(terms),'candidate_genes':sorted(gset),'variants':values})
out={'status':'fresh independent model reload of three outcome-known pair interactions; not new gate',
 'source_url':'https://bigg.ucsd.edu/static/models/iJN1463.json','source_sha256':expected,
 'rescue_predicate':'KO growth >= .95 * original WT growth .5861175448479794',
 'pairs':rows,'limits':'Objective coefficients are changed, which may change feasible growth non-monotonically; no biochemical epistasis or biological viability inferred.'}
Path('results/postresult_pair_interactions_fresh_check.json').write_text(json.dumps(out,indent=2)+'\n')
for row in rows:
 print(row['pair'],'wild types',[round(z['wt_growth'],6) for z in row['variants']],
       'pair rescued',len([r for r in row['variants'][-1]['gene_knockouts'].values() if r['growth'] is not None and r['growth']>=.95*.5861175448479794]))
