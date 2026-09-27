"""Compare saved single and pair gene sets; reread original saved outcome source."""
import json
from pathlib import Path
one=json.loads(Path('results/biomass_single_term_sensitivity.json').read_text())
two=json.loads(Path('results/postresult_equal_coefficient_pair_null.json').read_text())
assert one['model_sha256']==two['source_sha256']
assert len(one['rows'])==102 and len(two['pairs'])==253
singles={r['metabolite_id']:set(r['rescued_genes']) for r in one['rows']}
rows=[]
for pair in two['pairs']:
 a,b=pair['removed']; u=singles[a]|singles[b]; target=set(pair['rescued_ko_growth'])
 rows.append({'pair':[a,b],'single_a_count':len(singles[a]),'single_b_count':len(singles[b]),
              'union_count':len(u),'pair_count':len(target),
              'pair_only_genes':sorted(target-u),'single_union_only_genes':sorted(u-target),
              'shared_genes':sorted(u&target)})
mo=json.loads(Path('results/moco_subset_interaction_ijn1463.json').read_text())
mo_pair=next(r for r in mo['subsets'] if r['removed']==['bmocogdp_c','mocogdp_c'])
assert singles['bmocogdp_c']==singles['mocogdp_c']==set()
assert len(mo_pair['rescued_ko_growth'])==6
out={'status':'post-result saved single/pair set reconciliation; no prospective interaction test',
     'sources':['results/biomass_single_term_sensitivity.json','results/postresult_equal_coefficient_pair_null.json','results/moco_subset_interaction_ijn1463.json'],
     'model_sha256':one['model_sha256'],'n_pairs':len(rows),
     'n_pair_only':sum(bool(r['pair_only_genes']) for r in rows),
     'n_lost_single_union':sum(bool(r['single_union_only_genes']) for r in rows),
     'n_exact_single_union':sum(not r['pair_only_genes'] and not r['single_union_only_genes'] for r in rows),
     'mo_pair_saved':{'removed':mo_pair['removed'],'single_union_count':0,'pair_count':6,
                      'pair_only_genes':sorted(mo_pair['rescued_ko_growth'])},
     'rows':rows,
     'limits':'All outcomes and top patterns seen before protocol; objective-stoichiometry intervention, not measured biochemical epistasis or biological rescue; inspect unexpected loss with fresh LPs.'}
assert out['n_pair_only']==2 and out['n_lost_single_union']==1
Path('results/postresult_pair_interactions.json').write_text(json.dumps(out,indent=2)+'\n')
print({k:v for k,v in out.items() if k not in ('rows','sources')})
for r in rows:
 if r['pair_only_genes'] or r['single_union_only_genes']:print(r)
