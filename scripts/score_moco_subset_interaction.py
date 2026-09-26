"""Score locked exploratory MoCo coefficient subsets in BiGG iJN1463."""
import csv, hashlib, itertools, json
from pathlib import Path
import cobra
p=Path('data/bigg_survey/iJN1463.json')
row=next(z for z in csv.DictReader(open('results/bigg_model_survey.csv')) if z['accession']=='iJN1463')
assert hashlib.sha256(p.read_bytes()).hexdigest()==row['sha256']
previous=json.load(open('results/bigg_second_family_verification.json'))
assert previous['snapshot_sha256']==row['sha256']
m=cobra.io.load_json_model(str(p)); obj=m.reactions.get_by_id('BIOMASS_KT2440_WT3')
terms=tuple(sorted(previous['objective_negative_moco_coefficients']))
coef={z.id:float(c) for z,c in obj.metabolites.items() if z.id in terms}
assert set(terms)==set(coef) and all(coef[z]<0 for z in terms)
wt=float(m.slim_optimize()); assert abs(wt-previous['baseline_wt_growth'])<1e-8
base=cobra.flux_analysis.single_gene_deletion(m,processes=1)
essential=sorted(next(iter(ids)) for ids,gr in zip(base['ids'],base['growth']) if gr is not None and gr<1e-6)
assert len(essential)==previous['baseline_essential_genes']==262
out={'source_url':row['source_url'],'source_sha256':row['sha256'],'plan':'notes/prereg_moco_subset_interaction.md',
     'objective':obj.id,'baseline_growth':wt,'baseline_essential_count':len(essential),
     'rescue_threshold':'edited KO growth >= 0.95 * original WT growth','subsets':[]}
for n in range(1,4):
 for subset in itertools.combinations(terms,n):
  with m:
   ob=m.reactions.get_by_id(obj.id)
   ob.subtract_metabolites({m.metabolites.get_by_id(z):coef[z] for z in subset})
   patched_wt=float(m.slim_optimize())
   rescued={}
   for g in essential:
    with m:
     m.genes.get_by_id(g).knock_out();gr=float(m.slim_optimize() or 0)
     if gr>=.95*wt:rescued[g]=gr
   out['subsets'].append({'removed':subset,'edited_wt_growth':patched_wt,'n_rescued':len(rescued),'rescued_ko_growth':rescued})
  print(subset,len(rescued),flush=True)
triple=next(x for x in out['subsets'] if len(x['removed'])==3)
assert sorted(triple['rescued_ko_growth'])==previous['flips_from_fresh_manual_counterfactual']
assert all(x['n_rescued']==0 for x in out['subsets'] if len(x['removed'])==1)
pairs=[x for x in out['subsets'] if len(x['removed'])==2]
union=set().union(*(x['rescued_ko_growth'] for x in pairs))
out['triple_only_rescues']=sorted(set(triple['rescued_ko_growth'])-union)
out['scope']='post-observation exploratory interactions within a single reconstruction; no G2 repair, wet-lab claim or independent replication'
Path('results/moco_subset_interaction_ijn1463.json').write_text(json.dumps(out,indent=2)+'\n')
print('triple-only',out['triple_only_rescues'])
