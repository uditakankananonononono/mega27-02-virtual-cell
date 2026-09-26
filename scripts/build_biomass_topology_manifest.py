"""Build pre-perturbation features only; no target/outcome score loading."""
import csv,hashlib,json
from pathlib import Path
import cobra
p=Path('data/bigg_survey/iJN1463.json');row=next(r for r in csv.DictReader(open('results/bigg_model_survey.csv')) if r['accession']=='iJN1463')
sha=hashlib.sha256(p.read_bytes()).hexdigest();assert sha==row['sha256']
m=cobra.io.load_json_model(str(p));o=m.reactions.get_by_id('BIOMASS_KT2440_WT3');wt=m.optimize();assert wt.status=='optimal' and wt.objective_value>0
seen={z['removed'][0] for z in json.load(open('results/moco_specificity_null_ijn1463.json'))['controls']}
terms=[]
for met,c in sorted(o.metabolites.items(),key=lambda x:x[0].id):
 if c>=0:continue
 links=[r for r in met.reactions if r.id!=o.id]
 producer=[r for r in links if r.metabolites[met]*float(wt.fluxes[r.id])>1e-9]
 terms.append({'metabolite_id':met.id,'name':met.name,'compartment':met.compartment,
               'abs_biomass_coefficient':abs(float(c)),
               'non_biomass_reaction_count':len(links),
               'gene_associated_reaction_count':sum(bool(r.genes) for r in links),
               'exchange_uptake_allowed':any(r.boundary and r.lower_bound<0 for r in links),
               'active_producer_reaction_count_wt':len(producer),
               'active_producer_abs_flux_sum_wt':round(sum(abs(float(wt.fluxes[r.id])) for r in producer),8),
               'outcome_status':'previously_scored_control' if met.id in seen else ('previously_scored_moco_bundle_only' if met.id in {'bmocogdp_c','mobd_c','mocogdp_c'} else 'single_term_not_scored_in_prior_48_control')})
assert len(terms)==102 and sum(t['outcome_status']=='previously_scored_control' for t in terms)==48
out={'model':'iJN1463','source_url':row['source_url'],'sha256':sha,'objective':o.id,'wt_growth':round(float(wt.objective_value),10),
     'status':'feature inventory; no new individual-term outcome scoring','n_terms':len(terms),'terms':terms,
     'limitations':'reaction-degree and active-flux features are proxies, not causal redundancy; 48 single-term outcomes were already seen; no held-out claim'}
Path('results/biomass_topology_manifest.json').write_text(json.dumps(out,indent=2)+'\n')
print('model',out['model'],'terms',len(terms),'observed controls',len(seen),'unscored singles',sum(t['outcome_status']=='single_term_not_scored_in_prior_48_control' for t in terms))
