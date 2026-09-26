"""Exploratory single-coefficient outcomes missing from prior 48-control run.
Writes each completed result immediately so interrupted work resumes without duplicating outcomes.
"""
import csv,hashlib,json,time
from pathlib import Path
import cobra
p=Path('data/bigg_survey/iJN1463.json');src=next(z for z in csv.DictReader(open('results/bigg_model_survey.csv')) if z['accession']=='iJN1463');assert hashlib.sha256(p.read_bytes()).hexdigest()==src['sha256']
manifest=json.load(open('results/biomass_topology_manifest.json'));assert manifest['sha256']==src['sha256'] and manifest['n_terms']==102
prior=json.load(open('results/moco_specificity_null_ijn1463.json'))
known={z['removed'][0]:z for z in prior['controls']}
remaining=[z['metabolite_id'] for z in manifest['terms'] if z['metabolite_id'] not in known]
assert len(known)==48 and len(remaining)==54
m=cobra.io.load_json_model(str(p));o=m.reactions.get_by_id('BIOMASS_KT2440_WT3');coef={z.id:float(c) for z,c in o.metabolites.items() if c<0}
wt=float(m.slim_optimize());assert round(wt,10)==prior['original_wt_growth'];essential=sorted(set(prior['positive_bundle']['rescued_genes'])|set(g for z in known.values() for g in z['rescued_genes']))
# Reconstruct full original zero-growth essential set, never restrict to previously rescued genes.
base=cobra.flux_analysis.single_gene_deletion(m,processes=1);full=sorted(next(iter(ids)) for ids,gr in zip(base['ids'],base['growth']) if gr is not None and gr<1e-6);assert len(full)==262 and set(essential)<=set(full)
checkpoint=Path('results/biomass_remaining_singles_progress.json');done=json.loads(checkpoint.read_text()) if checkpoint.exists() else {'model_sha256':src['sha256'],'original_wt':round(wt,10),'baseline_essential_genes':len(full),'results':{}}
assert done['model_sha256']==src['sha256'] and done['original_wt']==round(wt,10)
for z in remaining:
 if z in done['results']:continue
 with m:
  m.reactions.get_by_id('BIOMASS_KT2440_WT3').subtract_metabolites({m.metabolites.get_by_id(z):coef[z]})
  patched_wt=float(m.slim_optimize());flips=[]
  for gene in full:
   with m:
    m.genes.get_by_id(gene).knock_out();v=float(m.slim_optimize() or 0)
    if v>=.95*wt:flips.append(gene)
 done['results'][z]={'coefficient':coef[z],'patched_wt_growth':round(patched_wt,10),'n_baseline_essential_rescued':len(flips),'rescued_genes':flips}
 checkpoint.write_text(json.dumps(done,indent=2)+'\n');print(len(done['results']),'/',len(remaining),z,len(flips),flush=True)
# All 102 terms with uniform single-term intervention and manifest features.
allrows=[]
for row in manifest['terms']:
 z=row['metabolite_id'];out=done['results'][z] if z in done['results'] else known[z]
 allrows.append({'metabolite_id':z,'prior_status':row['outcome_status'],
                 'n_baseline_essential_rescued':out['n_baseline_essential_rescued'],
                 'rescued_genes':out['rescued_genes']})
assert len(allrows)==102
final={'status':'exploratory all-term iJN1463 biomass single-coefficient sensitivity, not external validation',
 'model_sha256':src['sha256'],'n_terms':102,'n_original_essential':len(full),'n_prior_controls':48,
 'n_new_singles':54,'rows':allrows,
 'limits':'48 outcomes already seen; the model is one reconstruction; no preregistered predictive feature model, independent family test, wet-lab finding or G2 repair'}
Path('results/biomass_single_term_sensitivity.json').write_text(json.dumps(final,indent=2)+'\n')
print('COMPLETE terms',len(allrows),'nonzero',sum(z['n_baseline_essential_rescued']>0 for z in allrows),'max',max(z['n_baseline_essential_rescued'] for z in allrows))
