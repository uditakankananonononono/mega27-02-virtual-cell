"""Independent manual-style second-family counterfactual, without reading census result for targets.
Recomputes full baseline essential set and targeted Pseudomonas objective repair; no result tuning.
"""
import csv,hashlib,json
from pathlib import Path
import cobra
from scipy.optimize import linprog
p=Path('data/bigg_survey/iJN1463.json'); survey={r['accession']:r for r in csv.DictReader(open('results/bigg_model_survey.csv'))};row=survey['iJN1463']
assert hashlib.sha256(p.read_bytes()).hexdigest()==row['sha256']
m=cobra.io.load_json_model(str(p));obj=m.reactions.get_by_id('BIOMASS_KT2440_WT3')
original={z.id:float(c) for z,c in obj.metabolites.items() if z.id in ['bmocogdp_c','mobd_c','mocogdp_c']}
assert len(original)==3 and all(v<0 for v in original.values())
wt=float(m.slim_optimize());assert wt>1e-6
base=cobra.flux_analysis.single_gene_deletion(m,processes=1)
E=sorted(sorted(g)[0] for g,gr in zip(base['ids'],base['growth']) if gr is not None and gr<1e-6)
# Remove precisely the 3 objective metabolites, not all metabolites sharing name; independent explicit selection.
for name in original:
 obj.subtract_metabolites({m.metabolites.get_by_id(name):original[name]})
patched_wt=float(m.slim_optimize());assert patched_wt>0
flips=[];profiles={}
for gene in E:
 with m:
  m.genes.get_by_id(gene).knock_out();g=float(m.slim_optimize() or 0)
  if g>=.95*wt:
   flips.append(gene);profiles[gene]={'patched_growth':g,'vs_original_wt':g/wt}
# independently recompute baseline growth for flipped genes in a freshly loaded model, to avoid context contamination
b=cobra.io.load_json_model(str(p));orig_growth={}
for gene in flips:
 with b:
  b.genes.get_by_id(gene).knock_out();orig_growth[gene]=float(b.slim_optimize() or 0)
  assert orig_growth[gene]<1e-6
out={'model':'Pseudomonas putida KT2440 iJN1463','source_url':row['source_url'],'snapshot_sha256':row['sha256'],
 'objective':'BIOMASS_KT2440_WT3','objective_negative_moco_coefficients':original,
 'baseline_wt_growth':round(wt,10),'patched_wt_growth':round(patched_wt,10),'baseline_essential_genes':len(E),
 'flips_from_fresh_manual_counterfactual':flips,'baseline_ko_growth_for_flips':{g:round(v,10) for g,v in orig_growth.items()},
 'patched_ko_growth_for_flips':{g:{k:round(v,10) for k,v in z.items()} for g,z in profiles.items()},'independently_computed':True,
 'scope':'in-silico objective-dependence only, no in-vivo viability proof; separate organism family from E. coli iJO1366'}
Path('results/bigg_second_family_verification.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k not in ['patched_ko_growth_for_flips','baseline_ko_growth_for_flips']},indent=2))
