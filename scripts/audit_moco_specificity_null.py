"""Post-result one-model coefficient-matched single-term biomass null; no assay labels."""
import csv,hashlib,json,time
from pathlib import Path
import cobra
from cobra.flux_analysis import single_gene_deletion
p=Path('data/bigg_survey/iJN1463.json')
row=next(r for r in csv.DictReader(open('results/bigg_model_survey.csv')) if r['accession']=='iJN1463')
assert hashlib.sha256(p.read_bytes()).hexdigest()==row['sha256']
m=cobra.io.load_json_model(str(p));o=m.reactions.get_by_id('BIOMASS_KT2440_WT3')
wt=float(m.slim_optimize());base=single_gene_deletion(m,processes=1)
essential=sorted(next(iter(g)) for g,gr in zip(base['ids'],base['growth']) if gr is not None and gr<1e-6)
mo={'bmocogdp_c','mobd_c','mocogdp_c'};coeff={z.id:float(c) for z,c in o.metabolites.items() if c<0};assert mo<=coeff.keys()
sel=sorted(z for z,c in coeff.items() if z not in mo and any(abs(coeff[a])/2<=abs(c)<=2*abs(coeff[a]) for a in mo))
assert len(sel)==48 and len(essential)==262

def test(ids):
 with m:
  ob=m.reactions.get_by_id('BIOMASS_KT2440_WT3')
  for z in ids:ob.subtract_metabolites({m.metabolites.get_by_id(z):coeff[z]})
  patched_wt=float(m.slim_optimize());flips=[]
  for gene in essential:
   with m:
    m.genes.get_by_id(gene).knock_out();growth=float(m.slim_optimize() or 0)
    if growth>=.95*wt:flips.append(gene)
 return {'removed':ids,'patched_wt_growth':round(patched_wt,10),'n_baseline_essential_rescued':len(flips),'rescued_genes':flips}
t0=time.time();controls=[]
for z in sel:
 controls.append({'coefficient':coeff[z],**test([z])});print(z,controls[-1]['n_baseline_essential_rescued'],flush=True)
positive=test(sorted(mo))
out={'status':'post-result Pseudomonas iJN1463 single-model objective-term specificity scout',
 'source_url':row['source_url'],'source_sha256':row['sha256'],'objective':'BIOMASS_KT2440_WT3',
 'original_wt_growth':round(wt,10),'baseline_essential_genes':len(essential),
 'matching_rule':'each non-MoCo negative objective coefficient within factor 2 of at least one of three MoCo coefficients',
 'positive_bundle':positive,'controls':controls,
 'n_controls_with_rescue_ge_moco_bundle':sum(r['n_baseline_essential_rescued']>=positive['n_baseline_essential_rescued'] for r in controls),
 'limits':'single-term controls against three-term MoCo bundle; category/connectivity and cross-model control matching absent; in-silico only; cannot repair G2 or claim biological fitness',
 'elapsed_seconds':round(time.time()-t0,2)}
Path('results/moco_specificity_null_ijn1463.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v if k not in ['controls','positive_bundle'] else (len(v) if k=='controls' else {'n_rescued':v['n_baseline_essential_rescued'],'genes':v['rescued_genes']}) for k,v in out.items()},indent=2))
