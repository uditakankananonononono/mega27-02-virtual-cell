"""Post-result named-gene overlap and contradiction caution for iYO844 objective rescues."""
import csv,hashlib,json
from pathlib import Path
import cobra
src=Path('data/bigg_survey/iYO844.json');row=next(r for r in csv.DictReader(open('results/bigg_model_survey.csv')) if r['accession']=='iYO844')
assert hashlib.sha256(src.read_bytes()).hexdigest()==row['sha256']
m=cobra.io.load_json_model(str(src));res=json.load(open('results/biomass_transfer_iyo844.json'))
assert res['source_sha256']==row['sha256'] and res['n_terms']==60
source=Path('data/external/bsub/subtiwiki_essential_names.json');ess={z.casefold() for z in json.load(open(source))}
model_names={g.id:(g.name or '').strip() for g in m.genes};names_to_ids={}
for gid,name in model_names.items():
 if name:names_to_ids.setdefault(name.casefold(),[]).append(gid)
base=cobra.flux_analysis.single_gene_deletion(m,processes=1)
original={next(iter(ids)) for ids,gr in zip(base['ids'],base['growth']) if gr is not None and gr<1e-6}
assert len(original)==res['baseline_essential_count']
rescued=set().union(*(r['rescued_ko_growth'] for r in res['rows']))
assert rescued<=original

def map_ids(ids):
 unique=[];ambiguous=[];unnamed=[]
 for gid in sorted(ids):
  name=model_names[gid]
  if not name:unnamed.append(gid)
  elif len(names_to_ids[name.casefold()])!=1:ambiguous.append({'id':gid,'name':name,'ids':names_to_ids[name.casefold()]})
  else:unique.append({'id':gid,'name':name,'in_subtiwiki_essential_list':name.casefold() in ess})
 return {'n_source_gene_ids':len(ids),'n_unique_name_mapped':len(unique),'n_unnamed':len(unnamed),'unnamed_ids':unnamed,
         'n_ambiguous_name':len(ambiguous),'ambiguous':ambiguous,
         'n_in_subtiwiki_essential_list':sum(z['in_subtiwiki_essential_list'] for z in unique),'unique_rows':unique}
terms=[]
for r in res['rows']:
 if not r['rescued_ko_growth']:continue
 p=map_ids(set(r['rescued_ko_growth']));terms.append({'metabolite_id':r['metabolite_id'],
    'n_rescued':len(r['rescued_ko_growth']),'n_unique_names':p['n_unique_name_mapped'],
    'n_subtiwiki_essential_names':p['n_in_subtiwiki_essential_list'],
    'subtiwiki_essential_names':[z['name'] for z in p['unique_rows'] if z['in_subtiwiki_essential_list']],
    'ambiguous_names':p['ambiguous'],'unnamed_ids':p['unnamed_ids']})
assert len(terms)==res['n_positive']==17
out={'status':'post-result cross-list essential-name overlap; not independent or condition-matched validation',
 'protocol':'notes/postresult_bsub_objective_lab_check.md','bigg_source_url':row['source_url'],'bigg_sha256':row['sha256'],
 'subtiwiki_source_url':'https://subtiwiki.uni-goettingen.de/wiki/index.php?title=Essential_genes',
 'subtiwiki_local_html_sha256':hashlib.sha256(Path('data/external/bsub/subtiwiki_essential_genes.html').read_bytes()).hexdigest(),
 'subtiwiki_names_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'n_subtiwiki_names':len(ess),
 'baseline_model_essential':map_ids(original),'rescued_union':map_ids(rescued),
 'nonzero_terms':terms,'n_terms':res['n_terms'],
 'limits':'SubtiWiki category page combines sources; its stated 251 protein and 2 RNA essential genes do not equal the 260 parsed names. iYO844 reconstruction incorporated gene-essentiality evidence, so source independence is unproven. Strain, medium and assay may differ. An overlap warns against assuming model rescue is viability; absence does not prove rescue. Names can be incomplete or ambiguous; no p-value or gate inference is valid.'}
Path('results/bsub_rescue_subtiwiki_overlap.json').write_text(json.dumps(out,indent=2)+'\n')
for k in ['baseline_model_essential','rescued_union']:
 z=out[k];print(k,z['n_source_gene_ids'],z['n_unique_name_mapped'],z['n_in_subtiwiki_essential_list'],z['n_ambiguous_name'],z['n_unnamed'])
for t in terms:print(t['metabolite_id'],t['n_rescued'],t['n_subtiwiki_essential_names'],','.join(t['subtiwiki_essential_names']))
