"""Source-bound, post-result reaction/GPR annotation for all saved pair exceptions."""
import csv,hashlib,json
from pathlib import Path
import cobra
model_path=Path('data/bigg_survey/iJN1463.json')
sha=hashlib.sha256(model_path.read_bytes()).hexdigest()
ledger=next(row for row in csv.DictReader(open('results/bigg_model_survey.csv')) if row['accession']=='iJN1463')
assert sha==ledger['sha256']=='8232a1ea39af9b9b0313601a3cefca1028e50e53b2ca2c06002c1a9adf7568bf'
m=cobra.io.load_json_model(str(model_path))
saved=json.load(open('results/postresult_pair_interactions.json'))
exceptions=[x for x in saved['rows'] if x['pair_only_genes'] or x['single_union_only_genes']]
assert len(exceptions)==saved['n_pair_only']+saved['n_lost_single_union']==3
rows=[]
for x in exceptions:
 metabolites=[]
 for mid in x['pair']:
  met=m.metabolites.get_by_id(mid)
  metabolites.append({'metabolite_id':met.id,'name':met.name,'formula':met.formula,
    'gene_associated_reactions':[{'id':r.id,'name':r.name,'gene_reaction_rule':r.gene_reaction_rule,'stoichiometric_coefficient':float(r.metabolites[met])} for r in sorted(met.reactions,key=lambda r:r.id) if r.genes and met in r.metabolites]})
 geneids=sorted(set(x['pair_only_genes']+x['single_union_only_genes']))
 genes=[]
 for gid in geneids:
  g=m.genes.get_by_id(gid)
  genes.append({'gene_id':g.id,'source_model_name':g.name,'model_reactions':[{'id':r.id,'name':r.name,'gene_reaction_rule':r.gene_reaction_rule} for r in sorted(g.reactions,key=lambda r:r.id)]})
 rows.append({'pair':x['pair'],'single_counts':[x['single_a_count'],x['single_b_count']],
   'union_count':x['union_count'],'pair_count':x['pair_count'],'pair_only_genes':x['pair_only_genes'],
   'lost_union_genes':x['single_union_only_genes'],'removed_metabolites':metabolites,'affected_genes':genes})
out={'status':'post-result model-internal reaction/GPR annotation, not flux or phenotype evidence',
 'protocol':'notes/postresult_pair_exception_reaction_audit.md','model_url':ledger['source_url'],
 'model_sha256':sha,'saved_pair_result':'results/postresult_pair_interactions.json','n_exceptions':len(rows),
 'exceptions':rows,'limits':'No knockout-specific flux or dual proof; reaction association does not establish causal route; all three are selected after result; G2 still failed.'}
Path('results/postresult_pair_exception_reactions.json').write_text(json.dumps(out,indent=2)+'\n')
print([(x['pair'],[(g['gene_id'],g['source_model_name'],[r['id'] for r in g['model_reactions']]) for g in x['affected_genes']]) for x in rows])
