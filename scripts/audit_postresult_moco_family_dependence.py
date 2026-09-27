"""Post-result family-dependence sensitivity on saved MoCo model census."""
import csv,hashlib,json
from pathlib import Path
from collections import defaultdict,Counter
survey=Path('results/bigg_model_survey.csv')
rows=list(csv.DictReader(survey.open()))
selection=[r for r in rows if r['moco_biomass']=='True']
assert len(selection)==68 and len({r['accession'] for r in selection})==68
files={}; grouped=defaultdict(list); patterns=defaultdict(list)
for r in selection:
 p=Path('results/bigg_moco_growth')/(r['accession']+'.json')
 q=json.loads(p.read_text()); assert q['accession']==r['accession'] and q['status']=='scored'
 assert q['n_flips']==len(q['flips']) and len(set(q['flips']))==len(q['flips'])
 family=r['organism'].split()[0]
 x={'accession':r['accession'],'organism':r['organism'],'n_flips':q['n_flips'],'flip_gene_ids':q['flips']}
 grouped[family].append(x); patterns[tuple(q['flips'])].append(r['accession'])
 files[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest()
assert sum(len(v) for v in grouped.values())==68
assert sum(x['n_flips']>0 for v in grouped.values() for x in v)==63
assert sum(x['n_flips'] for v in grouped.values() for x in v)==511
summary=[]
for family,v in sorted(grouped.items()):
 summary.append({'survey_genus':family,'n_model_accessions':len(v),'n_with_at_least_one_flip':sum(x['n_flips']>0 for x in v),
                 'model_gene_flip_count':sum(x['n_flips'] for x in v),
                 'n_unique_gene_id_strings':len(set(g for x in v for g in x['flip_gene_ids'])),
                 'zero_flip_model_ids':[x['accession'] for x in v if not x['n_flips']],
                 'model_ids':[x['accession'] for x in v]})
other=[x for fam,v in grouped.items() if fam not in ('Escherichia','Shigella') for x in v]
patterns_out=[{'flipped_gene_ids':list(k),'n_model_accessions':len(v),'model_ids':sorted(v)} for k,v in sorted(patterns.items(),key=lambda kv:(-len(kv[1]),kv[0]))]
out={'status':'post-result crude genus/template-dependence audit; no biological replication claim',
     'protocol':'notes/postresult_moco_family_sensitivity.md',
     'source_url':'https://bigg.ucsd.edu/api/v2/models',
     'survey_sha256':hashlib.sha256(survey.read_bytes()).hexdigest(),
     'per_model_sha256':files,
     'n_model_accessions':68,'n_with_flip':63,'total_model_gene_flip_calls':511,
     'genus_groups':summary,
     'escherichia_shigella':{'n_model_accessions':sum(len(grouped[k]) for k in ('Escherichia','Shigella')),
       'n_with_flip':sum(x['n_flips']>0 for k in ('Escherichia','Shigella') for x in grouped[k])},
     'other_genera':{'n_model_accessions':len(other),'n_with_flip':sum(x['n_flips']>0 for x in other),
       'n_genera_with_flip':sum(any(x['n_flips']>0 for x in v) for k,v in grouped.items() if k not in ('Escherichia','Shigella')),
       'model_ids':[x['accession'] for x in other]},
     'repeated_exact_flip_sets':patterns_out,
     'limits':'Genus from stated organism label is a rough proxy, not a reconstruction ancestry tree. Repeated gene-ID strings may have shared naming and reaction templates; unique IDs across genera need not be independent genes or wet-lab calls. Outcomes and strata known before this audit. No condition-matched phenotype test or G2 repair.'}
Path('results/postresult_moco_family_dependence.json').write_text(json.dumps(out,indent=2)+'\n')
print('groups',[(x['survey_genus'],x['n_model_accessions'],x['n_with_at_least_one_flip']) for x in summary])
print('other',out['other_genera']);print('most common exact sets',[(x['n_model_accessions'],x['flipped_gene_ids']) for x in patterns_out[:5]])
