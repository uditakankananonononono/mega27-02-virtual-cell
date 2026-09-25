"""Exploratory independent annotations for the iML1515 AHCYSNS route.

QuickGO, Pfam/InterPro and PDB do not establish cellular uptake of exogenous
SAM or validate an in-vivo rescue. This is protein/reaction triangulation only.
"""
import json,urllib.request,hashlib
from pathlib import Path
P='P0AF12'; HEAD={'User-Agent':'VC2-virtual-cell/1.0 (research annotation check)','Accept':'application/json'}
URLS={
 'quickgo':f'https://www.ebi.ac.uk/QuickGO/services/annotation/search?geneProductId=UniProtKB:{P}&taxonId=83333&limit=100',
 'pfam_interpro':f'https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/{P}/?page_size=100',
 'rcsb_pdb':'https://data.rcsb.org/rest/v1/core/entry/1JYS',
 'panther':'https://www.pantherdb.org/services/oai/pantherdb/geneinfo?geneInputList=P0AF12&organism=83333'}
def get(u):
 with urllib.request.urlopen(urllib.request.Request(u,headers=HEAD),timeout=30) as r:return json.load(r)
q=get(URLS['quickgo']);pf=get(URLS['pfam_interpro']);pdb=get(URLS['rcsb_pdb']);pan=get(URLS['panther'])
assert q['numberOfHits']==len(q['results']) and q['pageInfo']['total']==1
assert {x['geneProductId'] for x in q['results']}=={f'UniProtKB:{P}'}
assert pf['count']==len(pf['results']) and not pf['next']
assert pdb['rcsb_id']=='1JYS' and 'nucleosidase' in pdb['struct']['title'].lower()
# Cross-check against saved UniProt-derived Pfam/PDB tables actually used in v6/v8.
pf_ids={x['metadata']['accession'] for x in pf['results']}
saved_pf=next(x.split('\t')[1].strip().strip(';').split(';') for x in Path('data/external/uniprot_ecoli_pfam.tsv').read_text().splitlines() if x.startswith(P+'\t'))
saved_pdb=next(x.split('\t')[1].strip().strip(';').split(';') for x in Path('data/external/uniprot_ecoli_pdb.tsv').read_text().splitlines() if x.startswith(P+'\t'))
assert pf_ids&set(saved_pf) and '1JYS' in saved_pdb
assert any(x.get('goAspect')=='molecular_function' for x in q['results'])
gene=pan['search']['mapped_genes']['gene']
assert gene['mapped_id_list']==P and 'UniProtKB='+P in gene['accession']
assert 'NUCLEOSIDASE' in gene['family_name'].upper()
result={'design':'exploratory annotation triangulation of gene mtnN/b0159, not independent phenotype validation',
 'source_urls':URLS,'protein_accession':P,
 'quickgo':{'annotation_count':len(q['results']), 'molecular_function_go_ids':sorted({x['goId'] for x in q['results'] if x['goAspect']=='molecular_function'}),
            'direct_assigned_by':sorted({str(x['assignedBy']) for x in q['results']})},
 'pfam':{'live_ids':sorted(pf_ids),'saved_feature_ids':saved_pf,'shared_ids':sorted(pf_ids&set(saved_pf))},
 'panther':{'version':pan['search']['product']['version'],'family_id':gene['family_id'],
            'subfamily_id':gene['sf_id'],'family_name':gene['family_name']},
 'pdb':{'live_id':pdb['rcsb_id'],'title':pdb['struct']['title'],'resolution_angstrom':pdb.get('rcsb_entry_info',{}).get('resolution_combined'),
        'saved_feature_contains_id':True},
 'input_sha256':{x:hashlib.sha256(Path(x).read_bytes()).hexdigest() for x in ['data/external/uniprot_ecoli_pfam.tsv','data/external/uniprot_ecoli_pdb.tsv']},
 'caveat':'These sources annotate the enzyme and structure. They do not show that exogenous SAM is imported or that model-predicted purine rescue occurs in cells.'}
Path('results/sah_orthogonal_annotations.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
