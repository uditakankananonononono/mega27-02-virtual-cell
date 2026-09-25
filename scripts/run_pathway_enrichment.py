"""Offline gseapy overrepresentation for false high-priority essentiality calls."""
import collections,hashlib,json
from pathlib import Path
import pandas as pd
import gseapy as gp
A=pd.read_csv('results/ensemble_v2_oof.csv')
assert len(A)==1249 and set(A.essential.unique())=={0,1}
universe=set(A.bnumber)
query=A[A.essential==0].nlargest(100,'v2_lr').bnumber.tolist()
sets=collections.defaultdict(set)
for line in open('data/external/kegg_eco_pathway.tsv'):
    g,p=line.rstrip('\n').split('\t');sets[p.replace('path:eco','')].add(g.replace('eco:',''))
gsets={p:sorted(g&universe) for p,g in sets.items() if 5<=len(g&universe)<=500}
assert len(gsets)>10
name={x.split('\t')[0].replace('eco',''):x.split('\t')[1].strip() for x in open('data/external/kegg_pathway_names.tsv') if '\t' in x}
r=gp.enrich(gene_list=query,gene_sets=gsets,background=sorted(universe),outdir=None,no_plot=True,verbose=False)
d=r.results.sort_values('Adjusted P-value')
# gseapy output column spellings are explicit here, avoiding tuple identifier rewrite.
out=[]
for _,x in d.iterrows():
    out.append({'pathway_id':str(x['Term']),'pathway_name':name.get(str(x['Term']),''),
                'overlap':str(x['Overlap']),'p_value':float(x['P-value']),
                'fdr_bh':float(x['Adjusted P-value']), 'genes':str(x['Genes'])})
summary={'question':'KEGG overrepresentation of top 100 v2 false essentiality priorities',
         'n_universe':len(universe),'n_query':len(query),'n_tested_pathways':len(gsets),
         'n_returned_pathways':len(out),'n_fdr_lt_0_05':sum(x['fdr_bh']<0.05 for x in out),
         'top10':out[:10],
         'input_sha256':{p:hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in ['results/ensemble_v2_oof.csv','data/external/kegg_eco_pathway.tsv']},
         'caveat':'Exploratory ranked false positives; pathway overlap and selection can inflate interpretation. Offline saved KEGG map, not a new accession.'}
Path('results/pathway_enrichment.json').write_text(json.dumps(summary,indent=2)+'\n')
pd.DataFrame(out).to_csv('results/pathway_enrichment.csv',index=False)
print(json.dumps(summary,indent=2))
