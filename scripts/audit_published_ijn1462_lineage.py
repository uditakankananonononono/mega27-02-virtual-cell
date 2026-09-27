"""Outcome-aware published-model provenance audit, not phenotype validation."""
import csv,hashlib,json
from pathlib import Path
from collections import Counter
from lxml import etree
import cobra,openpyxl
s=Path('/tmp/nogales-supp/EMI-22-255-s003.xml');b=Path('data/bigg_survey/iJN1463.json')
assert hashlib.sha256(s.read_bytes()).hexdigest()=='0c8113ee6412c21df1afbd0682d1e89e73a515b558a2cab29ee8213bb5afb5f0'
expected=next(z['sha256'] for z in csv.DictReader(open('results/bigg_model_survey.csv')) if z['accession']=='iJN1463')
assert hashlib.sha256(b.read_bytes()).hexdigest()==expected
root=etree.parse(str(s));model=root.getroot()[0];blank=[c for c in model.find('{*}listOfCompartments') if c.get('id')=='']
assert len(blank)==1 and not any(z.get('compartment')=='' for z in model.find('{*}listOfSpecies'))
model.find('{*}listOfCompartments').remove(blank[0]);tmp=Path('/tmp/nogales-supp/s3-fixed-temp.xml');root.write(str(tmp),encoding='utf-8',xml_declaration=True)
a=cobra.io.read_sbml_model(str(tmp));c=cobra.io.load_json_model(str(b))
def ids(x,t):return {z.id for z in getattr(x,t)}
out={'status':'post-result source-lineage audit, not prospective validation','protocol':'notes/postresult_published_ijn1462_lineage_audit.md',
     'published_s3_url':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7078882/supplementaryFiles',
     'published_s3_sha256':hashlib.sha256(s.read_bytes()).hexdigest(),'published_s4_sha256':'73afee37058fab4ef254f3d0acd83fd9eeac9f48939d832a7f94c78b5bd7e7ab','published_s7_sha256':'3788cb2009b295c54301fa939cb28d107e317b90f1835737d465860b49cca322','bigg_url':'https://bigg.ucsd.edu/static/models/iJN1463.json',
     'bigg_sha256':expected,'s3_read_fix':'Removed only one empty unreferenced XML compartment from temporary copy; original source unmodified.',
     'models':[{'id':x.id,'genes':len(x.genes),'reactions':len(x.reactions),'metabolites':len(x.metabolites)} for x in [a,c]],'sets':{}}
for t in ('genes','reactions','metabolites'):
 x,y=ids(a,t),ids(c,t);out['sets'][t]={'shared':len(x&y),'published_only':sorted(x-y),'bigg_only':sorted(y-x)}
shared=ids(a,'reactions')&ids(c,'reactions');bad=[]
for id in sorted(shared):
 x=a.reactions.get_by_id(id);y=c.reactions.get_by_id(id)
 diff=[]
 if {z.id:float(v) for z,v in x.metabolites.items()}!={z.id:float(v) for z,v in y.metabolites.items()}:diff.append('stoichiometry')
 if x.bounds!=y.bounds:diff.append('bounds')
 if x.gene_reaction_rule!=y.gene_reaction_rule:diff.append('gpr')
 if diff:bad.append({'id':id,'fields':diff,'published_bounds':x.bounds,'bigg_bounds':y.bounds,'published_gpr':x.gene_reaction_rule,'bigg_gpr':y.gene_reaction_rule})
out['shared_reaction_changes']=bad;out['shared_reaction_change_counts']=dict(Counter(f for r in bad for f in r['fields']))
obj1=next(z for z in a.reactions if z.objective_coefficient);obj2=next(z for z in c.reactions if z.objective_coefficient)
q=lambda r:{z.id:float(v) for z,v in r.metabolites.items()}
out['biomass_objective']={'published_id':obj1.id,'bigg_id':obj2.id,'exact_same_stoichiometry':q(obj1)==q(obj2),'published_bounds':obj1.bounds,'bigg_bounds':obj2.bounds,'diff':{k:[q(obj1).get(k),q(obj2).get(k)] for k in set(q(obj1))|set(q(obj2)) if q(obj1).get(k)!=q(obj2).get(k)}}
mo=json.load(open('results/moco_subset_interaction_ijn1463.json'))
mo_ids=set(next(z for z in mo['subsets'] if z['removed']==['bmocogdp_c','mocogdp_c'])['rescued_ko_growth'])
ex=json.load(open('results/postresult_pair_exception_reactions.json'))
pair_ids=set(next(z for z in ex['exceptions'] if z['pair']==['2fe2s_c','btamp_c'])['lost_union_genes'])
assert len(mo_ids)==6 and len(pair_ids)==8
out['target_gene_presence']={'moco_six':sorted(mo_ids),'suppressive_eight':sorted(pair_ids),'all_fourteen_in_published':(mo_ids|pair_ids)<=ids(a,'genes')}
for fn,label in [('EMI-22-255-s007.xlsx','S4_essentiality'),('EMI-22-255-s010.xlsx','S7_barseq')]:
 p=Path('/tmp/nogales-supp')/fn;w=openpyxl.load_workbook(p,read_only=True,data_only=True);ws=w.active
 target={id:[] for id in mo_ids|pair_ids}
 for row in ws.iter_rows(values_only=True):
  for id in target:
   if id in row:target[id].append([i+1 for i,z in enumerate(row) if z==id])
 out[label]={'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'sheet_rows':ws.max_row,'sheet_columns':ws.max_column,'target_occurrence_columns':{k:v for k,v in target.items() if v}}
out['limits']='Reused published validation tables and close but nonidentical model lineage; no independent phenotype score; no G2 gate change.'
Path('results/postresult_published_ijn1462_lineage.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print('models',out['models'],'shared genes',out['sets']['genes']['shared'],'shared reactions',out['sets']['reactions']['shared'],'changes',out['shared_reaction_change_counts'])
print('biomass same',out['biomass_objective']['exact_same_stoichiometry'],'targets in published',out['target_gene_presence']['all_fourteen_in_published'])
print('S4 occurrences',out['S4_essentiality']['target_occurrence_columns']);print('S7 occurrences',out['S7_barseq']['target_occurrence_columns'])
