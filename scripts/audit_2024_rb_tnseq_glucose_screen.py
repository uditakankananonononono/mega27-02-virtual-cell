"""Frozen metadata-conditioned source screen of 2024 compendium; exploratory, not G2."""
import hashlib,json,statistics
from pathlib import Path
import openpyxl
source=Path('/tmp/fModule_Metadata.xlsx');assert hashlib.sha256(source.read_bytes()).hexdigest()=='4d649385ac06684482396a125f135df22a2a5060da73485b2cd14468f8cc8be1'
w=openpyxl.load_workbook(source,read_only=True,data_only=True)
meta=list(w['metadata'].values); h=meta[0]; s=[dict(zip(h,row)) for row in meta[1:] if row[5]=='carbon source' and row[23]=='D-Glucose' and row[26] is None]
assert len(s)==13
mo=json.load(open('results/moco_subset_interaction_ijn1463.json'));mo_ids=sorted(next(z for z in mo['subsets'] if z['removed']==['bmocogdp_c','mocogdp_c'])['rescued_ko_growth'])
ex=json.load(open('results/postresult_pair_exception_reactions.json'));pair_ids=sorted(next(z for z in ex['exceptions'] if z['pair']==['2fe2s_c','btamp_c'])['lost_union_genes'])
assert len(mo_ids)==6 and len(pair_ids)==8 and len(set(mo_ids+pair_ids))==14
f=w['fitness_measurements'];t=w['T-like_statistics']; fh=next(f.values);th=next(t.values); assert fh==th and len(fh)==337
idx={x['expName']:next(i for i,v in enumerate(fh) if isinstance(v,str) and v.startswith(x['expName']+' ')) for x in s}
chosen={}
for sheet,key in [(f,'fitness'),(t,'t_like')]:
 for r in sheet.iter_rows(values_only=True):
  g=r[1]
  if g in mo_ids+pair_ids:chosen.setdefault(g,{})[key]={exp:r[i] for exp,i in idx.items()}
assert len(chosen)==11 and all(len(d)==2 for d in chosen.values())
missing=sorted(set(mo_ids+pair_ids)-set(chosen));assert missing==['PP_0437','PP_1602','PP_3457']
out={'status':'exploratory post-result cohort screen, NOT independent validation or G2 pass','protocol':'notes/prereg_2024_rb_tnseq_condition_ledger.md','source':'https://github.com/beckham-lab/fModule','source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'missing_panel_genes':missing,'samples':[{k:str(x[k]) for k in ['expName','dateStarted','timeZeroSet','set','expDesc','expDescLong','mutantLibrary','person','Inoculum media type','media','rep','total_rep','condition_1','concentration_1','condition_2']} for x in s],'panels':{},'limits':'2020 assay independence from 2019 reconstruction curation/validation unverified; 2022 set101 has internal date conflict; 2017 set5 predates reconstruction; no growth survival endpoint; library polarity and gene coverage remain considerations. -2.7 was outcome-tuned in older published S7 and must not be ported as a validated threshold.'}
for name,ids in [('moco',mo_ids),('suppressive_pair',pair_ids)]:
 out['panels'][name]={}
 for g in ids:
  d=chosen.get(g);
  if d is None:out['panels'][name][g]={'status':'not represented in workbook'};continue
  groups={}
  for x in s:
   label=x['set']+'_'+str(x['concentration_1'])+'mM_'+x['media'];groups.setdefault(label,[]).append(x['expName'])
  out['panels'][name][g]={'measurements':d,'groups':{k:{'exp_names':v,'fitness_median':statistics.median(d['fitness'][n] for n in v) if all(isinstance(d['fitness'][n],(float,int)) for n in v) else None,'t_like_values':[d['t_like'][n] for n in v]} for k,v in groups.items()}}
Path('results/postresult_2024_rb_tnseq_glucose_screen.json').write_text(json.dumps(out,indent=2,sort_keys=True,default=str)+'\n')
print('SOURCE',out['source_sha256'],'SAMPLES',len(s),'PANEL',len(chosen))
for panel,genes in out['panels'].items():
 print(panel)
 for g,d in genes.items():print(g,[(k,round(v['fitness_median'],3) if v['fitness_median'] is not None else None) for k,v in d.get('groups',{}).items()] if 'groups' in d else 'MISSING')
