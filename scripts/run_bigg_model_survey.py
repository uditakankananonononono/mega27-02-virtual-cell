"""Public BiGG v2 model snapshot survey with reproducible accession+hash ledger.

Run with --fetch to make bounded network requests; reruns reuse saved snapshots.
No inferred or failed accession is counted. Stored-medium FBA is descriptive.
"""
import concurrent.futures as cf
import csv, hashlib, json, sys, urllib.request
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
URL='https://bigg.ucsd.edu/api/v2/models'
BASE='https://bigg.ucsd.edu/static/models/'
DATA=ROOT/'data'/'bigg_survey';DATA.mkdir(parents=True,exist_ok=True)
OUT=ROOT/'results';OUT.mkdir(exist_ok=True)
def fetch(url):
    with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'VC2-research/1.0 (public data reproducibility)'}),timeout=45) as r:
        return r.read()
def digest(buf):return hashlib.sha256(buf).hexdigest()
listing=DATA/'model_list.json'
if '--fetch' in sys.argv and not listing.exists():
    buf=fetch(URL);json.loads(buf);listing.write_bytes(buf)
if not listing.exists():raise RuntimeError('missing public model list; run with --fetch')
items=json.loads(listing.read_bytes())['results']
assert len(items)==len({i['bigg_id'] for i in items})
def obtain(it):
    ident=it['bigg_id'];url=BASE+ident+'.json';p=DATA/(ident+'.json')
    if not p.exists():
        if '--fetch' not in sys.argv:return (it,url,None,'snapshot missing')
        try:
            b=fetch(url);j=json.loads(b)
            if j.get('id')!=ident or not j.get('reactions'):raise ValueError('model identity or reactions mismatch')
            p.write_bytes(b)
        except Exception as e:return (it,url,None,type(e).__name__+': '+str(e)[:120])
    return (it,url,p,'')
with cf.ThreadPoolExecutor(max_workers=4) as pool:got=list(pool.map(obtain,items))
rows=[]
for it,url,p,error in got:
    id=it['bigg_id'];rec={'accession':id,'organism':it['organism'],'source_url':url,'list_url':URL,
                           'path':str(p.relative_to(ROOT)) if p else '',
                           'sha256':digest(p.read_bytes()) if p else '', 'bytes':p.stat().st_size if p else 0,
                           'status':'failed' if error else 'fetched','error':error,
                           'genes':'','reactions':'','metabolites':'','biomass_reactions':'',
                           'moco_biomass':'','moco_metabolites':'','growth_on_stored_medium':'not_scored'}
    if p:
        try:
            j=json.loads(p.read_bytes());assert j.get('id')==id
            rec.update(genes=len(j['genes']),reactions=len(j['reactions']),metabolites=len(j['metabolites']))
            mm={m['id']:m.get('name','') for m in j['metabolites']}
            bio=[r for r in j['reactions'] if float(r.get('objective_coefficient',0))>0]
            ids=sorted({m for r in bio for m,coef in r['metabolites'].items()
                        if float(coef)<0 and any(t in (m+' '+mm.get(m,'')).lower() for t in ('moco','molybdo','bmoco','mobd'))})
            rec['biomass_reactions']=';'.join(r['id'] for r in bio)
            rec['moco_biomass']=bool(ids)
            rec['moco_metabolites']=';'.join(ids)
            # JSON structure is the primary survey. FBA by organism-specific
            # stored media is not comparable and is not required for this gate.
            rec['status']='analysed'
        except Exception as e:
            rec['status']='failed_analysis';rec['error']=type(e).__name__+': '+str(e)[:180]
    rows.append(rec)
fields=list(rows[0]);path=OUT/'bigg_model_survey.csv'
with path.open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
analysed=[r for r in rows if r['status']=='analysed']
summary={'source_url':URL,'source_list_sha256':digest(listing.read_bytes()),
         'listed':len(items),'fetched_and_analysed':len(analysed),
         'failed':[{k:r[k] for k in ('accession','status','error')} for r in rows if r['status']!='analysed'],
         'moco_biomass_models':[r['accession'] for r in analysed if r['moco_biomass']],
         'organism_count':len({r['organism'] for r in analysed}),
         'caveat':'Distinct BiGG model accession IDs are model datasets, not independent wet-lab observations. Objective matching uses annotated metabolite name or ID. Model growth was not scored; media and reaction conventions differ by model.'}
(OUT/'bigg_model_survey.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
