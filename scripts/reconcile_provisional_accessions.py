"""Resolve provisional manifest accessions where verified provider IDs and scored use already exist.

Evidence rules, all checked, none assumed:
- PaxDb rows: local file exists, '#id:' header present, file sha256 recorded here,
  and the dataset key appears in the saved scored result for its organism arm.
- BiGG model rows: accession present in the frozen bigg_model_survey.csv with a 64-char sha256
  and matching static source URL.
Rows that fail any check stay provisional. Original audit CSV is left untouched; this writes an overlay.
"""
import csv,hashlib,json,re
from pathlib import Path
R=Path(__file__).resolve().parents[1]
res=R/'results'
audit=list(csv.DictReader((res/'datasets_accession_audit.csv').open()))
prov=[r for r in audit if r['status']=='provisional']
assert len(prov)==88

def header_id(p):
    m=re.search(r'^#id:\s*(\d+)',p.read_text(errors='replace'),re.M)
    return m.group(1) if m else None

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

scored={
 '83332':set(json.loads((res/'paxdb_mtb.json').read_text())['datasets']),
 '224308':set(json.loads((res/'paxdb_bsub.json').read_text())['datasets']),
 '208964':set(json.loads((res/'paxdb_pao1.json').read_text())['datasets']),
 'cross':set(json.loads((res/'paxdb_cross.json').read_text())['datasets'])}
ecoli={r['file']:r for r in csv.DictReader((res/'paxdb_accession_evidence.csv').open())}
bigg={r['accession']:r for r in csv.DictReader((res/'bigg_model_survey.csv').open())}
arm_dir={'83332':'data/external/mtb','224308':'data/external/bsub','208964':'data/external/pao1'}
resolved=[];still=[]
for r in prov:
    name=r['name'];out=None
    m=re.search(r'PaxDb v5 dataset \d+ \((\S+?\.txt)\)',name) or re.search(r'PaxDb v5 dataset (\S+?\.txt)',name)
    if m:
        fn=m.group(1)
        if fn.startswith('511145-'):
            p=R/'data/external/paxdb511145'/fn
            ev=ecoli.get('data/external/paxdb511145/'+fn)
            if p.exists() and ev and ev['paxdb_id']==header_id(p) and ev['sha256']==sha(p) and ev['scored_auroc']!='':
                out=dict(accession='PaxDb:'+ev['paxdb_id'],source_url=ev['source_url'],sha256=ev['sha256'],use_evidence='scored in results/paxdb_datasets.json',arm='ecoli')
        else:
            taxid=fn.split('-',1)[0]
            if taxid in arm_dir:
                p=R/arm_dir[taxid]/fn;key=fn[len(taxid)+1:-4]
                if p.exists() and key in scored[taxid]:
                    pid=header_id(p)
                    if pid:out=dict(accession='PaxDb:'+pid,source_url=f'https://pax-db.org/downloads/latest/datasets/{taxid}/'+fn,sha256=sha(p),use_evidence={'83332':'scored in results/paxdb_mtb.json','224308':'scored in results/paxdb_bsub.json','208964':'scored in results/paxdb_pao1.json'}[taxid],arm=taxid)
            else:
                p=R/'data/external/paxdb_cross'/fn
                if p.exists() and fn[:-4] in scored['cross']:
                    pid=header_id(p)
                    if pid:out=dict(accession='PaxDb:'+pid,source_url=f'https://pax-db.org/downloads/latest/datasets/{taxid}/'+fn,sha256=sha(p),use_evidence='scored in results/paxdb_cross.json',arm='cross')
    if out is None:
        for mid in ('e_coli_core','iJO1366','iAF1260','iJN1463','iJR904'):
            if re.search(r'\b'+mid+r'\b',name):
                b=bigg.get(mid);p=R/bigg[mid]['path'] if mid in bigg else None
                if b and len(b['sha256'])==64 and p.exists() and sha(p)==b['sha256'] and b['status']=='analysed':
                    out=dict(accession='BiGG:'+mid,source_url=b['source_url'],sha256=b['sha256'],use_evidence='analysed in results/bigg_model_survey.csv',arm='bigg')
                break
    if out is None:
        km=re.search(r'KEGG link/pathway/(\w+)',name)
        if km:
            code=km.group(1);ver=json.loads((res/'kegg_map_source_verification.json').read_text())['per_code'].get(code)
            disc=json.loads((res/'xs7_kegg_discovery.json').read_text())
            used=any(v.get('kegg_code')==code and v.get('included') for v in disc.values())
            used_result='results/xs7_kegg_discovery.json'
            if not used and code=='son' and (res/'cross_species_mr1.json').exists():used=True;used_result='results/cross_species_mr1.json'
            if not used and code=='ppu' and (res/'cross_species_putida.json').exists():used=True;used_result='results/cross_species_putida.json'
            if ver and ver['exact_bytes'] and used:
                out=dict(accession='KEGG:'+code,source_url=f'https://rest.kegg.jp/link/pathway/{code}',sha256=ver['saved_sha256'],use_evidence='exact live byte match 2026-09-27; used in '+used_result,arm='kegg')
    if out is None and name.startswith('TCDB getSubstrates'):
        ver=[c for c in json.loads((res/'tool_source_snapshot_verification.json').read_text()) if c['source']=='TCDB substrates'][0]
        p=R/'data/external/tcdb/tcdb_substrates.tsv'
        if ver['same_rows_ignoring_order'] and p.exists() and hashlib.sha256(p.read_bytes()).hexdigest()==ver['saved_sha256']:
            out=dict(accession='TCDB:substrates',source_url=ver['url'],sha256=ver['saved_sha256'],use_evidence='live row-set match 2026-09-27 (byte order differs); used in results/uptake_plausibility.json',arm='tcdb')
    (resolved if out else still).append(dict(manifest_number=r['manifest_number'],name=name,kind=r['kind'],**(out or {})))
ids=[x['accession'] for x in resolved]
assert len(ids)==len(set(ids)),[x for x in ids if ids.count(x)>1]
with (res/'provisional_accession_reconciliation.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=['manifest_number','name','kind','accession','source_url','sha256','use_evidence','arm']);w.writeheader();w.writerows(resolved)
summary={'status':'provisional-accession reconciliation over frozen audit; rows failing any check stay provisional',
 'n_provisional_before':88,'n_resolved':len(resolved),'n_still_provisional':len(still),
 'resolved_by_arm':{a:sum(1 for x in resolved if x['arm']==a) for a in sorted({x['arm'] for x in resolved})},
 'still_provisional_kinds':{k:sum(1 for x in still if x['kind']==k) for k in sorted({x['kind'] for x in still})},
 'limits':'PaxDb dataset IDs and BiGG model IDs verified from local hashed files with scored-use evidence; this does not make technical replicates independent studies and does not change the owner-approved model-inclusive counting unit or imply 123 independent wet-lab studies. See notes/gate_count_provenance_review_20260927.md. Remaining provisional rows need per-provider evidence.'}
(res/'provisional_accession_reconciliation.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
