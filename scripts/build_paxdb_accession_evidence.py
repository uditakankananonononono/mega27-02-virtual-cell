"""Hash E. coli PaxDb accession files actually scored; distinguish studies from file IDs."""
import csv,hashlib,json,re
from pathlib import Path
R=Path(__file__).resolve().parents[1]
per=json.loads((R/'results/paxdb_datasets.json').read_text())['per_dataset']
rows=[]
for p in sorted((R/'data/external/paxdb511145').glob('511145-*.txt')):
    text=p.read_text(errors='replace')
    pid=re.search(r'^#id:\s*(\d+)',text,re.M)
    if not pid:continue
    name=p.name;d=per.get(name,{})
    pub=re.search(r'^#link:\s*(.*)',text,re.M)
    rows.append(dict(paxdb_id=pid.group(1),source_url='https://pax-db.org/downloads/latest/datasets/511145/'+name,
                     file=str(p.relative_to(R)),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),bytes=p.stat().st_size,
                     scored_auroc=d.get('auroc',''),n_genes_scored=d.get('n_covered',''),publication=pub.group(1) if pub else ''))
assert len(rows)==19 and len({r['paxdb_id'] for r in rows})==19
with (R/'results/paxdb_accession_evidence.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print('PaxDb accessions, hashed and scored',len(rows),'all scored',all(r['scored_auroc']!='' for r in rows))
