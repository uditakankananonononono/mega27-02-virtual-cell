"""Read-only source identity checks for two already analysed public tables."""
import hashlib,json,urllib.request
from pathlib import Path
checks=[('KEGG eco pathway','https://rest.kegg.jp/link/pathway/eco','data/external/kegg_eco_pathway.tsv'),('TCDB substrates','https://tcdb.org/cgi-bin/substrates/getSubstrates.py','data/external/tcdb/tcdb_substrates.tsv')]
out=[]
for name,url,path in checks:
 with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'VC2-research/1.0'}),timeout=45) as r:raw=r.read()
 saved=Path(path).read_bytes();live_rows=raw.splitlines();saved_rows=saved.splitlines()
 result={'source':name,'url':url,'saved_path':path,'live_sha256':hashlib.sha256(raw).hexdigest(),'saved_sha256':hashlib.sha256(saved).hexdigest(),'live_rows':len(live_rows),'saved_rows':len(saved_rows),'exact_bytes':raw==saved,'same_rows_ignoring_order':sorted(live_rows)==sorted(saved_rows),'caveat':'Current public source identity check; historical acquisition time cannot be independently reconstructed.'}
 assert result['same_rows_ignoring_order'],name
 out.append(result)
Path('results/tool_source_snapshot_verification.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
