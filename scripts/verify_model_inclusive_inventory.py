"""Read-only local audit of the 108+15 model-inclusive inventory.

Not an online provider-currentness or biological-independence assertion.
"""
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

R=Path(__file__).resolve().parents[1]
models=list(csv.DictReader((R/'results/bigg_model_survey.csv').open()))
pax=list(csv.DictReader((R/'results/paxdb_accession_evidence.csv').open()))
score=json.loads((R/'results/paxdb_datasets.json').read_text())['per_dataset']
assert len(models)==len({x['accession'] for x in models})==108
assert all(x['status']=='analysed' and x['source_url']=='https://bigg.ucsd.edu/static/models/'+x['accession']+'.json' for x in models)
assert all(hashlib.sha256((R/x['path']).read_bytes()).hexdigest()==x['sha256'] for x in models)
assert len(pax)==len({x['paxdb_id'] for x in pax})==19
assert all(hashlib.sha256((R/x['file']).read_bytes()).hexdigest()==x['sha256'] for x in pax)
assert all(x['source_url'].startswith('https://pax-db.org/downloads/latest/datasets/511145/') for x in pax)
assert all(float(x['scored_auroc'])==score[Path(x['file']).name]['auroc'] for x in pax)
# Provider #link metadata is a grouping lead, not a guarantee of independent studies.
links=Counter(x['publication'] for x in pax)
assert sorted(n for n in links.values() if n>1)==[2,3]
assert len(links)==16
# The older conservative figure 15/19 is retained as a lower bound: the
# 16 source-link strings could hide an additional shared study/source.
out={'status':'LOCAL_LEDGER_CHECK_PASSED; owner-approved counting unit, not a global completion gate',
     'biGG_unique_model_ids':len(models),'biGG_hash_matched_local_files':len(models),
     'biGG_status_analysed':len(models),'paxdb_unique_dataset_ids':len(pax),
     'paxdb_hash_matched_local_files':len(pax),'paxdb_saved_scores_matched':len(pax),
     'paxdb_distinct_provider_publication_links':len(links),
     'technical_variant_link_multiplicities':sorted(n for n in links.values() if n>1),
     'conservative_paxdb_study_floor_carried_forward':15,
     'conservative_model_inclusive_arithmetic':108+15,
     'caveats':'Local files/hashes and saved numerical use verified; provider freshness not rechecked; no proof that 15 links represent independent wet-lab studies or model snapshots independent organisms; eight other provisional manifest rows remain; 50+ substantive body pages and original G2 still fail.'}
(R/'results/model_inclusive_inventory_local_verification.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
