import json
from pathlib import Path

def test_local_model_inclusive_inventory():
    root=Path(__file__).resolve().parents[1]
    x=json.loads((root/'results/model_inclusive_inventory_local_verification.json').read_text())
    assert x['biGG_unique_model_ids']==x['biGG_hash_matched_local_files']==108
    assert x['paxdb_unique_dataset_ids']==x['paxdb_hash_matched_local_files']==x['paxdb_saved_scores_matched']==19
    assert x['paxdb_distinct_provider_publication_links']==16
    assert x['technical_variant_link_multiplicities']==[2,3]
    assert x['conservative_model_inclusive_arithmetic']==123
    assert x['status'].startswith('LOCAL_LEDGER_CHECK_PASSED')
