import json
from pathlib import Path
from collections import Counter

ROOT = Path(__file__).resolve().parents[1]

def test_complete_equal_coefficient_pair_null():
    x=json.loads((ROOT/'results/postresult_equal_coefficient_pair_null.json').read_text())
    manifest=json.loads((ROOT/'results/biomass_topology_manifest.json').read_text())
    previous=json.loads((ROOT/'results/moco_subset_interaction_ijn1463.json').read_text())
    mo=next(z for z in previous['subsets'] if z['removed']==['bmocogdp_c','mocogdp_c'])
    assert x['source_sha256']==manifest['sha256']==previous['source_sha256']
    assert x['n_original_zero_growth_genes']==previous['baseline_essential_count']==262
    assert x['n_pairs']==len(x['pairs'])==253
    assert len(x['eligible_non_moco_terms'])==23
    assert x['mo_pair_fresh_checksum']['matches_prior']
    assert sorted(x['mo_pair_fresh_checksum']['rescued_ko_growth'])==sorted(mo['rescued_ko_growth'])
    assert abs(x['mo_pair_fresh_checksum']['edited_wt_growth']-mo['edited_wt_growth'])<1e-8
    counts=Counter(z['n_rescued'] for z in x['pairs'])
    assert sum(v for k,v in counts.items() if k>=6)==x['summary']['n_pairs_ge_mo_six']==104
    assert sum(v for k,v in counts.items() if k>6)==x['summary']['n_pairs_gt_mo_six']==88
    assert counts[0]==x['summary']['zero_rescue_pairs']==78
    assert max(counts)==x['summary']['max_non_moco_rescues']==35
    assert all(not z['solver_failure_genes'] for z in x['pairs'])
    assert all(not z['borderline_genes'] for z in x['pairs'])
    assert all(not z['overlap_with_mo_pair'] for z in x['pairs'])
    assert all(v==0 for v in x['summary']['mo_gene_overlap_counts'].values())
    assert len({tuple(z['removed']) for z in x['pairs']})==253
    assert all(len(z['rescued_ko_growth'])==z['n_rescued'] for z in x['pairs'])
