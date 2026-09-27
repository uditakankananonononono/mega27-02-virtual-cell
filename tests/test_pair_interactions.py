import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]/'results'
def test_saved_pair_interactions_and_fresh_check():
    x=json.loads((ROOT/'postresult_pair_interactions.json').read_text())
    f=json.loads((ROOT/'postresult_pair_interactions_fresh_check.json').read_text())
    assert x['model_sha256']==f['source_sha256']
    assert x['n_pairs']==253 and x['n_exact_single_union']==250
    assert x['n_pair_only']==2 and x['n_lost_single_union']==1
    special={tuple(r['pair']):r for r in x['rows'] if r['pair_only_genes'] or r['single_union_only_genes']}
    assert special[('hemeO_c','pheme_c')]['pair_only_genes']==['PP_0189','PP_0744','PP_5074']
    assert special[('adocbl_c','fad_c')]['pair_only_genes']==['PP_0602']
    assert len(special[('2fe2s_c','btamp_c')]['single_union_only_genes'])==8
    for row in f['pairs']:
        assert all(g['status']=='optimal' for variant in row['variants'] for g in variant['gene_knockouts'].values())
    antagonistic=next(r for r in f['pairs'] if r['pair']==['2fe2s_c','btamp_c'])
    assert antagonistic['variants'][3]['wt_growth']<antagonistic['variants'][2]['wt_growth']
