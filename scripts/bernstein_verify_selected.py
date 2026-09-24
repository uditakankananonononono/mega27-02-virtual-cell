"""Pre-registered test (notes/prereg_cofactor_split.md): full-run selected set vs allcorr on held-out carbons."""
import json, sys, numpy as np; sys.path.insert(0, '.')
from vcell.bernstein import _pr, paired_gene_bootstrap
B = 'results/bernstein/'
a, s = 'iML1515_bernstein_allcorr', 'iML1515_allcorr_selected'
ia, is_ = json.load(open(B + a + '_ids.json')), json.load(open(B + s + '_ids.json'))
assert ia == is_, 'gene/carbon order differs'
sel = json.load(open('results/bernstein_split_select.json'))
test = [ia['carbon'].index(c) for c in sel['test_carbons']]
fa, fs = np.load(B + a + '_fit.npy'), np.load(B + s + '_fit.npy'); assert np.allclose(fa, fs)
ba = (np.load(B + a + '_sim.npy') > 0.001).astype(int); bs = (np.load(B + s + '_sim.npy') > 0.001).astype(int)
allc = list(range(fa.shape[1]))
out = {'selected': sel['selected'], 'test_carbons': sel['test_carbons'],
       'heldout_prauc_selected': _pr(bs[:, test], fa[:, test]), 'heldout_prauc_allcorr': _pr(ba[:, test], fa[:, test]),
       'all25_prauc_selected': _pr(bs, fa), 'all25_prauc_allcorr': _pr(ba, fa),
       'bootstrap_heldout': paired_gene_bootstrap(bs, ba, fa, test, n=2000, seed=0),
       'bootstrap_all25': paired_gene_bootstrap(bs, ba, fa, allc, n=2000, seed=0),
       'n_calls_flipped_to_growth': int(((bs == 1) & (ba == 0)).sum()), 'n_calls_flipped_to_nogrowth': int(((bs == 0) & (ba == 1)).sum())}
lo = out['bootstrap_heldout']['ci95'][0]
out['verdict'] = ('BREAK: held-out PR-AUC above published all-corrections model, paired gene bootstrap 95% CI excludes 0'
                  if out['heldout_prauc_selected'] > out['heldout_prauc_allcorr'] and lo > 0 else 'NEGATIVE: pre-registered criterion not met')
json.dump(out, open('results/bernstein_selected_verdict.json', 'w'), indent=1); print(json.dumps(out, indent=1))
