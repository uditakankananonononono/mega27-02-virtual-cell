"""Pre-registered split-half greedy cofactor selection (notes/prereg_cofactor_split.md).
Uses single-cofactor rescue attribution (approximation); the selected set is then verified with a full run."""
import json, numpy as np
from sklearn.metrics import precision_recall_curve as pre_rec, auc as sk_auc
B = 'results/bernstein/'
import sys
base = sys.argv[1] if len(sys.argv) > 1 else 'iML1515_bernstein_allcorr'
atag = sys.argv[2] if len(sys.argv) > 2 else 'allcorr'
OUT = sys.argv[3] if len(sys.argv) > 3 else 'results/bernstein_split_select.json'
sim = np.load(B + base + '_sim.npy'); fit = np.load(B + base + '_fit.npy'); ids = json.load(open(B + base + '_ids.json'))
A = np.load(B + 'attrib_' + atag + '.npz'); pairs, R, cofs = A['pairs'], A['R'], list(A['cofs'])
order = sorted(range(len(ids['carbon'])), key=lambda i: (ids['carbon'][i], i))
sel_c = set(order[0::2]); test_c = set(order[1::2])
def prauc(bmat, cols):
    cols = sorted(cols); b = bmat[:, cols].flatten(); f = fit[:, cols].flatten()
    p, r, _ = pre_rec(b, -f, pos_label=0); return float(sk_auc(r, p))
b0 = (sim > 0.001).astype(int)
def apply(S):
    b = b0.copy()
    if S:
        idx = [cofs.index(c) for c in S]; hit = R[:, idx].max(1) > 0
        for (gi, ci) in pairs[hit]: b[gi, ci] = 1
    return b
S = []; cur = prauc(b0, sel_c); trace = [{'set': [], 'sel_prauc': cur}]
while True:
    best = None
    for c in cofs:
        if c in S: continue
        v = prauc(apply(S + [c]), sel_c)
        if v > cur + 1e-9 and (best is None or v > best[1]): best = (c, v)
    if best is None: break
    S.append(best[0]); cur = best[1]; trace.append({'set': list(S), 'sel_prauc': cur})
out = {'selected': S, 'trace': trace, 'sel_carbons': [ids['carbon'][i] for i in sorted(sel_c)],
       'test_carbons': [ids['carbon'][i] for i in sorted(test_c)],
       'approx_test_prauc_selected': prauc(apply(S), test_c), 'test_prauc_allcorr': prauc(b0, test_c),
       'approx_all_prauc_selected': prauc(apply(S), range(fit.shape[1])), 'all_prauc_allcorr': prauc(b0, range(fit.shape[1])),
       'per_cof_rescues': dict(zip(cofs, R.sum(0).tolist()))}
json.dump(out, open(OUT, 'w'), indent=1); print(json.dumps({k: v for k, v in out.items() if k not in ('trace',)}, indent=1))
