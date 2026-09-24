"""Per-cofactor rescue attribution on the Bernstein 2023 benchmark (their pipeline, their metric).

For a base model variant (default: Bernstein all-corrections iML1515), every model no-growth
(gene, carbon) pair is re-simulated with each candidate intracellular supplement opened alone.
Output: results/bernstein/attrib_<tag>.npz with the rescue matrix (pairs x cofactors).
Supplements can only add growth, so only no-growth pairs need re-simulation.

usage: bernstein_cofactor_attribution.py MODEL_XML_OR_NAME BASE_TAG ATTRIB_TAG
"""
import sys, json, os, time
import numpy as np
BASE = os.path.abspath('external/E_coli_GEM_validation') + '/'
nb = json.load(open(BASE + 'Analysis/Analysis_Notebook.ipynb'))
src = {i: ''.join(c['source']) for i, c in enumerate(nb['cells']) if c['cell_type'] == 'code'}
import cobra, pandas as pd, copy
from sklearn.metrics import precision_recall_curve as pre_rec, auc as sk_auc, roc_auc_score
g = {'cobra': cobra, 'pd': pd, 'np': np, 'copy': copy, 'pre_rec': pre_rec, 'sk_auc': sk_auc,
     'roc_auc_score': roc_auc_score}
for i in (9, 11, 13, 15, 17, 19, 21, 23):
    exec(src[i], g)
model_arg, base_tag, tag = sys.argv[1], sys.argv[2], sys.argv[3]
COF = json.load(open('results/auto_cofactor_set.json'))
t0 = time.time()
if model_arg.endswith('.xml'):
    model = cobra.io.read_sbml_model(model_arg)
    for ex in model.exchanges:
        ex.lower_bound = 0; ex.upper_bound = 1000
else:
    model = g['load_model'](model_arg, BASE)
med, carb, carb_exp = g['load_environment'](BASE)
dexp, dgenes, dfit = g['load_data'](BASE)
gm, cem, cmm, fm = g['match_model_data'](model, carb, carb_exp, dexp, dgenes, dfit)
g['name_genes_matched'] = gm
m, gma, cema, cmma, fma = g['model_adjustments'](1, 1, 1, model, gm, cem, cmm, fm)
mei, cei = g['check_environment'](m, med, cmma)
ids = json.load(open(f'results/bernstein/{base_tag}_ids.json'))
assert list(gma) == ids['genes'] and list(cmma) == ids['carbon'], 'gene/carbon order mismatch vs base run'
sim0 = np.load(f'results/bernstein/{base_tag}_sim.npy')
med_rx = [m.exchanges[i] for i in mei if i != -1]
carb_rx = [m.exchanges[i] if i != -1 else None for i in cei]
for r in med_rx:
    r.lower_bound = -1000
sup = {}
for c in COF:
    if c + '_c' in m.metabolites:
        r = cobra.Reaction('SUP_' + c); r.lower_bound = 0; r.upper_bound = 0
        r.add_metabolites({m.metabolites.get_by_id(c + '_c'): -1.0}); m.add_reactions([r]); sup[c] = r
cofs = list(sup)
ng = np.argwhere(sim0 <= 0.001)
R = np.zeros((len(ng), len(cofs)), dtype=np.int8)
allrescue = np.zeros(len(ng), dtype=np.int8)
for k, (gi, ci) in enumerate(ng):
    cr = carb_rx[ci]
    if cr is None:
        continue
    with m:
        cr.lower_bound = -10
        m.genes.get_by_id(gma[gi]).knock_out()
        for r in sup.values():
            r.lower_bound = -1000
        v = m.slim_optimize()
        if np.isnan(v) or v <= 0.001:
            continue
        allrescue[k] = 1
        for r in sup.values():
            r.lower_bound = 0
        for j, c in enumerate(cofs):
            sup[c].lower_bound = -1000
            v = m.slim_optimize()
            R[k, j] = int(not np.isnan(v) and v > 0.001)
            sup[c].lower_bound = 0
np.savez(f'results/bernstein/attrib_{tag}.npz', pairs=ng, R=R, allrescue=allrescue, cofs=np.array(cofs))
print(json.dumps({'tag': tag, 'n_nogrowth': int(len(ng)), 'n_rescued_by_all': int(allrescue.sum()),
                  'per_cof': dict(zip(cofs, R.sum(0).tolist())), 'seconds': round(time.time() - t0, 1)}, indent=1))
