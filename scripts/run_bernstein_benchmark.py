"""Head-to-head on the Bernstein et al. 2023 (Mol Syst Biol) benchmark, using THEIR released code
(external/E_coli_GEM_validation, MIT-licensed notebook) for loading, matching, strain/essential/carbon
adjustments, simulation and the precision-recall AUC metric. Only the model variant changes.

usage: run_bernstein_benchmark.py MODEL_NAME_OR_XML VARIANT_TAG [intracellular-supplement ids ...]
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

model_arg, tag, supplements = sys.argv[1], sys.argv[2], sys.argv[3:]
t0 = time.time()
if model_arg.endswith('.xml'):
    model = cobra.io.read_sbml_model(model_arg)
    # keep exchanges the saved model supplies (Bernstein saved allcorr with EX_btn_e/EX_thm_e/EX_pnto__R_e open);
    # closing them (bug fixed 10:24 PM) silently removed their vitamin correction
    for ex in model.exchanges:
        ex.lower_bound = min(ex.lower_bound, 0) if ex.lower_bound < 0 else 0; ex.upper_bound = 1000
else:
    model = g['load_model'](model_arg, BASE)
med, carb, carb_exp = g['load_environment'](BASE)
dexp, dgenes, dfit = g['load_data'](BASE)
gm, cem, cmm, fm = g['match_model_data'](model, carb, carb_exp, dexp, dgenes, dfit)
# their model_adjustments body reads the notebook-global name_genes_matched (param is names_genes_matched)
g['name_genes_matched'] = gm
model_adj, gma, cema, cmma, fma = g['model_adjustments'](1, 1, 1, model, gm, cem, cmm, fm)
added = []
for mid in supplements:  # intracellular exchange, same construction as Bernstein Part 6
    if mid + '_c' in model_adj.metabolites:
        r = cobra.Reaction('EX_' + mid + '_c'); r.lower_bound = -1000; r.upper_bound = 1000
        r.add_metabolites({model_adj.metabolites.get_by_id(mid + '_c'): -1.0})
        model_adj.add_reactions([r]); added.append(mid)
mei, cei = g['check_environment'](model_adj, med, cmma)
sim = g['simulate_phenotype'](model_adj, gma, cmma, mei, cei)
# THEIR metric (plot_precision_recall_curve, cell 29): labels = model growth call, score = -fitness, pos_label=0
b = (sim > 0.001).astype(int).flatten(); f = fma.flatten()
pre, rec, _ = pre_rec(b, f * -1, pos_label=0)
out = {'model': model_arg, 'variant': tag, 'supplements_added': added, 'n_genes': len(gma),
       'n_carbon': len(cmma), 'n_pairs': int(b.size), 'n_model_nogrowth': int((b == 0).sum()),
       'pr_auc_bernstein_metric': float(sk_auc(rec, pre)), 'roc_auc': float(roc_auc_score(b, f)),
       'seconds': round(time.time() - t0, 1)}
os.makedirs('results/bernstein', exist_ok=True)
np.save(f'results/bernstein/{tag}_sim.npy', sim); np.save(f'results/bernstein/{tag}_fit.npy', fma)
json.dump({'genes': list(gma), 'carbon': list(cmma)}, open(f'results/bernstein/{tag}_ids.json', 'w'))
json.dump(out, open(f'results/bernstein/{tag}.json', 'w'), indent=2)
print(json.dumps(out, indent=2))
