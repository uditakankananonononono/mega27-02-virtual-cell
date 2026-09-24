"""Pre-registered multi-organism replication with CarveMe models (notes/prereg_cross_species.md, amendments 3-4)."""
import json, logging, os, warnings
import numpy as np, pandas as pd, cobra
from scipy.stats import mannwhitneyu
from vcell.rescue import rescue_audit, load_kegg_links
warnings.filterwarnings('ignore'); logging.disable(logging.CRITICAL)
D = 'data/cross/'
ORGS = {'Caulo': 'ccs', 'Cola': 'evi', 'PS': 'dsu', 'Ponti': 'pact', 'Smeli': 'sme', 'SyringaeB728a': 'psb'}
M9 = ['ca2', 'cl', 'cobalt2', 'cu2', 'fe2', 'fe3', 'glc__D', 'h2o', 'h', 'k', 'mg2', 'mn2', 'mobd', 'na1', 'nh4', 'ni2', 'o2', 'pi', 'so4', 'zn2']
SUP = {'btn_c': {'00780'}, 'thm_c': {'00730'}, 'pnto__R_c': {'00770'}, 'thf_c': {'00790', '00670'}, 'nad_c': {'00760'},
       'amet_c': {'00270'}, 'pydx5p_c': {'00750'}}
per, pooled = {}, []
for org, kc in ORGS.items():
    path = D + f'carve/{org}.xml'
    if not os.path.exists(path):
        per[org] = {'excluded': 'CarveMe model not built'}; continue
    m = cobra.io.read_sbml_model(path)
    med = {f'EX_{c}_e': (10.0 if c == 'glc__D' else 1000.0) for c in M9 if f'EX_{c}_e' in m.reactions}
    m.medium = med
    wt = m.slim_optimize(error_value=0.0)
    if wt < 1e-6:
        per[org] = {'excluded': 'criterion (b): no growth on M9 glucose', 'wt': wt}; continue
    kg = load_kegg_links(D + f'kegg_{kc}_pathway.tsv')
    rows = rescue_audit(m, SUP, {g.id: kg.get(g.id, set()) for g in m.genes})
    pd.DataFrame(rows).to_csv(f'results/rescue_audit_carveme_{org}.csv', index=False)
    lab, sups = {}, {}
    for r in rows:
        g = r['gene']; lab[g] = 'on_pathway' if (lab.get(g) == 'on_pathway' or r['label'] == 'on_pathway') else 'off_pathway'
        sups.setdefault(g, set()).add(r['supplement'])
    fit = pd.read_csv(D + f'{org}_gene_median_fitness.tsv', sep='\t', dtype={'locusId': str}).dropna(subset=['median_fitness'])
    pct = dict(zip(fit.locusId, fit.median_fitness.rank(pct=True)))
    medf = dict(zip(fit.locusId, fit.median_fitness))
    on = [medf[g] for g in lab if lab[g] == 'on_pathway' and g in medf]; off = [medf[g] for g in lab if lab[g] == 'off_pathway' and g in medf]
    rec = {'wt': wt, 'n_genes_model': len(m.genes), 'missing_supplements': [k for k in SUP if k not in m.metabolites],
           'n_rescue_pairs': len(rows), 'genes_on': sum(v == 'on_pathway' for v in lab.values()), 'genes_off': sum(v == 'off_pathway' for v in lab.values()),
           'with_fitness_on': len(on), 'with_fitness_off': len(off),
           'genes': {g: {'label': lab[g], 'supplements': sorted(sups[g]), 'median_fitness': medf.get(g), 'percentile': pct.get(g)} for g in lab}}
    if len(on) < 5 or len(off) < 5: rec['verdict'] = 'UNDERPOWERED'
    else:
        d = float(np.median(on) - np.median(off)); p = float(mannwhitneyu(on, off, alternative='greater').pvalue)
        rec.update(median_on=float(np.median(on)), median_off=float(np.median(off)), median_diff=d, mwu_p=p,
                   verdict='REPLICATED' if (p < 0.05 and d > 0) else ('FALSIFIED' if d <= 0 else 'NOT SIGNIFICANT (direction consistent)'))
    per[org] = rec
    for g, v in rec['genes'].items():
        if v['percentile'] is not None: pooled.append((org, g, v['label'], v['percentile'], v['supplements']))
def ptest(items, min_n=10):
    a = [x[3] for x in items if x[2] == 'on_pathway']; b = [x[3] for x in items if x[2] == 'off_pathway']
    out = {'n_on': len(a), 'n_off': len(b)}
    if len(a) < min_n or len(b) < min_n: out['verdict'] = 'UNDERPOWERED'; return out
    d = float(np.median(a) - np.median(b)); p = float(mannwhitneyu(a, b, alternative='greater').pvalue)
    out.update(median_pct_on=float(np.median(a)), median_pct_off=float(np.median(b)), diff=d, mwu_p=p,
               verdict='REPLICATED' if (p < 0.05 and d > 0) else ('FALSIFIED' if d <= 0 else 'NOT SIGNIFICANT (direction consistent)'))
    return out
res = {'per_organism': per, 'primary_pooled_percentile': ptest(pooled),
       'secondary_pooled_excluding_SAM': ptest([x for x in pooled if x[4] != ['amet_c']])}
json.dump(res, open('results/cross_species_carveme.json', 'w'), indent=1)
print(json.dumps({'primary': res['primary_pooled_percentile'], 'noSAM': res['secondary_pooled_excluding_SAM'],
                  'per': {k: {kk: vv for kk, vv in v.items() if kk != 'genes'} for k, v in per.items()}}, indent=1))
