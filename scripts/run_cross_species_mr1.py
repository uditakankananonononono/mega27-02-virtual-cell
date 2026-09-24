"""Pre-registered cross-species replication (notes/prereg_cross_species.md): S. oneidensis MR-1."""
import json, warnings
import numpy as np, pandas as pd, cobra
from scipy.stats import mannwhitneyu
from vcell.rescue import rescue_audit, load_kegg_links
warnings.filterwarnings('ignore')
D = 'data/cross/'
m = cobra.io.read_sbml_model(D + 'iMR1_799.xml')
BASE = ['EX_O2_E', 'EX_NH4_E', 'EX_PI_E', 'EX_SO4_E', 'EX_H2O_E', 'EX_H_E', 'EX_CO2_E', 'EX_FE2_E', 'EX_FE3_E', 'EX_K_E',
        'EX_MG2_E', 'EX_NA1_E', 'EX_CL_E', 'EX_CA2_E', 'EX_MN2_E', 'EX_CU2_E', 'EX_COBALT2_E', 'EX_MOBD_E', 'EX_NI2_E']
for r in m.reactions:
    if not r.id.startswith('EX_'): continue
    if r not in m.boundary: r.bounds = (0, 0)   # multi-metabolite pseudo-exchanges (casamino acids, gelatin, Tween 20)
    else: r.lower_bound = -(10 if r.id == 'EX_LAC_L_E' else 1000 if r.id in BASE else 0)
wt = m.slim_optimize(error_value=0.0)
SUP = {'btn': {'00780'}, 'thm': {'00730'}, 'pnto_r': {'00770'}, 'thf': {'00790', '00670'}, 'nad': {'00760'},
       'amet': {'00270'}, 'pydx5p': {'00750'}}
norm = lambda g: g.replace('g_', '').replace('SO_', 'SO')
kg = {norm(k): v for k, v in load_kegg_links(D + 'kegg_son_pathway.tsv').items()}
gp = {g.id: kg.get(norm(g.id), set()) for g in m.genes}
rows = rescue_audit(m, SUP, gp)
pd.DataFrame(rows).to_csv('results/rescue_audit_MR1_lactate.csv', index=False)
lab = {}
for r in rows:
    g = norm(r['gene']); lab[g] = 'on_pathway' if (lab.get(g) == 'on_pathway' or r['label'] == 'on_pathway') else 'off_pathway'
fit = pd.read_csv(D + 'MR1_fit_logratios_good.tab', sep='\t')
expcols = [c for c in fit.columns if c not in ('locusId', 'sysName', 'desc', 'comb')]
med = dict(zip(fit.sysName, fit[expcols].median(axis=1)))
on = [med[g] for g in lab if lab[g] == 'on_pathway' and g in med]
off = [med[g] for g in lab if lab[g] == 'off_pathway' and g in med]
n_on, n_off = sum(v == 'on_pathway' for v in lab.values()), sum(v == 'off_pathway' for v in lab.values())
out = {'wt_growth_lactate': wt, 'n_rescue_pairs': len(rows), 'n_experiments': len(expcols),
       'genes_on': n_on, 'genes_off': n_off, 'with_fitness_on': len(on), 'with_fitness_off': len(off),
       'absent_frac_on': 1 - len(on) / n_on if n_on else None, 'absent_frac_off': 1 - len(off) / n_off if n_off else None,
       'labels': lab, 'gene_median_fitness': {g: med.get(g) for g in lab}}
if len(on) < 5 or len(off) < 5:
    out['verdict'] = 'UNDERPOWERED (fewer than 5 genes with fitness in a group)'
else:
    d = float(np.median(on) - np.median(off)); p = float(mannwhitneyu(on, off, alternative='greater').pvalue)
    out.update(median_on=float(np.median(on)), median_off=float(np.median(off)), median_diff=d, mwu_p_one_sided=p)
    out['verdict'] = 'REPLICATED' if (p < 0.05 and d > 0) else ('FALSIFIED in MR-1' if d <= 0 else 'NOT SIGNIFICANT (direction consistent)')
json.dump(out, open('results/cross_species_mr1.json', 'w'), indent=1)
print(json.dumps({k: v for k, v in out.items() if k not in ('labels', 'gene_median_fitness')}, indent=1))
