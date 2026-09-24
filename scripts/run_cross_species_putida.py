"""Pre-registered cross-species replication, amendment 1 (notes/prereg_cross_species.md): P. putida KT2440."""
import gzip, io, json, warnings
import numpy as np, pandas as pd, cobra
from scipy.stats import mannwhitneyu
from vcell.rescue import rescue_audit, load_kegg_links
warnings.filterwarnings('ignore')
D = 'data/cross/'
m = cobra.io.load_json_model(io.StringIO(gzip.open(D + 'iJN1463.json.gz', 'rt').read()))
INORG = {'EX_ca2_e', 'EX_cl_e', 'EX_co2_e', 'EX_cobalt2_e', 'EX_cu2_e', 'EX_fe2_e', 'EX_h2o_e', 'EX_h_e', 'EX_hco3_e', 'EX_k_e',
         'EX_mg2_e', 'EX_mn2_e', 'EX_mobd_e', 'EX_na1_e', 'EX_nh4_e', 'EX_ni2_e', 'EX_o2_e', 'EX_pi_e', 'EX_sel_e', 'EX_so4_e',
         'EX_tungs_e', 'EX_zn2_e'}
med = {r: 1000.0 for r in INORG if r in m.reactions}; med['EX_glc__D_e'] = 10.0
m.medium = med
wt = m.slim_optimize(error_value=0.0)
SUP = {'btn_c': {'00780'}, 'thm_c': {'00730'}, 'pnto__R_c': {'00770'}, 'thf_c': {'00790', '00670'}, 'nad_c': {'00760'},
       'amet_c': {'00270'}, 'pydx5p_c': {'00750'}}
missing_sup = [k for k in SUP if k not in m.metabolites]
kg = load_kegg_links(D + 'kegg_ppu_pathway.tsv')
gp = {g.id: kg.get(g.id, set()) for g in m.genes}
rows = rescue_audit(m, SUP, gp)
pd.DataFrame(rows).to_csv('results/rescue_audit_putida_glucose.csv', index=False)
lab = {}
for r in rows:
    g = r['gene']; lab[g] = 'on_pathway' if (lab.get(g) == 'on_pathway' or r['label'] == 'on_pathway') else 'off_pathway'
fit = pd.read_csv(D + 'Putida_gene_median_fitness.tsv', sep='\t')
medf = dict(zip(fit.locusId, fit.median_fitness))
on = [medf[g] for g in lab if lab[g] == 'on_pathway' and g in medf]
off = [medf[g] for g in lab if lab[g] == 'off_pathway' and g in medf]
n_on, n_off = sum(v == 'on_pathway' for v in lab.values()), sum(v == 'off_pathway' for v in lab.values())
out = {'wt_growth_glucose': wt, 'supplements_missing_from_model': missing_sup, 'n_rescue_pairs': len(rows),
       'genes_on': n_on, 'genes_off': n_off, 'with_fitness_on': len(on), 'with_fitness_off': len(off),
       'absent_frac_on': 1 - len(on) / n_on if n_on else None, 'absent_frac_off': 1 - len(off) / n_off if n_off else None,
       'labels': lab, 'gene_median_fitness': {g: medf.get(g) for g in lab},
       'rescue_supplements': {g: sorted({r['supplement'] for r in rows if r['gene'] == g}) for g in lab}}
if len(on) < 5 or len(off) < 5:
    out['verdict'] = 'UNDERPOWERED (fewer than 5 genes with fitness in a group)'
else:
    d = float(np.median(on) - np.median(off)); p = float(mannwhitneyu(on, off, alternative='greater').pvalue)
    out.update(median_on=float(np.median(on)), median_off=float(np.median(off)), median_diff=d, mwu_p_one_sided=p)
    out['verdict'] = 'REPLICATED' if (p < 0.05 and d > 0) else ('FALSIFIED in P. putida' if d <= 0 else 'NOT SIGNIFICANT (direction consistent)')
json.dump(out, open('results/cross_species_putida.json', 'w'), indent=1)
print(json.dumps({k: v for k, v in out.items() if k not in ('labels', 'gene_median_fitness', 'rescue_supplements')}, indent=1))
