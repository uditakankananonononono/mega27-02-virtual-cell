"""Corrected benchmark: does the provenance audit sort FBA essential calls correctly?
Falsifiable test: biomass_forced calls should be enriched for experimentally
NON-essential genes (Gerdes 2003) relative to network_forced calls (Fisher exact).
Then re-score iJO1366 FBA with biomass_forced genes flipped to non-essential."""
import sys, json, ast
sys.path.insert(0, '.')
import pandas as pd
from scipy.stats import fisher_exact
from vcell.benchmark import classification_report, bootstrap_ci

aud = pd.read_csv('results/essentiality_provenance_audit.csv')
al = pd.read_csv('results/aligned_predictions_all_models.csv')
m = al.merge(aud[['gene', 'verdict', 'unproducible']], left_on='bnumber', right_on='gene', how='left')
called = m[m.verdict.notna()].copy()
from vcell.compound_classes import gene_class
called['uclass'] = [gene_class(ast.literal_eval(u)) if isinstance(u, str) else 'none' for u in called.unproducible]
called['stratum'] = called.verdict + ':' + called.uclass
strat = called.groupby('stratum').essential.agg(['count', 'sum']).reset_index()
cof = called[(called.verdict == 'biomass_forced') & (called.uclass == 'cofactor_only')]
rest = called[~called.index.isin(cof.index)]
tbl2 = [[int((cof.essential == 0).sum()), int((cof.essential == 1).sum())],
        [int((rest.essential == 0).sum()), int((rest.essential == 1).sum())]]
odds2, p2 = fisher_exact(tbl2, alternative='greater') if len(cof) else (float('nan'), float('nan'))
tab = called.groupby('verdict').essential.agg(['count', 'sum']).rename(columns={'sum': 'gerdes_essential'})
tab['gerdes_essential_frac'] = tab.gerdes_essential / tab['count']
bf = called[called.verdict == 'biomass_forced']; nf = called[called.verdict != 'biomass_forced']
tbl = [[int((bf.essential == 0).sum()), int((bf.essential == 1).sum())],
       [int((nf.essential == 0).sum()), int((nf.essential == 1).sum())]]
odds, p = fisher_exact(tbl, alternative='greater') if len(bf) else (float('nan'), float('nan'))
m['fba_min_pred'] = (m.fba_min > 0.5).astype(int)
m = m.merge(called[['bnumber','uclass']], on='bnumber', how='left')
m['fba_min_corr'] = m.fba_min.where(~((m.verdict == 'biomass_forced') & (m.uclass == 'cofactor_only')), 0.0)
m['fba_min_corr_pred'] = (m.fba_min_corr > 0.5).astype(int)
base = classification_report(m, 'fba_min', pred_col='fba_min_pred')
corr = classification_report(m, 'fba_min_corr', pred_col='fba_min_corr_pred')
corr.update(bootstrap_ci(m, 'fba_min_corr'))
# compound census
comp = {}
for u in bf.unproducible:
    for c in ast.literal_eval(u):
        comp[c] = comp.get(c, 0) + 1
out = {'strata': strat.to_dict(orient='records'), 'fisher_cofactor_only_vs_rest': tbl2,
       'fisher_cofactor_odds': odds2, 'fisher_cofactor_p_one_sided': p2,
       'verdict_table': tab.reset_index().to_dict(orient='records'),
       'fisher_2x2_[bf_noness,bf_ess],[nf_noness,nf_ess]': tbl,
       'fisher_odds_ratio': odds, 'fisher_p_one_sided': p,
       'fba_baseline_on_aligned': base, 'fba_corrected': corr,
       'biomass_forced_compound_census': dict(sorted(comp.items(), key=lambda x: -x[1])),
       'biomass_forced_genes': bf[['bnumber', 'gene_name', 'essential', 'unproducible']].to_dict(orient='records')}
json.dump(out, open('results/corrected_benchmark.json', 'w'), indent=2, default=float)
print(json.dumps({k: v for k, v in out.items() if k != 'biomass_forced_genes'}, indent=1, default=float))
print(bf[['bnumber', 'gene_name', 'essential', 'unproducible']].to_string())
