"""iJO1366 deletion screen under a rich (all-importable) aerobic medium -
matches the rich-LB conditions of the Gerdes footprinting experiment."""
import sys
sys.path.insert(0, '.')
import warnings; warnings.filterwarnings('ignore')
from vcell import metabolism as vm
m = vm.load_model('data/iJO1366.json')
# rich medium: every importable exchange open, oxygen open
for rxn in m.exchanges:
    rxn.lower_bound = -1000.0
m._vcell_default_medium = {r.id: -r.lower_bound for r in m.exchanges if r.lower_bound < 0}
ess = vm.gene_essentiality_scan(m, processes=2)
ess.to_csv('results/ijo1366_gene_essentiality_rich.csv', index=False)
print('rich scan:', len(ess), 'genes,', int(ess.predicted_essential.sum()), 'essential, wt growth', round(float(ess.wt_growth_h.iloc[0]),3))
