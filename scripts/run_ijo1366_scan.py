"""Genome-scale single-gene-deletion screen on iJO1366 (1367 genes)."""
import sys, time
sys.path.insert(0, '.')
import warnings; warnings.filterwarnings('ignore')
from vcell import metabolism as vm
m = vm.load_model('data/iJO1366.json')
t0 = time.time()
ess = vm.gene_essentiality_scan(m, processes=2)
ess.to_csv('results/ijo1366_gene_essentiality.csv', index=False)
print('iJO1366 scan done:', len(ess), 'genes,', int(ess.predicted_essential.sum()),
      'predicted essential,', round(time.time()-t0, 1), 's')
