"""Train + evaluate the metabolic-graph GNN."""
import sys, json, time
sys.path.insert(0, '.')
import warnings; warnings.filterwarnings('ignore')
import numpy as np
from vcell import metabolism as vm, data as vd, graphgnn as vg

model = vm.load_model('data/iJO1366.json')
labels = vd.essentiality_labels(vd.parse_gerdes_s1('data/raw/gerdes_table_s1.txt'))
model_genes = [g.id for g in model.genes if g.id.startswith('b')]
meta = labels[labels.bnumber.isin(model_genes)].reset_index(drop=True)
genes = meta['bnumber'].tolist()
y = meta['essential'].to_numpy()
cds_index = __import__('pandas').read_csv('data/cds_index.csv')
print('graph genes:', len(genes), 'essential:', int(y.sum()), flush=True)
t0 = time.time()
A = vg.build_gene_graph(model, genes)
print('edges:', int(A.sum() // 2), '| build secs', round(time.time()-t0, 1), flush=True)
np.save('data/gene_graph_adj.npy', A)
X = vg.node_features(model, genes, cds_index)
t0 = time.time()
res = vg.train_eval_gnn(A, X, y)
print('GNN:', round(res['oof_auroc'],4), 'AUROC |', round(time.time()-t0,1), 's', flush=True)
np.save('results/gnn_oof_scores.npy', res['oof_scores'])
meta.to_csv('data/gnn_meta.csv', index=False)
out = {k: v for k, v in res.items() if k != 'oof_scores'}
json.dump(out, open('results/gnn_results.json', 'w'), indent=2)
print(json.dumps(out, indent=2), flush=True)
