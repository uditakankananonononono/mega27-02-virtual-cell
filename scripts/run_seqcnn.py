"""Train + evaluate the sequence CNN and k-mer baseline (background job)."""
import sys, time, json
sys.path.insert(0, '.')
import warnings; warnings.filterwarnings('ignore')
import numpy as np
from vcell import data as vd, seqcnn as vs

cds = vs.extract_cds('data/U00096.3.gb')
print('CDS extracted:', len(cds), flush=True)
cds.drop(columns=['cds']).to_csv('data/cds_index.csv', index=False)
labels = vd.essentiality_labels(vd.parse_gerdes_s1('data/raw/gerdes_table_s1.txt'))
X_seq, X_kmer, y, meta = vs.build_arrays(cds, labels, max_len=1200)
print('labeled CDS:', len(y), '| essential:', int(y.sum()), flush=True)
meta.to_csv('data/seq_dataset_meta.csv', index=False)
np.savez_compressed('data/seq_arrays.npz', X_kmer=X_kmer, y=y)

t0 = time.time()
base = vs.train_eval_kmer_baseline(X_kmer, y)
print('kmer baseline:', round(base['oof_auroc'],4), 'AUROC |', round(time.time()-t0,1), 's', flush=True)

t0 = time.time()
cnn = vs.train_eval_cnn(X_seq, y)
print('CNN:', round(cnn['oof_auroc'],4), 'AUROC |', round(time.time()-t0,1), 's', flush=True)

np.save('results/seqcnn_oof_scores.npy', cnn['oof_scores'])
np.save('results/kmer_oof_scores.npy', base['oof_scores'])
out = {'kmer_baseline': {k: v for k, v in base.items() if k != 'oof_scores'},
       'cnn': {k: v for k, v in cnn.items() if k != 'oof_scores'}}
json.dump(out, open('results/seqcnn_results.json', 'w'), indent=2)
print(json.dumps(out, indent=2), flush=True)
