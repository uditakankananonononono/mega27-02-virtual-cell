"""ESM-2 (t6, 8M params; Lin et al. 2023, fair-esm) mean-pooled protein embeddings for every U00096.3 CDS."""
import numpy as np, torch, esm, json
from Bio import SeqIO
rec = SeqIO.read('data/U00096.3.gb', 'genbank')
ids, seqs = [], []
for f in rec.features:
    if f.type == 'CDS' and 'translation' in f.qualifiers and 'locus_tag' in f.qualifiers:
        ids.append(f.qualifiers['locus_tag'][0]); seqs.append(f.qualifiers['translation'][0][:1022])
model, alphabet = esm.pretrained.esm2_t6_8M_UR50D(); model.eval(); bc = alphabet.get_batch_converter()
torch.set_num_threads(2)
E = np.zeros((len(ids), 320), dtype=np.float32)
order = np.argsort([len(s) for s in seqs])
B = 16
with torch.no_grad():
    for k in range(0, len(order), B):
        idx = order[k:k + B]
        _, _, toks = bc([(ids[i], seqs[i]) for i in idx])
        rep = model(toks, repr_layers=[6])['representations'][6]
        for j, i in enumerate(idx):
            L = len(seqs[i]); E[i] = rep[j, 1:L + 1].mean(0).numpy()
        if k % 800 == 0: print(k, len(ids), flush=True)
np.save('results/esm2_t6_embeddings.npy', E); json.dump(ids, open('results/esm2_t6_ids.json', 'w'))
print('done', E.shape)
