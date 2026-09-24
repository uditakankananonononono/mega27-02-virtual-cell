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
import os
CK = 'results/esm2_t6_partial.npz'
E = np.zeros((len(ids), 320), dtype=np.float32); done = np.zeros(len(ids), bool)
if os.path.exists(CK):
    z = np.load(CK); E, done = z['E'], z['done']
order = [i for i in np.argsort([len(s) for s in seqs]) if not done[i]]
batches, cur, toks_budget = [], [], 3000  # token-budget batching (memory-safe on 2 GB)
for i in order:
    if cur and (len(cur) + 1) * (len(seqs[i]) + 2) > toks_budget:
        batches.append(cur); cur = []
    cur.append(i)
if cur: batches.append(cur)
with torch.no_grad():
    for k, idx in enumerate(batches):
        _, _, toks = bc([(ids[i], seqs[i]) for i in idx])
        rep = model(toks, repr_layers=[6])['representations'][6]
        for j, i in enumerate(idx):
            L = len(seqs[i]); E[i] = rep[j, 1:L + 1].mean(0).numpy()
            done[i] = True
        if k % 25 == 0:
            np.savez(CK, E=E, done=done); print(int(done.sum()), len(ids), flush=True)
np.save('results/esm2_t6_embeddings.npy', E); json.dump(ids, open('results/esm2_t6_ids.json', 'w'))
print('done', E.shape)
