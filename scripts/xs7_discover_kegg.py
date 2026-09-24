"""Amendment 7 organism discovery: KEGG find/genes on 3 sample FB locusIds per organism; criterion (a)."""
import glob, json, os, time, urllib.request, pandas as pd
SKIP = {'Caulo', 'Cola', 'PS', 'Ponti', 'Smeli', 'SyringaeB728a'}
out = {}
def get(u):
    for _ in range(3):
        try: return urllib.request.urlopen(u, timeout=60).read().decode()
        except Exception: time.sleep(2)
    return ''
for f in sorted(glob.glob('results/xs5/*.json')):
    org = json.load(open(f))['org']
    if org in SKIP: continue
    loc = pd.read_csv(f'data/cross/{org}_gene_median_fitness.tsv', sep='\t', dtype={'locusId': str}).locusId.tolist()
    codes = {}
    for l in [loc[0], loc[len(loc) // 2], loc[-1]]:
        for line in get(f'https://rest.kegg.jp/find/genes/{l}').splitlines():
            gid = line.split('\t')[0]
            if ':' in gid and gid.split(':', 1)[1] == l: codes[gid.split(':')[0]] = codes.get(gid.split(':')[0], 0) + 1
        time.sleep(0.4)
    rec = {'sample_hits': codes}
    if codes:
        kc = max(codes, key=codes.get); rec['kegg_code'] = kc
        path = f'data/cross/kegg_{kc}_pathway.tsv'
        if not os.path.exists(path):
            txt = get(f'https://rest.kegg.jp/link/pathway/{kc}'); open(path, 'w').write(txt); time.sleep(0.4)
        ids = {l.split('\t')[0].split(':', 1)[1] for l in open(path) if '\t' in l}
        frac = len(ids & set(loc)) / max(1, len(ids)); rec.update(n_kegg_pathway_genes=len(ids), frac_kegg_ids_in_fb=frac, included=frac >= 0.5)
    else: rec['included'] = False; rec['reason'] = 'no exact KEGG gene id match for sample loci'
    out[org] = rec; print(org, rec.get('kegg_code'), rec.get('frac_kegg_ids_in_fb'), rec['included'], flush=True)
json.dump(out, open('results/xs7_kegg_discovery.json', 'w'), indent=1)
