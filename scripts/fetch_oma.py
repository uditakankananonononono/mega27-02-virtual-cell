"""Count OMA 1:1 orthologs per model gene (resume-safe) -> data/external/oma/oma_1to1_counts.tsv"""
import time, gzip, json, os, re, urllib.request, pandas as pd
up = pd.read_csv('data/external/uniprot_ecoli.tsv', sep='\t'); b2u = {}
for _, x in up.iterrows():
    for t in str(x['Gene Names (ordered locus)']).split():
        if re.fullmatch(r'b\d{4}', t): b2u.setdefault(t, x['Entry'])
genes = pd.read_csv('results/aligned_predictions_all_models.csv').bnumber
os.makedirs('data/external/oma', exist_ok=True); out = 'data/external/oma/oma_1to1_counts.tsv'
done = set()
if os.path.exists(out):
    t = pd.read_csv(out, sep='\t'); t = t[t.status.astype(str).str.startswith('ok:')]; t.to_csv(out, sep='\t', index=False); done = set(t.bnumber)
with open(out, 'a') as f:
    if not done: f.write('bnumber\tuniprot\tn_orth_1to1\tn_species_codes\tstatus\n')
    for b in genes:
        if b in done: continue
        u = b2u.get(b); row = [b, u or '', 0, 0, 'no_uniprot']
        if u:
            try:
                H = {'Accept-Encoding': 'gzip', 'User-Agent': 'curl/8.5 (mega27 research)'}
                def get(url):
                    r = urllib.request.urlopen(urllib.request.Request(url, headers=H), timeout=60); raw = r.read()
                    return json.loads(gzip.decompress(raw) if r.headers.get('Content-Encoding') == 'gzip' else raw)
                xr = [x for x in get(f'https://omabrowser.org/api/xref/?search={u}') if str(x.get('omaid', '')).startswith('ECOLI')]
                if not xr: row = [b, u, 0, 0, 'no_ECOLI_entry']
                else:
                    d = get(f"https://omabrowser.org/api/protein/{xr[0]['entry_nr']}/orthologs/?rel_type=1:1")
                    row = [b, u, len(d), len({x['species']['code'] for x in d}), 'ok:' + xr[0]['omaid']]
            except urllib.error.HTTPError as e:
                row = [b, u, 0, 0, f'http{e.code}']
            except Exception as e:
                row = [b, u, '', '', 'ERR:' + type(e).__name__]
        f.write('\t'.join(map(str, row)) + '\n'); f.flush(); time.sleep(0.5)
