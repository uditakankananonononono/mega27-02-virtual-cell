"""Fetch AlphaFold DB summary metrics for model genes (resume-safe) -> data/external/alphafold/af_metrics.tsv"""
import json, os, re, time, urllib.request, pandas as pd
up = pd.read_csv('data/external/uniprot_ecoli.tsv', sep='\t'); b2u = {}
for _, x in up.iterrows():
    for t in str(x['Gene Names (ordered locus)']).split():
        if re.fullmatch(r'b\d{4}', t): b2u.setdefault(t, x['Entry'])
genes = pd.read_csv('results/aligned_predictions_all_models.csv').bnumber
out = 'data/external/alphafold/af_metrics.tsv'
done = set(pd.read_csv(out, sep='\t').bnumber) if os.path.exists(out) else set()
with open(out, 'a') as f:
    if not done: f.write('bnumber\tuniprot\tmodel_id\tplddt_mean\tplddt_frac_vlow\n')
    for b in genes:
        if b in done: continue
        u = b2u.get(b)
        row = [b, u or '', '', '', '']
        if u:
            try:
                d = json.load(urllib.request.urlopen(f'https://alphafold.ebi.ac.uk/api/prediction/{u}', timeout=30))[0]
                row = [b, u, d['modelEntityId'], d['globalMetricValue'], d['fractionPlddtVeryLow']]
            except Exception as e:
                row = [b, u, 'ERR:' + type(e).__name__, '', '']
        f.write('\t'.join(map(str, row)) + '\n'); f.flush()
