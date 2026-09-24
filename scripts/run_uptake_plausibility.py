"""Exploratory: TCDB/ChEBI uptake plausibility of rescue-audit supplements (notes/prereg_uptake_plausibility.md)."""
import json, re, urllib.request, urllib.parse, pandas as pd
FOLLOW = ('is tautomer of', 'is conjugate acid of', 'is conjugate base of', 'is enantiomer of', 'has part')
m = {x['id']: x for x in json.load(open('data/iJO1366.json'))['metabolites']}
aud = pd.read_csv('results/rescue_audit_iML1515_M9.csv')
def graph(cid):
    iri = urllib.parse.quote(urllib.parse.quote(f'http://purl.obolibrary.org/obo/{cid.replace(":", "_")}', safe=''), safe='')
    try: return json.load(urllib.request.urlopen(f'https://www.ebi.ac.uk/ols4/api/ontologies/chebi/terms/{iri}/graph', timeout=30))
    except Exception: return {'nodes': [], 'edges': []}
iri2c = lambda s: s.rsplit('/', 1)[-1].replace('_', ':')
sub = pd.read_csv('data/external/tcdb/tcdb_substrates.tsv', sep='\t', header=None, names=['tc', 'subs'])
tc2chebi = {r.tc: set(re.findall(r'CHEBI:\d+', str(r.subs))) for r in sub.itertuples()}
up = pd.read_csv('data/external/uniprot_ecoli_tcdb.tsv', sep='\t').fillna('')
k12 = {}
for r in up.itertuples():
    for tc in [t for t in r.TCDB.split(';') if t]: k12.setdefault(tc, []).append(r._2 or r.Entry)
out = {}
for s in sorted(aud.supplement.unique()):
    base = set(m[s]['annotation'].get('chebi', [])); exp = set(base)
    for c in base:
        g = graph(c)
        for e in g['edges']:
            if e['label'] in FOLLOW:
                exp |= {iri2c(e['source']), iri2c(e['target'])}
    for c in list(exp):
        iri = urllib.parse.quote(urllib.parse.quote(f'http://purl.obolibrary.org/obo/{c.replace(":", "_")}', safe=''), safe='')
        try:
            t = json.load(urllib.request.urlopen(f'https://www.ebi.ac.uk/ols4/api/ontologies/chebi/terms/{iri}', timeout=30))
            exp |= {a if a.startswith('CHEBI:') else 'CHEBI:' + a.split(':')[-1] for a in t.get('annotation', {}).get('alt_id', []) + t.get('annotation', {}).get('has_alternative_id', [])}
        except Exception: pass
    hits = {tc: k12[tc] for tc in k12 if tc2chebi.get(tc, set()) & exp}
    lab = aud[aud.supplement == s].label.value_counts().to_dict()
    out[s] = {'n_chebi_model': len(base), 'n_chebi_expanded': len(exp), 'k12_tcdb_systems': hits,
              'uptake_supported': bool(hits), 'rescues': lab}
exp_pre = {'thm_c': True, 'pnto__R_c': True, 'btn_c': True, 'amet_c': False, 'nad_c': False, 'thf_c': False, 'pydx5p_c': False}
res = {'per_supplement': out, 'n_k12_tc_systems': len(k12), 'n_tcdb_systems_with_substrates': len(tc2chebi),
       'matches_expectation': {s: out[s]['uptake_supported'] == exp_pre[s] for s in out},
       'off_pathway_rescues_without_uptake': int(sum(v['rescues'].get('off_pathway', 0) for v in out.values() if not v['uptake_supported'])),
       'off_pathway_rescues_total': int((aud.label == 'off_pathway').sum()),
       'note': 'exploratory, descriptive, n=7 supplements; TCDB substrate curation is incomplete, so absence is weaker evidence than presence'}
json.dump(res, open('results/uptake_plausibility.json', 'w'), indent=1)
for s, v in out.items(): print(s, v['uptake_supported'], v['n_chebi_expanded'], list(v['k12_tcdb_systems'].items())[:3], v['rescues'])
print(res['matches_expectation'], res['off_pathway_rescues_without_uptake'], '/', res['off_pathway_rescues_total'])
