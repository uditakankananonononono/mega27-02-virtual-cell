"""GO enrichment (goatools, Fisher + BH) of the v2 ensemble's hard errors against the benchmark gene universe.
FN = essential genes in the bottom half of v2_lr OOF scores; FP = non-essential genes in the top 10%."""
import gzip, json, collections, numpy as np, pandas as pd
from Bio import SeqIO
import statsmodels.sandbox.stats.multicomp as _smc; from statsmodels.stats.multitest import multipletests as _mt; _smc.multipletests = _mt  # goatools compat shim
from goatools.obo_parser import GODag
from goatools.go_enrichment import GOEnrichmentStudy
rec = SeqIO.read('data/U00096.3.gb', 'genbank'); sym2b = {}
for f in rec.features:
    if f.type == 'CDS' and 'locus_tag' in f.qualifiers:
        b = f.qualifiers['locus_tag'][0]
        for key in ('gene', 'gene_synonym'):
            for v in f.qualifiers.get(key, []):
                for s in v.replace(';', ' ').split():
                    sym2b.setdefault(s.lower(), b)
assoc = collections.defaultdict(set); unmapped = 0
for line in gzip.open('data/external/ecocyc.gaf.gz', 'rt'):
    if line.startswith('!'): continue
    c = line.rstrip('\n').split('\t')
    if 'NOT' in c[3]: continue
    b = sym2b.get(c[2].lower())
    if b is None:
        for s in c[10].split('|'):
            b = sym2b.get(s.lower())
            if b: break
    if b is None: unmapped += 1; continue
    assoc[b].add(c[4])
dag = GODag('data/external/go-basic.obo', load_obsolete=False, prt=None)
oof = pd.read_csv('results/ensemble_v2_oof.csv')
pop = [b for b in oof.bnumber if b in assoc]
o = oof.set_index('bnumber').loc[pop]
r = o.v2_lr.rank(pct=True)
fn = list(o.index[(o.essential == 1) & (r < 0.5)]); fp = list(o.index[(o.essential == 0) & (r > 0.9)])
g = GOEnrichmentStudy(pop, {k: assoc[k] for k in pop}, dag, propagate_counts=True, alpha=0.05, methods=['fdr_bh'], log=None)
out = {'n_universe': len(pop), 'n_annotated_genes_total': len(assoc), 'gaf_rows_unmapped': unmapped, 'n_fn': len(fn), 'n_fp': len(fp)}
for name, study in (('false_negatives', fn), ('false_positives', fp)):
    res = [x for x in g.run_study(study, prt=None) if x.enrichment == 'e' and x.p_fdr_bh < 0.05]
    res.sort(key=lambda x: x.p_fdr_bh)
    out[name] = [{'go': x.GO, 'name': x.name, 'ns': x.NS, 'study': f'{x.study_count}/{x.study_n}', 'pop': f'{x.pop_count}/{x.pop_n}',
                  'p_fdr_bh': float(x.p_fdr_bh)} for x in res[:25]]
json.dump(out, open('results/go_error_enrichment.json', 'w'), indent=1)
print(json.dumps({k: (v if not isinstance(v, list) else v[:8]) for k, v in out.items()}, indent=1))
