"""Full essentiality-provenance audit of iJO1366 essential genes (parallel, resumable).
Model loaded once per worker; rows appended to a partial CSV as they finish."""
import sys, time, json, os, csv
sys.path.insert(0, '.')
import warnings; warnings.filterwarnings('ignore')
import pandas as pd
from multiprocessing import Pool

_M = None
def _init():
    global _M, _BM, _BID
    from vcell import metabolism as vm, audit as va
    _M = vm.load_model('data/iJO1366.json')
    _BM, _BID = va.biomass_metabolites(_M)

def one(gene):
    from vcell import audit as va
    return va.audit_gene(_M, gene, _BID, _BM)

if __name__ == '__main__':
    part = 'results/essentiality_provenance_audit.partial.csv'
    ess = pd.read_csv('results/ijo1366_gene_essentiality.csv')
    genes = ess[ess.predicted_essential].gene.tolist()
    done = set(pd.read_csv(part).gene) if os.path.exists(part) else set()
    todo = [g for g in genes if g not in done]
    print('auditing', len(genes), 'genes;', len(todo), 'remaining', flush=True)
    t0 = time.time()
    new = not os.path.exists(part)
    with open(part, 'a', newline='') as fh, Pool(2, initializer=_init) as p:
        w = csv.DictWriter(fh, fieldnames=['gene', 'verdict', 'rescued_growth', 'wt_growth', 'unproducible'])
        if new: w.writeheader()
        for i, r in enumerate(p.imap_unordered(one, todo), 1):
            r.setdefault('wt_growth', ''); r['unproducible'] = repr(r['unproducible'])
            w.writerow(r); fh.flush()
            if i % 10 == 0: print(f'{i}/{len(todo)} {time.time()-t0:.0f}s', flush=True)
    df = pd.read_csv(part)
    df.to_csv('results/essentiality_provenance_audit.csv', index=False)
    summ = df.verdict.value_counts().to_dict(); summ['seconds'] = round(time.time() - t0, 1)
    json.dump(summ, open('results/essentiality_provenance_summary.json', 'w'), indent=2)
    print(json.dumps(summ, indent=2), flush=True)
