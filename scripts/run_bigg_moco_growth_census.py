"""Prereg VC2-R2: MoCo biomass-constituent removal essentiality-flip census.

For each MoCo-positive model in results/bigg_model_survey.csv:
  1. verify SHA-256 of the stored snapshot against the survey CSV
  2. baseline FBA; scorable iff growth >= 1e-6
  3. full single-gene deletion at baseline -> essential gene set (growth < 1e-6)
  4. zero MoCo-named constituents in the positive-objective reaction(s);
     retest ONLY the baseline-essential genes; flip iff growth >= 0.95 * WT
Positive control: iJO1366 must reproduce the moaD-cluster rescue.
Resumable: per-model JSON in results/bigg_moco_growth/ skips finished models.
"""
import csv, hashlib, json, sys, warnings
from pathlib import Path
warnings.filterwarnings('ignore')
ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'data'/'bigg_survey'
OUTD = ROOT/'results'/'bigg_moco_growth'; OUTD.mkdir(parents=True, exist_ok=True)
MOCO_TERMS = ('moco', 'molybdo', 'bmoco', 'mobd')
GROWTH_EPS = 1e-6

def moco_mets(model):
    out = []
    for m in model.metabolites:
        s = (m.id + ' ' + (m.name or '')).lower()
        if any(t in s for t in MOCO_TERMS):
            out.append(m)
    return out

def main(limit=None, only=None):
    import cobra
    survey = list(csv.DictReader((ROOT/'results'/'bigg_model_survey.csv').open()))
    targets = [r for r in survey if r['moco_biomass'] == 'True']
    if only:
        targets = [r for r in targets if r['accession'] in only]
    print(f'{len(targets)} MoCo-positive models targeted')
    for r in targets:
        acc = r['accession']
        out = OUTD/f'{acc}.json'
        if out.exists():
            continue
        snap = DATA/f'{acc}.json'
        rec = {'accession': acc, 'status': 'pending'}
        try:
            digest = hashlib.sha256(snap.read_bytes()).hexdigest()
            if digest != r['sha256']:
                rec.update(status='hash_mismatch', expected=r['sha256'], got=digest)
                out.write_text(json.dumps(rec, indent=2)); continue
            model = cobra.io.load_json_model(str(snap))
            wt = model.slim_optimize()
            rec['wt_growth'] = None if wt is None else float(wt)
            if wt is None or wt < GROWTH_EPS:
                rec['status'] = 'non_scorable_no_baseline_growth'
                out.write_text(json.dumps(rec, indent=2)); continue
            # baseline essentiality
            base = cobra.flux_analysis.single_gene_deletion(model, processes=1)
            def gid(v):
                return sorted(v)[0] if isinstance(v, (set, frozenset)) else v
            ess = sorted(gid(v) for v, gr in zip(base['ids'], base['growth'])
                         if gr is not None and gr < GROWTH_EPS)
            rec['n_genes'] = len(model.genes)
            rec['n_essential_baseline'] = len(ess)
            # MoCo-removed objective
            mets = moco_mets(model)
            rec['moco_metabolite_ids'] = sorted(m.id for m in mets)
            obj_rxns = [rx for rx in model.reactions if rx.objective_coefficient and rx.objective_coefficient > 0]
            rec['objective_reactions'] = [rx.id for rx in obj_rxns]
            removed = 0
            for rx in obj_rxns:
                delta = {}
                for met in mets:
                    if met in rx.metabolites and rx.metabolites[met] < 0:
                        delta[met] = rx.metabolites[met]; removed += 1
                if delta:
                    rx.subtract_metabolites(delta)
            rec['moco_constituents_zeroed'] = removed
            flips, retained = [], {}
            if ess:
                with model:
                    # objective constituents already zeroed in this context
                    for g in ess:
                        with model:
                            model.genes.get_by_id(g).knock_out()
                            gr = model.slim_optimize()
                        gr = 0.0 if gr is None or gr < 0 else float(gr)
                        if gr >= 0.95 * rec['wt_growth']:
                            flips.append(g)
                        retained[g] = gr
            rec['status'] = 'scored'
            rec['flips'] = flips
            rec['n_flips'] = len(flips)
            rec['retained_growth_rescued'] = {g: retained[g] for g in flips}
        except Exception as e:
            rec['status'] = 'error'
            rec['error'] = type(e).__name__ + ': ' + str(e)[:300]
        out.write_text(json.dumps(rec, indent=2))
        print(acc, rec['status'], 'ess=', rec.get('n_essential_baseline'), 'flips=', rec.get('n_flips'), flush=True)
        if limit:
            limit -= 1
            if limit <= 0:
                break

if __name__ == '__main__':
    lim = int(sys.argv[1]) if len(sys.argv) > 1 else None
    only = set(sys.argv[2].split(',')) if len(sys.argv) > 2 else None
    main(lim, only)
