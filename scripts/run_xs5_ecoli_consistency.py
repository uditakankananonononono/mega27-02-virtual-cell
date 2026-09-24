"""Amendment 5 secondary: model-internal (production-blocked) label applied to E. coli iML1515 on the H.7 Bernstein
rescue rows (consistency check with results/pathway_concordance.json). Label per (gene, supplement) from
rescue_audit_production on iML1515, M9 glucose, same supplement ids as the cross-species runs."""
import json, collections, logging, warnings, numpy as np, cobra
from scipy.stats import mannwhitneyu
from vcell.rescue import rescue_audit_production
warnings.filterwarnings('ignore'); logging.disable(logging.CRITICAL)
B = 'results/bernstein/'
kg = collections.defaultdict(set)
for line in open('data/external/kegg_eco_pathway.tsv'):
    g, p = line.split('\t'); kg[g.replace('eco:', '')].add(p.strip().replace('path:eco', ''))
PW = {'btn': {'00780'}, 'thm': {'00730'}, 'pnto__R': {'00770'}, 'thf': {'00790', '00670'}, 'nad': {'00760'},
      'amet': {'00270'}, 'pydx5p': {'00750'}, '10fthf': {'00670', '00790'}, 'mlthf': {'00670', '00790'}}
VIT = set().union(*[PW[v] for v in ('btn', 'thm', 'pnto__R', 'thf', 'nad')])
rows = []  # (gene, carbon_idx, fitness, supplement_group, on_pathway)
sb = np.load(B + 'iML1515_base_sim.npy'); sv = np.load(B + 'iML1515_bernstein_vitamins_sim.npy'); fb = np.load(B + 'iML1515_base_fit.npy')
gb = json.load(open(B + 'iML1515_base_ids.json'))['genes']
for gi, ci in np.argwhere((sb <= 0.001) & (sv > 0.001)):
    rows.append((gb[gi], ci, fb[gi, ci], 'vitamins5', bool(kg[gb[gi]] & VIT)))
A = np.load(B + 'attrib_allcorr_fixed.npz'); pairs, R, cofs = A['pairs'], A['R'], list(A['cofs'])
fa = np.load(B + 'iML1515_bernstein_allcorr_fixed_fit.npy'); ga = json.load(open(B + 'iML1515_bernstein_allcorr_fixed_ids.json'))['genes']
for k, (gi, ci) in enumerate(pairs):
    hit = [c for j, c in enumerate(cofs) if R[k, j]]
    if not hit: continue
    on = any(kg[ga[gi]] & PW.get(c, set()) for c in hit)
    rows.append((ga[gi], ci, fa[gi, ci], '+'.join(hit), on))

M9 = ['ca2', 'cl', 'cobalt2', 'cu2', 'fe2', 'fe3', 'glc__D', 'h2o', 'h', 'k', 'mg2', 'mn2', 'mobd', 'na1', 'nh4', 'ni2', 'o2', 'pi', 'so4', 'zn2']
m = cobra.io.read_sbml_model('external/E_coli_GEM_validation/Models/iML1515.xml')
m.medium = {f'EX_{c}_e': (10.0 if c == 'glc__D' else 1000.0) for c in M9 if f'EX_{c}_e' in m.reactions}
VIT5 = ['btn', 'thm', 'pnto__R', 'thf', 'nad']
needed = collections.defaultdict(set)
for g, ci, fit, grp, _ in rows:
    for c in (VIT5 if grp == 'vitamins5' else grp.split('+')): needed[g].add(c + '_c')
sups = sorted({c for v in needed.values() for c in v if c in m.metabolites})
lab = {}
for g in sorted(needed):
    if g not in m.genes: continue
    for r in rescue_audit_production(m, [s for s in sups if s in needed[g]], genes=[g]):
        lab[(g, r['supplement'])] = r['label']
f, on, keep = [], [], []
for g, ci, fit, grp, _ in rows:
    ls = [lab.get((g, c + '_c')) for c in (VIT5 if grp == 'vitamins5' else grp.split('+'))]
    ls = [x for x in ls if x]
    if not ls: continue  # the supplement does not rescue this knockout on glucose M9 in this model
    f.append(fit); on.append('on_pathway' in ls); keep.append(g)
f, on = np.array(f), np.array(on)
a, b = f[on], f[~on]
res = {'n_rows_H7': len(rows), 'n_rows_labelled': int(len(f)), 'n_on_rows': int(on.sum()), 'n_off_rows': int((~on).sum()),
       'n_on_genes': len({g for g, o in zip(keep, on) if o}), 'n_off_genes': len({g for g, o in zip(keep, on) if not o})}
if len(a) and len(b):
    res.update(median_on=float(np.median(a)), median_off=float(np.median(b)), median_diff=float(np.median(a) - np.median(b)),
               mwu_p_one_sided=float(mannwhitneyu(a, b, alternative='greater').pvalue))
res['consistent_with_H7'] = bool(res.get('median_diff', 0) > 0 and res.get('mwu_p_one_sided', 1) < 0.05)
json.dump(res, open('results/xs5_ecoli_consistency.json', 'w'), indent=1); print(json.dumps(res, indent=1))
