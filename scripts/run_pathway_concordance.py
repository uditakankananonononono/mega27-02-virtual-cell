"""Pre-registered test of pathway concordance of supplement rescue (notes/prereg_pathway_concordance.md)."""
import json, collections, numpy as np
from scipy.stats import mannwhitneyu
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
f = np.array([r[2] for r in rows]); on = np.array([r[4] for r in rows]); genes = np.array([r[0] for r in rows]); grp = np.array([r[3] for r in rows])
def test(mask):
    a, b = f[mask & on], f[mask & ~on]
    if len(a) == 0 or len(b) == 0: return {'n_on': int(len(a)), 'n_off': int(len(b)), 'note': 'one side empty'}
    p = mannwhitneyu(a, b, alternative='greater').pvalue
    ug = np.unique(genes[mask]); rng = np.random.default_rng(0); ds = []
    idx = {g: np.where((genes == g) & mask)[0] for g in ug}
    for _ in range(2000):
        s = np.concatenate([idx[g] for g in rng.choice(ug, len(ug))])
        x, y = f[s][on[s]], f[s][~on[s]]
        if len(x) and len(y): ds.append(np.median(x) - np.median(y))
    return {'n_on_pairs': int(len(a)), 'n_off_pairs': int(len(b)), 'n_on_genes': int(len(set(genes[mask & on]))), 'n_off_genes': int(len(set(genes[mask & ~on]))),
            'median_on': float(np.median(a)), 'median_off': float(np.median(b)), 'frac_gt_-2_on': float((a > -2).mean()), 'frac_gt_-2_off': float((b > -2).mean()),
            'mwu_p_one_sided': float(p), 'median_diff': float(np.median(a) - np.median(b)), 'cluster_boot_ci95': [float(np.percentile(ds, 2.5)), float(np.percentile(ds, 97.5))]}
allm = np.ones(len(f), bool)
out = {'primary_all_rescued': test(allm), 'secondary_within_vitamins': test(grp == 'vitamins5'),
       'secondary_within_nonvitamin': test(grp != 'vitamins5'),
       'off_pathway_vitamin_genes': sorted(set(genes[(grp == 'vitamins5') & ~on])),
       'on_pathway_nonvitamin_genes': sorted(set(genes[(grp != 'vitamins5') & on]))}
p = out['primary_all_rescued']
out['verdict'] = ('SUPPORTED (pre-registered primary): on-pathway rescues have higher fitness, CI excludes 0'
                  if p.get('mwu_p_one_sided', 1) < 0.05 and p['cluster_boot_ci95'][0] > 0 else 'FALSIFIED on primary criterion')
json.dump(out, open('results/pathway_concordance.json', 'w'), indent=1); print(json.dumps(out, indent=1))
