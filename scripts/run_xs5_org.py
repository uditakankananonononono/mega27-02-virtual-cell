"""Amendment 5 (notes/prereg_cross_species.md): per-organism CarveMe + production-blocked rescue labels. Resume-safe."""
import json, logging, os, subprocess, sys, warnings
import numpy as np, pandas as pd, cobra
from vcell.rescue import rescue_audit_production
warnings.filterwarnings('ignore'); logging.disable(logging.CRITICAL)
D = 'data/cross/'
M9 = ['ca2', 'cl', 'cobalt2', 'cu2', 'fe2', 'fe3', 'glc__D', 'h2o', 'h', 'k', 'mg2', 'mn2', 'mobd', 'na1', 'nh4', 'ni2', 'o2', 'pi', 'so4', 'zn2']
SUP = ['btn_c', 'thm_c', 'pnto__R_c', 'thf_c', 'nad_c', 'amet_c', 'pydx5p_c']
FIG = json.load(open('/tmp/fig.json'))
env = dict(os.environ, PATH=os.path.expanduser('~/bin') + ':' + os.environ['PATH'])
for org in sys.argv[1:]:
    out = f'results/xs5/{org}.json'
    if os.path.exists(out): continue
    xml = D + f'carve/{org}.xml'
    if not os.path.exists(xml):
        faa = f'/tmp/{org}.faa'
        subprocess.run(f"zcat {D}FB_aaseqs.gz | awk -v o='{org}' '/^>/{{split(substr($0,2),a,\":\"); p=(a[1]==o); if(p) print \">\"a[2]; next}} p' > {faa}", shell=True, check=True)
        subprocess.run(['carve', faa, '-o', xml, '-g', 'M9', '-i', 'M9', '--fbc2', '--solver', 'scip'], env=env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        os.remove(faa)
    if not os.path.exists(xml):
        json.dump({'org': org, 'excluded': 'CarveMe failed'}, open(out, 'w')); continue
    fitf = D + f'{org}_gene_median_fitness.tsv'
    if not os.path.exists(fitf):
        fid = [f['id'] for f in FIG['files'] if f['name'] == f'db.StrainFitness.{org}.gz'][0]
        subprocess.run(['curl', '-sL', '-m', '600', '-A', 'Mozilla/5.0', '-o', D + f'{org}_strainfit.gz', f'https://ndownloader.figshare.com/files/{fid}'], check=True)
        subprocess.run([sys.executable, 'scripts/aggregate_fb_strainfit.py', org], check=True)
        os.remove(D + f'{org}_strainfit.gz')
    m = cobra.io.read_sbml_model(xml)
    m.medium = {f'EX_{c}_e': (10.0 if c == 'glc__D' else 1000.0) for c in M9 if f'EX_{c}_e' in m.reactions}
    wt = m.slim_optimize(error_value=0.0)
    if wt < 1e-6:
        json.dump({'org': org, 'excluded': 'no growth on M9 glucose', 'wt': wt}, open(out, 'w')); continue
    rows = rescue_audit_production(m, SUP)
    fit = pd.read_csv(fitf, sep='\t', dtype={'locusId': str}).dropna(subset=['median_fitness'])
    pct = dict(zip(fit.locusId, fit.median_fitness.rank(pct=True))); medf = dict(zip(fit.locusId, fit.median_fitness))
    genes = {}
    for r in rows:
        g = genes.setdefault(r['gene'], {'label': 'off_pathway', 'supplements': []})
        g['supplements'].append(r['supplement'])
        if r['label'] == 'on_pathway': g['label'] = 'on_pathway'
    for gid, g in genes.items():
        g['median_fitness'] = medf.get(gid); g['percentile'] = pct.get(gid)
    json.dump({'org': org, 'wt': wt, 'n_genes_model': len(m.genes), 'n_rescue_pairs': len(rows),
               'missing_supplements': [s for s in SUP if s not in m.metabolites], 'genes': genes}, open(out, 'w'), indent=1)
    print(org, 'done', len(rows), flush=True)
