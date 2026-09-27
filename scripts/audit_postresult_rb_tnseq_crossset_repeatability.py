"""Protocol-locked descriptive cross-set fitness-median concordance, not validation."""
import hashlib,json,math,statistics
from pathlib import Path
from scipy.stats import pearsonr,spearmanr
source=Path('/tmp/fModule_Metadata.xlsx')
x=json.loads(Path('results/postresult_2024_rb_tnseq_glucose_screen.json').read_text())
assert source.exists() and hashlib.sha256(source.read_bytes()).hexdigest()==x['source_sha256']=='4d649385ac06684482396a125f135df22a2a5060da73485b2cd14468f8cc8be1'
assert len(x['samples'])==13
samples={a['expName']:a for a in x['samples']}
assert {k for k in samples if k.startswith('set100')}=={'set100IT008','set100IT009','set100IT010'}
assert {k for k in samples if k.startswith('set101')}=={'set101IT007','set101IT008','set101IT009'}
assert all(samples[k]['media']=='M9_medium' and samples[k]['concentration_1']=='20' and samples[k]['condition_1']=='D-Glucose' for k in samples if k.startswith(('set100','set101')))
rows=[];missing=[]
for panel,genes in x['panels'].items():
 for gene,p in sorted(genes.items()):
  if not p.get('groups'):
   missing.append({'panel':panel,'gene':gene,'reason':p['status']});continue
  a=p['groups']['set100_20mM_M9_medium'];b=p['groups']['set101_20mM_M9_medium']
  assert len(a['exp_names'])==len(b['exp_names'])==3
  assert set(a['exp_names'])=={z for z in samples if z.startswith('set100')}
  assert set(b['exp_names'])=={z for z in samples if z.startswith('set101')}
  va,vb=a['fitness_median'],b['fitness_median']
  assert va is not None and vb is not None and all(math.isfinite(z) for z in (va,vb))
  rows.append({'panel':panel,'gene':gene,'set100_median':va,'set101_median':vb,'abs_median_difference':abs(va-vb),'strict_sign_agree':va*vb>0,'both_below_minus_one':va<-1 and vb<-1,'either_below_minus_one':va<-1 or vb<-1})
assert len(rows)==11 and len(missing)==3
left=[r['set100_median'] for r in rows];right=[r['set101_median'] for r in rows]
out={'status':'post-result source-table cross-set concordance only, no phenotype pass','protocol':'notes/postresult_rb_tnseq_crossset_repeatability_protocol.md','source_url':x['source'],'source_sha256':x['source_sha256'],'n_complete_genes':len(rows),'n_missing_genes':len(missing),'missing':missing,'pairs':rows,'strict_sign_agreements':sum(r['strict_sign_agree'] for r in rows),'pearson_r':float(pearsonr(left,right).statistic),'spearman_rho':float(spearmanr(left,right).statistic),'median_absolute_difference':statistics.median(r['abs_median_difference'] for r in rows),'both_below_minus_one':sum(r['both_below_minus_one'] for r in rows),'either_below_minus_one':sum(r['either_below_minus_one'] for r in rows),'limits':'11 selected gene summaries across two related source assay sets, each summary uses three samples; neither 11 independent experiments nor direct knockout survival. Set100 accession unresolved, set101 IT-index/run mapping unresolved. No p-value, no essentiality threshold and original G2 FAIL.'}
Path('results/postresult_rb_tnseq_crossset_repeatability.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print({k:out[k] for k in ('n_complete_genes','strict_sign_agreements','pearson_r','spearman_rho','median_absolute_difference','both_below_minus_one','either_below_minus_one')})
