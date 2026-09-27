"""Retrospective complete two-coefficient in-model null, fixed protocol in notes/."""
import csv, hashlib, itertools, json, sys
from pathlib import Path
import cobra
from cobra.flux_analysis import single_gene_deletion
p=Path('data/bigg_survey/iJN1463.json')
row=next(z for z in csv.DictReader(open('results/bigg_model_survey.csv')) if z['accession']=='iJN1463')
sha=hashlib.sha256(p.read_bytes()).hexdigest();assert sha==row['sha256']=='8232a1ea39af9b9b0313601a3cefca1028e50e53b2ca9c73a9adf7568bf'
m=cobra.io.load_json_model(str(p)); obj=m.reactions.get_by_id('BIOMASS_KT2440_WT3')
coef={z.id:float(c) for z,c in obj.metabolites.items() if c<0}
mo=('bmocogdp_c','mocogdp_c');target=coef[mo[0]]
assert all(abs(coef[z]-target)<1e-12 for z in mo)
eligible=sorted(z for z,c in coef.items() if z not in ('bmocogdp_c','mobd_c','mocogdp_c') and abs(c-target)<1e-12)
pairs=list(itertools.combinations(eligible,2)); assert len(eligible)==23 and len(pairs)==253
wt=float(m.slim_optimize()); base=single_gene_deletion(m,processes=1)
essential=sorted(next(iter(ids)) for ids,gr in zip(base['ids'],base['growth']) if gr is not None and gr<1e-6)
assert len(essential)==262 and abs(wt-.5861175448479794)<1e-8
prior=json.load(open('results/moco_subset_interaction_ijn1463.json'))
mo_prior=next(x for x in prior['subsets'] if tuple(x['removed'])==mo)
assert mo_prior['n_rescued']==6
outpath=Path('results/postresult_equal_coefficient_pair_null.json')
header={'status':'retrospective equal-coefficient pair null; not exchangeable biological control',
 'protocol':'notes/postresult_equal_coefficient_pair_null.md','source_url':row['source_url'],'source_sha256':sha,
 'objective':obj.id,'original_wt_growth':wt,'coefficient':target,
 'n_original_zero_growth_genes':len(essential),'eligible_non_moco_terms':eligible,'n_pairs':len(pairs),
 'rescue_threshold':'edited KO growth >= .95 * unedited WT growth',
 'mo_pair_ids':list(mo),'mo_pair_prior_rescued':sorted(mo_prior['rescued_ko_growth']),
 'pairs':[]}
if outpath.exists():
 priorout=json.load(open(outpath)); assert {k:priorout[k] for k in header}==header
 assert [x['removed'] for x in priorout['pairs']]==[list(t) for t in pairs[:len(priorout['pairs'])]]]
 out=priorout
else: out=header
start=len(out['pairs']); n=int(sys.argv[1]) if len(sys.argv)>1 else 25
for pair in pairs[start:start+n]:
 with m:
  ob=m.reactions.get_by_id(obj.id)
  ob.subtract_metabolites({m.metabolites.get_by_id(z):coef[z] for z in pair})
  patched=float(m.slim_optimize()); rescued={}; failures=[]; borderline=[]
  for gene in essential:
   with m:
    m.genes.get_by_id(gene).knock_out()
    val=m.slim_optimize()
    if val is None: failures.append(gene); continue
    growth=float(val)
    if abs(growth-.95*wt)<=1e-7*wt: borderline.append([gene,growth])
    if growth>=.95*wt: rescued[gene]=growth
  out['pairs'].append({'removed':list(pair),'edited_wt_growth':patched,'n_rescued':len(rescued),
   'rescued_ko_growth':rescued,'solver_failure_genes':failures,'borderline_genes':borderline,
   'overlap_with_mo_pair':sorted(set(rescued)&set(header['mo_pair_prior_rescued']))})
 outpath.write_text(json.dumps(out,indent=2)+'\n')
 print(f"{len(out['pairs'])}/{len(pairs)} {pair} rescued {len(rescued)} failures {len(failures)}",flush=True)
if len(out['pairs'])==len(pairs):
 counts=[r['n_rescued'] for r in out['pairs']]; target_count=len(header['mo_pair_prior_rescued'])
 out['summary']={'n_pairs_ge_mo_six':sum(x>=target_count for x in counts),
  'n_pairs_gt_mo_six':sum(x>target_count for x in counts),
  'zero_rescue_pairs':counts.count(0),'max_non_moco_rescues':max(counts),
  'any_solver_failure_pairs':sum(bool(r['solver_failure_genes']) for r in out['pairs']),
  'mo_gene_overlap_counts':{gene:sum(gene in r['rescued_ko_growth'] for r in out['pairs']) for gene in header['mo_pair_prior_rescued']},
  'interpretation':'descriptive enumeration after MoCo-pair outcome exposure; not a randomization p-value or biological validation'}
 outpath.write_text(json.dumps(out,indent=2)+'\n')
 print('SUMMARY',json.dumps(out['summary'],indent=2),flush=True)
