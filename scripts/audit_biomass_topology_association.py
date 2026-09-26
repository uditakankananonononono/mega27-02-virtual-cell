"""Post-result descriptive feature/flip association with permutation negative control.
No fitted predictor, independent test, or prospective hypothesis claim.
"""
import json
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr
m=json.load(open('results/biomass_topology_manifest.json'));x=json.load(open('results/biomass_single_term_sensitivity.json'))
lookup={z['metabolite_id']:z for z in x['rows']}; assert set(lookup)=={z['metabolite_id'] for z in m['terms']}
y=np.array([lookup[z['metabolite_id']]['n_baseline_essential_rescued'] for z in m['terms']])
keys=['abs_biomass_coefficient','non_biomass_reaction_count','gene_associated_reaction_count','active_producer_reaction_count_wt','active_producer_abs_flux_sum_wt']
rng=np.random.default_rng(20260926);feat=[]
for k in keys:
 v=np.array([z[k] for z in m['terms']]); rho=float(spearmanr(v,y).statistic);null=np.array([spearmanr(rng.permutation(v),y).statistic for _ in range(2000)])
 feat.append({'feature':k,'spearman_rho':rho,'random_feature_permutation_p_two_sided':float((1+(np.abs(null)>=abs(rho)).sum())/(len(null)+1)),'unique_values':len(np.unique(v))})
out={'status':'exploratory in-sample association only; features and 48 targets known before full 102-term scoring',
 'n_terms':len(y),'nonzero_flip_terms':int((y>0).sum()),'zero_flip_terms':int((y==0).sum()),
 'median_flips_all':float(np.median(y)),'max_flips':int(y.max()),'moco_individual_flips':{z['metabolite_id']:lookup[z['metabolite_id']]['n_baseline_essential_rescued'] for z in m['terms'] if z['outcome_status']=='previously_scored_moco_bundle_only'},
 'feature_correlations':feat,'permutation_draws_each':2000,'permutation_seed':20260926,
 'limits':'feature selection and target exposure post-result; correlated reaction features, multiple testing, zeros and one model; permutation checks numeric association but does not validate prediction or biological mechanism'}
Path('results/biomass_topology_exploratory.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
