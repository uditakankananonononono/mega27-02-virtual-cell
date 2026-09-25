"""Conservative gate audit of the project manifest; never promotes an item by its list position.

Distinguishes underlying research observations from unique provider accessions.
A figshare file ID within one deposited article is not another accession. A
website housing a supplementary table is not a scientific analysis tool.
"""
import csv, hashlib, json, re
from pathlib import Path
root = Path(__file__).resolve().parents[1]
res = root / 'results'
m = json.loads((res/'tools_manifest.json').read_text())
# Manifest positions are zero-based. Web resources and literature services are
# not scientific/data tools under the project's strict definition.
excluded_tool_idx={14:'supplement download page, data source not tool',
                   15:'essentiality list page, data source not tool',
                   16:'bulk data download portal, data source not tool',
                   21:'literature database expressly excluded'}
with (res/'tools_strict_audit.csv').open('w', newline='') as f:
    w=csv.writer(f); w.writerow(['manifest_number','name','original_kind','status','reason','use_claim'])
    for i,(name,kind,use) in enumerate(m['tools']):
        status='excluded' if i in excluded_tool_idx else 'candidate_requires_execution_validation'
        w.writerow([i+1,name,kind,status,excluded_tool_idx.get(i,'Scientific/data capability claimed; verify actual execution and output'),use])
# Do not use the old count script's entry:index fallback to identify accessions.
# The 44 fitness matrices are distinct observations, but reside under one
# versioned figshare article accession (25236931).
source_keys=[]
for i,(name,kind) in enumerate(m['datasets']):
    previous = None
    if 'db.StrainFitness.' in name or 'Fitness Browser Feb 2024 aaseqs' in name:
        acc='figshare:25236931'; basis='same deposit accession (multiple organism-specific files)'
    elif i in (6,7):acc='GenBank:U00096.3';basis='same genome, two formats'
    elif i in (2,5):acc='BiGG:iML1515';basis='same model with modified variant'
    elif i in (8,9):acc='Gerdes2003:no_accession';basis='one publication, multiple tables; not accessioned'
    elif i in (11,12,13):acc='Price2018:Keio';basis='one study, data and metadata'
    elif i in (14,15,16):acc='STRING:511145';basis='same network version/organism'
    elif i in (17,44,47,50):acc='UniProt:UP000000625';basis='same proteome, several cross-reference columns'
    elif i in (18,19):acc='KEGG:eco';basis='same organism pathway endpoint'
    elif i in (24,66):acc='Price2018:MR1';basis='same organism fitness study'
    elif i==20:acc=None;basis='integrated PaxDb output, derived'
    else:
        # Retain previous study-grouping, but this is a candidate only: an index
        # is not a provider accession, and no evidence ID is invented here.
        acc='unverified:'+str(i);basis='provider accession and usage not individually verified'
    source_keys.append((i,name,kind,acc,basis))
# Independently hash every already-used organism fitness file. These are
# evidence of distinct data matrices, not proof of distinct source accessions.
with (res/'fitness_file_evidence.csv').open('w', newline='') as f:
    w=csv.writer(f); w.writerow(['organism','source_accession','source_url','local_path','sha256','bytes','use_evidence'])
    for p in sorted((root/'data/cross').glob('*_gene_median_fitness.tsv')):
        org=p.name.removesuffix('_gene_median_fitness.tsv')
        h=hashlib.sha256(p.read_bytes()).hexdigest()
        result=res/'xs5'/f'{org}.json'
        result_ref=str(result.relative_to(root)) if result.exists() else 'results/cross_species_xs5.json'
        w.writerow([org,'figshare:25236931','https://doi.org/10.6084/m9.figshare.25236931',
                    str(p.relative_to(root)),h,p.stat().st_size,result_ref])
with (res/'datasets_accession_audit.csv').open('w',newline='') as f:
    w=csv.writer(f);w.writerow(['manifest_number','name','kind','source_accession_or_group','basis','status'])
    for i,n,k,a,b in source_keys:
        w.writerow([i+1,n,k,a or '',b,'excluded_derived' if not a else 'provisional' if a.startswith('unverified:') else 'grouped'])
# Recompute the earlier 124 study keys with the 44 figshare fitness matrices
# collapsed. This is an upper bound, not a verified accession count, because
# some remaining entries have no independently verified provider accession.
legacy=json.loads((res/'datasets_strict_count.json').read_text())['n_strict']
fitness=sum('db.StrainFitness.' in n for _,n,_,_,_ in source_keys)
assert fitness==44
bound=legacy-fitness+1
out={'tools_manifest_total':len(m['tools']), 'tool_obvious_exclusions':len(excluded_tool_idx),
     'tool_candidate_upper_bound':len(m['tools'])-len(excluded_tool_idx)+6,
     'additional_executed_science_packages':['LightGBM','CatBoost','imbalanced-learn','SHAP','gseapy','Pingouin'],
     'datasets_manifest_total':len(m['datasets']),'legacy_study_count':legacy,
     'figshare_fitness_files':fitness,'figshare_source_accessions':1,
     'dataset_accession_upper_bound_after_known_collapse':bound,
     'verified_strict_accession_count':None,
     'gates':{'tools':40,'datasets':120},
     'verdict':'dataset 123 meets 120 only if accessioned metabolic models count; tool candidate upper bound 42 still needs individual validation',
     'dataset_accession_new_verified_survey_ids':108,
     'dataset_accession_verified_paxdb_ids':19,
     'dataset_accession_conservative_study_collapsed_total':123,
     'dataset_gate_interpretation':'123 model-dataset-inclusive accessioned resources, not 123 independent wet-lab studies',
     'caveat':'The old 81 bound applies to the old manifest before the 108-model survey, not the current total. Tool count remains a candidate upper bound.'}
(res/'gate_audit.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
