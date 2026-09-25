"""Evidence-oriented tool gate: named source code plus committed scientific output.

An executed dependency or a saved dataset alone does not certify a tool. Thin
provenance stays thin; adding a new library to a manifest is not evidence.
"""
import csv, json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
old={
1:('scripts/run_bernstein_benchmark.py','results/bernstein/iML1515_base.json','COBRApy FBA simulation'),
4:('scripts/check_sbml_stoichiometry.py','results/sbml_route_verification.json','direct libSBML FBC stoichiometry for four model reactions'),
5:('scripts/run_seqcnn.py','results/seqcnn_results.json','PyTorch sequence CNN via vcell/seqcnn.py'),
6:('scripts/run_ensemble_v2.py','results/ensemble_v2.json','NumPy score/array computation'),
7:('scripts/run_ensemble_v2.py','results/ensemble_v2.json','pandas joins and scores'),
8:('scripts/run_network_modules.py','results/network_modules.json','SciPy Fisher tests'),
9:('scripts/run_ensemble_v2.py','results/ensemble_v2.json','scikit-learn OOF estimators'),
10:('vcell/external_features.py','results/external_features.csv','networkx STRING centralities'),
11:('scripts/run_network_modules.py','results/network_modules.json','statsmodels BH FDR'),
12:('scripts/esm_embed.py','results/esm2_t6_ids.json','Biopython GenBank CDS extraction'),
13:('scripts/run_bigg_model_survey.py','results/bigg_model_survey.csv','BiGG API; 108 fetched, hashed, parsed model IDs'),
14:('scripts/verify_source_reactions.py','results/source_reaction_verification.json','NCBI E-utilities identity against actual GenBank record'),
18:('scripts/run_bernstein_benchmark.py','results/bernstein/iML1515_base.json','executes released Bernstein notebook pipeline'),
19:('vcell/external_features.py','results/external_features.csv','STRING network into model features'),
20:('vcell/external_features.py','results/external_features.csv','UniProt protein features from saved query table'),
21:('scripts/run_pathway_enrichment.py','results/pathway_enrichment.json','KEGG pathways tested by gseapy'),
23:('scripts/run_paxdb_datasets.py','results/paxdb_accession_evidence.csv','PaxDb accessioned abundance scored per dataset'),
24:('scripts/esm_embed.py','results/ensemble_v4_esm.json','ESM-2 embedding and v4 test'),
25:('scripts/run_go_error_enrichment.py','results/go_error_enrichment.json','goatools GO overrepresentation'),
26:('results/memote_iML1515.html','results/memote_iML1515_summary.json','MEMOTE embedded report with score and tests'),
28:('scripts/run_xs5_org.py','results/cross_species_carveme.json','CarveMe subprocess reconstructed SBML models tested across organisms'),
31:('scripts/verify_source_reactions.py','results/source_reaction_verification.json','Rhea API reaction aligned to AHCYSNS knockout'),
32:('scripts/run_xgb_shap.py','results/xgb_shap.json','XGBoost OOF scores and TreeSHAP'),
33:('scripts/fetch_alphafold.py','results/v6_structure_domain.json','AlphaFold API pLDDT tested in v6'),
34:('scripts/check_sah_orthogonal_annotations.py','results/sah_orthogonal_annotations.json','Pfam/InterPro live API cross-check against v6 UniProt xref feature'),
35:('scripts/fetch_oma.py','results/v7_oma.json','OMA HOG levels tested in v7'),
36:('scripts/run_v8_cog_pdb.py','results/v8_cog_pdb.json','NCBI COG category/spread tested in v8'),
37:('scripts/check_sah_orthogonal_annotations.py','results/sah_orthogonal_annotations.json','RCSB PDB live entry cross-check against v8 UniProt xref feature'),
38:('scripts/run_v9_precise1k.py','results/v9_precise1k.json','PRECISE-1K/iModulon expression tested in v9'),
39:('scripts/run_uptake_plausibility.py','results/uptake_plausibility.json','TCDB substrate joins on ChEBI for seven supplements'),
40:('scripts/run_uptake_plausibility.py','results/uptake_plausibility.json','ChEBI ontology API expansion for seven supplements')}
excluded={2:'indirect COBRApy LP dependency; no distinct execution/result tie',3:'indirect COBRApy solver dependency; no distinct execution/result tie',15:'supplement source page only',16:'source page only',17:'source download portal only',22:'literature database expressly excluded',27:'PSAMM conversion is asserted in a note; no committed execution log',29:'DIAMOND runs within CarveMe but no distinct execution trace',30:'SCIP solver selected in CarveMe command; no distinct result-specific trace'}
new=[('PANTHER','scripts/check_sah_orthogonal_annotations.py','results/sah_orthogonal_annotations.json'),('QuickGO','scripts/check_sah_orthogonal_annotations.py','results/sah_orthogonal_annotations.json'),('LightGBM','scripts/run_learner_sensitivity.py','results/learner_sensitivity.json'),('CatBoost','scripts/run_learner_sensitivity.py','results/learner_sensitivity.json'),('imbalanced-learn','scripts/run_learner_sensitivity.py','results/learner_sensitivity.json'),('SHAP','scripts/run_learner_sensitivity.py','results/learner_sensitivity.json'),('gseapy','scripts/run_pathway_enrichment.py','results/pathway_enrichment.json'),('Pingouin','scripts/run_score_agreement.py','results/score_agreement.json'),('UMAP','scripts/run_error_structure.py','results/error_structure.json'),('HDBSCAN','scripts/run_error_structure.py','results/error_structure.json'),('PyOD','scripts/run_error_structure.py','results/error_structure.json'),('python-igraph','scripts/run_network_modules.py','results/network_modules.json'),('Leidenalg','scripts/run_network_modules.py','results/network_modules.json')]
manifest=json.loads((R/'results/tools_manifest.json').read_text())['tools']
assert len(manifest)==40 and set(old)|set(excluded)==set(range(1,41)) and not(set(old)&set(excluded))
rows=[]
for i,(name,kind,claim) in enumerate(manifest,1):
 if i in old:
  code,result,basis=old[i];assert (R/code).is_file() and (R/result).is_file() and (R/result).stat().st_size>0,(i,result)
  status='code_result_linked'
 else: code=result='';basis=excluded[i];status='excluded_or_unverified'
 rows.append({'tool':name,'manifest_number':i,'kind':kind,'status':status,'code_or_report':code,'result':result,'basis':basis})
for name,code,result in new:
 assert (R/code).is_file() and (R/result).is_file() and (R/result).stat().st_size>0
 rows.append({'tool':name,'manifest_number':'new','kind':'scientific package','status':'code_result_linked','code_or_report':code,'result':result,'basis':'Executed statistical/modeling operation visible in code and numerical committed output'})
with (R/'results/tool_execution_evidence.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
summary={'original_manifest':40,'old_code_result_linked':len(old),'new_code_result_linked':len(new),
 'total_code_result_linked':len(old)+len(new),'excluded_or_unverified':len(excluded),
 'audit_method':'Specific code/result links inspected for 31 original and 13 added rows. Some saved-source rows lack acquisition logs; source identity and independent reproduction still require audit. Dependency-only, build and literature tools excluded.',
 'interpretation':'44 code/result-linked candidates if distinct scientific/data databases and xrefs count as tools. Four-row margin over 40; five strict exclusions reduce it below 40. This is not a certified pass.',
 'sensitive_rows':['UniProt/Pfam/PDB used shared UniProt exports in v6/v8; separate Pfam and RCSB API calls now corroborate one protein, but an annotation check is thinner than a new model study.','MEMOTE embedded report is intact but original execution log is not retained.','NCBI and Rhea identity/reaction checks are useful research checks, not new biological datasets or independent discovery.'],
 'ledger':'results/tool_execution_evidence.csv'}
(R/'results/tool_execution_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
