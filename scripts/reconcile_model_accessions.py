"""Add verified BiGG model accession analysis to conservative gate accounting.

The result is a lower bound built only from distinct model IDs with source
snapshots, hashes, and computed rows. Other legacy data remain provisional.
"""
import csv,json
from pathlib import Path
R=Path(__file__).resolve().parents[1]/'results'
rows=list(csv.DictReader((R/'bigg_model_survey.csv').open()))
used=[r for r in rows if r['status']=='analysed' and len(r['sha256'])==64 and r['source_url'].endswith('/'+r['accession']+'.json')]
assert len(used)==len({r['accession'] for r in used})==108
# Previously used unique model IDs shown explicitly in the project manifest.
manifest=json.loads((R/'tools_manifest.json').read_text())['datasets']
prior={'e_coli_core','iJO1366','iML1515','iJR904','iAF1260','iJN1463'}
for id in prior:
    assert any(id in name for name,_ in manifest),id
new=[r for r in used if r['accession'] not in prior]
with (R/'bigg_accession_delta.csv').open('w',newline='') as f:
    w=csv.writer(f);w.writerow(['accession','source_url','sha256','bytes','organism','moco_biomass','status'])
    for r in rows:w.writerow([r['accession'],r['source_url'],r['sha256'],r['bytes'],r['organism'],r['moco_biomass'],'previously used' if r['accession'] in prior else 'new distinct model ID'])
out={'survey_analysed_unique_model_ids':len(used),'previous_manifest_overlap_explicit_ids':sorted(prior),
     'new_distinct_bigg_model_ids_after_overlap':len(new),'dataset_gate_target':120,
     'provider_model_accessions_are_datasets':True,
     'independent_wet_lab_datasets_from_this_survey':0,
     'combined_gate_verdict':'Model-dataset-inclusive criterion: 108 distinct BiGG model IDs plus at least 15 PaxDb accessioned/scored E. coli source studies = 123, over 120. Not 123 independent wet-lab screens. Tool gate remains unverified.',
     'note':'Do not add the old 124 study keys; many lacked distinct accessions. The 102 newly used BiGG IDs are individually hashed in bigg_model_survey.csv.',
     'verified_distinct_provider_ids':127,
     'conservative_gate_floor_after_technical_variant_collapse':123,
     'technical_variant_collapse_note':'PaxDb 19 accession files include three Arike 2012 and two Krug 2013 quantification variants; collapse these to 15 source studies, each with accessioned and scored files.'}
(R/'bigg_accession_delta.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
