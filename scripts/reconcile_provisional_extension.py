"""Second-pass provisional-accession overlay: essentiality labels, structure, orthology, compendia.

Evidence rule (same as first pass): local file exists, sha256 computed here, and scored/scripted
use in this repo is cited. Rows failing any check stay provisional. Overlay only; audit untouched.
"""
import csv, hashlib, json
from pathlib import Path
R = Path(__file__).resolve().parents[1]
res = R / 'results'
audit = list(csv.DictReader((res / 'datasets_accession_audit.csv').open()))
done = {r['manifest_number'] for r in csv.DictReader((res / 'provisional_accession_reconciliation.csv').open())}
still = {r['manifest_number']: r for r in audit if r['status'] == 'provisional' and r['manifest_number'] not in done}

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def row(num, accession, url, path, evidence):
    r = still.get(num)
    assert r, num
    p = R / path
    if not p.exists(): return None
    return dict(manifest_number=num, name=r['name'], kind=r['kind'], accession=accession,
                source_url=url, sha256=sha(p), use_evidence=evidence, arm='extension')

cands = [
 row('23', 'EcoCyc:83333-GAF', 'https://ecocyc.org/ (Gene Ontology annotations, taxon 83333)', 'data/external/ecocyc.gaf.gz', 'used by scripts/run_go_error_enrichment.py and scripts/run_ensemble_v5_go.py'),
 row('24', 'MEMOTE-model:iMR1_799', 'https://github.com/zhanglab/psamm-model-collection (Pinchuk 2010)', 'data/cross/iMR1_799.xml', 'used in results/cross_species_mr1.json via scripts/run_cross_species_mr1.py'),
 row('43', 'DOI:10.1371/journal.pgen.1007749.s012', 'https://doi.org/10.1371/journal.pgen.1007749.s012', 'data/external/rousset2018/pgen.1007749.s012.csv', 'used in results/rousset_validation.json'),
 row('44', 'AlphaFoldDB:AF-UP-F1-set', 'https://alphafold.ebi.ac.uk/ (AF-<UniProt>-F1 per-gene predictions)', 'data/external/alphafold/af_metrics.tsv', 'used by scripts/run_v6_structure_domain.py'),
 row('46', 'OMA:ECOLI', 'https://omabrowser.org/ (E. coli K-12 MG1655 genome + HOG levels)', 'data/external/oma/oma_hog_levels.tsv', 'used by scripts/run_v7_oma.py'),
 row('47', 'NCBI:COG-2020', 'https://www.ncbi.nlm.nih.gov/research/cog-project/ (cog-20 release)', 'data/external/cog/cog_genome_spread.tsv', 'used by scripts/run_v8_cog_pdb.py'),
 row('49', 'GitHub:SBRG/precise1k@' + (R/'data/external/precise1k/SOURCE_COMMIT.txt').read_text().strip().split()[0], 'https://github.com/SBRG/precise1k', 'data/external/precise1k/M.csv', 'used by scripts/run_v9_precise1k.py'),
 row('131', 'DOI:10.1128/mBio.02133-16', 'https://doi.org/10.1128/mBio.02133-16 (Table S3)', 'data/external/mtb/DeJesus_mbio.xlsx', 'used in results/paxdb_mtb.json'),
 row('137', 'SubtiWiki:essential', 'https://subtiwiki.uni-goettingen.de/ (Essential genes page, Koo et al. 2017 based)', 'data/external/bsub/subtiwiki_essential_genes.html', 'used in results/bsub_rescue_subtiwiki_overlap.json and scripts/run_paxdb_bsub.py'),
 row('148', 'DOI:10.1073/pnas.1900570116', 'https://doi.org/10.1073/pnas.1900570116 (Dataset S5)', 'data/external/pao1/pnas.1900570116.sd05.xlsx', 'used by scripts/run_paxdb_pao1.py'),
 row('149', 'DOI:10.1371/journal.ppat.1002251', 'https://doi.org/10.1371/journal.ppat.1002251 (Table 2)', 'data/external/mtb/griffin2011_table2.xlsx', 'used in results/paxdb_mtb_griffin.json'),
 row('150', 'DOI:10.1073/pnas.1900570116', 'https://doi.org/10.1073/pnas.1900570116 (Dataset S6, Turner 2015 calls)', 'data/external/pao1/pnas.1900570116.sd06.xlsx', 'used by scripts/run_paxdb_pao1_labels.py'),
 row('151', 'DOI:10.1073/pnas.1900570116', 'https://doi.org/10.1073/pnas.1900570116 (Dataset S6, Lee 2015 calls)', 'data/external/pao1/pnas.1900570116.sd06.xlsx', 'used by scripts/run_paxdb_pao1_labels.py'),
 row('152', 'DOI:10.1073/pnas.1900570116', 'https://doi.org/10.1073/pnas.1900570116 (Dataset S6, Skurnik 2013 calls)', 'data/external/pao1/pnas.1900570116.sd06.xlsx', 'used by scripts/run_paxdb_pao1_labels.py'),
 row('153', 'DOI:10.1073/pnas.1900570116', 'https://doi.org/10.1073/pnas.1900570116 (Dataset S6, Liberati 2006 calls)', 'data/external/pao1/pnas.1900570116.sd06.xlsx', 'used by scripts/run_paxdb_pao1_labels.py'),
]
new = [c for c in cands if c]
with (res / 'provisional_accession_reconciliation.csv').open('a', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['manifest_number','name','kind','accession','source_url','sha256','use_evidence','arm'])
    w.writerows(new)
summary = json.loads((res / 'provisional_accession_reconciliation.json').read_text())
summary['n_resolved'] += len(new)
summary['n_still_provisional'] -= len(new)
summary['extension_pass'] = {'n_added': len(new), 'still_provisional_numbers': sorted(still.keys() - {c['manifest_number'] for c in new}),
 'note': 'EcoliWiki list and go-basic.obo have no local hashed file; unused KEGG codes ccs/evi/dsu/pact/sme/psb have no scored use in-repo. These stay provisional.'}
(res / 'provisional_accession_reconciliation.json').write_text(json.dumps(summary, indent=2) + '\n')
print(json.dumps({'added': len(new), 'still': summary['n_still_provisional']}, indent=1))
