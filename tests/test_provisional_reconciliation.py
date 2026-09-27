import csv,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]/'results'
def test_reconciliation_counts_and_uniqueness():
    s=json.loads((ROOT/'provisional_accession_reconciliation.json').read_text())
    rows=list(csv.DictReader((ROOT/'provisional_accession_reconciliation.csv').open()))
    audit=list(csv.DictReader((ROOT/'datasets_accession_audit.csv').open()))
    prov=[r for r in audit if r['status']=='provisional']
    assert len(prov)==s['n_provisional_before']==88
    assert len(rows)==s['n_resolved']==80
    assert s['n_still_provisional']==8
    ids=[r['manifest_number'] for r in rows];assert len(ids)==len(set(ids))
    # Multiple label sets reuse one Poulsen 2019 supplemental DOI; a DOI is not a study count.
    accession=[r['accession'] for r in rows]
    assert accession.count('DOI:10.1073/pnas.1900570116')==5
    assert len(set(accession))==76
    assert s['extension_pass']['n_added']==15
    assert all(len(r['sha256'])==64 for r in rows)
    assert all(r['source_url'].startswith('https://') for r in rows)
    assert all(r['use_evidence'] for r in rows)
    assert s['resolved_by_arm']=={'208964':10,'224308':5,'83332':14,'bigg':5,'cross':9,'ecoli':19,'kegg':2,'tcdb':1}
    assert sum(r['arm']=='extension' for r in rows)==15
    kv=json.loads((ROOT/'kegg_map_source_verification.json').read_text())['per_code']
    assert kv['sme']['exact_bytes'] is False
    assert sum(1 for v in kv.values() if v['exact_bytes'])==7
