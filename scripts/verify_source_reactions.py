"""Independent source checks for genome accession and the SAM-bypass reaction.

These API calls are scientific checks, not an added wet-lab dataset or validation
of essentiality. Store exact source responses and their relationship to the
saved analysis so the services are traceable research tools.
"""
import hashlib, json, re, urllib.request
from pathlib import Path
from Bio import SeqIO
import cobra

HEAD = {'User-Agent': 'VC2-virtual-cell/1.0 (research provenance check)'}
def get(url):
    with urllib.request.urlopen(urllib.request.Request(url, headers=HEAD), timeout=30) as r:
        return r.read().decode('utf-8')

ncbi = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=U00096.3&rettype=acc&retmode=text'
accession_response = get(ncbi).strip()
genome = SeqIO.read('data/U00096.3.gb', 'genbank')
assert accession_response == genome.id == 'U00096.3'
assert len(genome.seq) == 4641652

rhea = 'https://www.rhea-db.org/rhea/?query=id:17805&columns=rhea-id,equation,ec&format=tsv&limit=5'
rows = get(rhea).strip().splitlines()
assert rows[0].split('\t') == ['Reaction identifier', 'Equation', 'EC number']
hits = [r.split('\t') for r in rows[1:] if r.startswith('RHEA:17805\t')]
assert len(hits) == 1
rid, equation, ec = hits[0]
model = cobra.io.read_sbml_model('external/E_coli_GEM_validation/Models/iML1515.xml')
rx = model.reactions.get_by_id('AHCYSNS')
name_by_id = {m.id: m.name.lower() for m in rx.metabolites}
assert any('adenine' in name_by_id[m.id] and v > 0 for m, v in rx.metabolites.items())
assert any('adenosyl-l-homocysteine' in name_by_id[m.id] and v < 0 for m, v in rx.metabolites.items())
assert ec == 'EC:3.2.2.9'
knockout = json.loads(Path('results/sam_bypass_mechanism.json').read_text())['rescued_growth']
blocked = sum(float(x['AHCYSNS (SAH nucleosidase)']) <= 0.01 for x in knockout.values())
assert blocked == len(knockout) == 5
out = {
    'scope':'Two research API checks: record identity and reaction mechanism. Neither is independent experimental replication.',
    'ncbi_entrez': {'url':ncbi,'response':accession_response,'local_genbank':'data/U00096.3.gb',
      'local_genbank_sha256':hashlib.sha256(Path('data/U00096.3.gb').read_bytes()).hexdigest(),
      'sequence_length':len(genome.seq),'n_cds':sum(f.type=='CDS' for f in genome.features)},
    'rhea_rest': {'url':rhea,'reaction_id':rid,'equation':equation,'ec':ec,
      'iML1515_model_reaction':'AHCYSNS','model_equation':rx.reaction,
      'model_stoichiometry':{m.id:v for m,v in rx.metabolites.items()},
      'knockout_source':'results/sam_bypass_mechanism.json',
      'sam_rescued_purine_genes_blocked_when_AHCYSNS_removed':blocked,
      'caveat':'Rhea confirms chemical annotation and model alignment, not in-vivo transport or a validated rescue.'}
}
Path('results/source_reaction_verification.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
