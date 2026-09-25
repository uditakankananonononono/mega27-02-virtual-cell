"""Direct SBML/FBC checks of the saved iML1515 SAM-to-adenine route.

libSBML is used as a parser and biochemical-model inspection tool here; it
is not counted solely for being installed beneath COBRApy.
"""
import json,hashlib
from pathlib import Path
import libsbml
p=Path('external/E_coli_GEM_validation/Models/iML1515.xml')
doc=libsbml.readSBMLFromFile(str(p));assert doc.getNumErrors()==0
m=doc.getModel();assert m and m.getId()=='iML1515'
want=['AHCYSNS','HCYSMT','ADD','ADPT']
out={}
for rid in want:
 x=m.getReaction('R_'+rid);assert x is not None,rid
 a={r.getSpecies().removeprefix('M_'):-float(r.getStoichiometry()) for r in x.getListOfReactants()}
 b={r.getSpecies().removeprefix('M_'):float(r.getStoichiometry()) for r in x.getListOfProducts()}
 out[rid]={'reactants':a,'products':b,'reversible':bool(x.getReversible())}
assert out['AHCYSNS']['reactants']=={'ahcys_c':-1.0,'h2o_c':-1.0}
assert out['AHCYSNS']['products']=={'ade_c':1.0,'rhcys_c':1.0}
assert out['HCYSMT']['reactants'].get('amet_c')==-1
assert out['ADPT']['products'].get('amp_c')==1
# Element balance is a model assertion, not measured flux.
elements={}
for sid in ['M_ahcys_c','M_h2o_c','M_ade_c','M_rhcys_c']:
 sp=m.getSpecies(sid);assert sp is not None
 plugin=sp.getPlugin('fbc');assert plugin is not None
 elements[sid.removeprefix('M_')]=plugin.getChemicalFormula()
r=json.loads(Path('results/source_reaction_verification.json').read_text())['rhea_rest']
assert r['reaction_id']=='RHEA:17805' and r['model_stoichiometry']==dict(out['AHCYSNS']['reactants'],**out['AHCYSNS']['products'])
rec={'scope':'Direct independent SBML parse of saved iML1515 cofactor-to-purine route; no new phenotype or wet-lab measurement.',
     'libsbml_version':libsbml.getLibSBMLDottedVersion(),'model_path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),
     'model_id':m.getId(),'n_reactions':m.getNumReactions(),'route':out,'fbc_formulas':elements,
     'rhea_check':'results/source_reaction_verification.json',
     'caveat':'Correct SBML stoichiometry and Rhea reaction identity do not show that external SAM enters E. coli.'}
Path('results/sbml_route_verification.json').write_text(json.dumps(rec,indent=2)+'\n')
print(json.dumps(rec,indent=2))
