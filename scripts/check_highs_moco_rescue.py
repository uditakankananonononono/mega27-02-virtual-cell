"""Independent HiGHS LP check of the exact MoCo biomass-edit rescue.

COBRApy reads and edits model/GPR; highspy constructs and solves the numeric
mass-balance LP directly, not via COBRApy's default GLPK interface. This
checks solver dependence, not biological validity of the model edit.
"""
import hashlib,json,math
from pathlib import Path
import cobra,highspy,numpy as np
from scipy.sparse import coo_matrix
from cobra.util.array import create_stoichiometric_matrix
SOURCE=Path('data/iJO1366.json')
model=cobra.io.load_json_model(str(SOURCE))
obj=model.reactions.get_by_id('BIOMASS_Ec_iJO1366_core_53p95M')
cof={x:v for x,v in obj.metabolites.items() if x.id in {'bmocogdp_c','mobd_c'}}
assert {x.id for x in cof}=={'bmocogdp_c','mobd_c'}

def solve(m):
 rx=list(m.reactions);met=list(m.metabolites);n=len(rx);k=len(met)
 c=np.array([-r.objective_coefficient for r in rx],np.float64)
 lb=np.array([r.lower_bound for r in rx],np.float64);ub=np.array([r.upper_bound for r in rx],np.float64)
 s=create_stoichiometric_matrix(m,array_type='lil').tocsr().tocsc()
 lp=highspy.Highs();lp.setOptionValue('output_flag',False);lp.setOptionValue('primal_feasibility_tolerance',1e-8)
 zero=np.zeros(k,np.float64);st=np.zeros(k+1,np.int32)
 assert lp.addRows(k,zero,zero,0,st,np.array([],np.int32),np.array([],np.float64))==highspy.HighsStatus.kOk
 a=lp.addCols(n,c,lb,ub,s.nnz,s.indptr.astype(np.int32),s.indices.astype(np.int32),s.data.astype(np.float64))
 assert a==highspy.HighsStatus.kOk
 assert lp.run()==highspy.HighsStatus.kOk
 status=lp.getModelStatus();assert status==highspy.HighsModelStatus.kOptimal,status
 flux=lp.getSolution().col_value;v=np.array(flux)
 imbalance=np.max(np.abs(s.tocsr()@v));assert imbalance<1e-6,imbalance
 value=float(v[rx.index(m.reactions.get_by_id('BIOMASS_Ec_iJO1366_core_53p95M'))])
 return {'growth_per_h':value,'max_mass_balance_residual':float(imbalance),'n_reactions':n,'n_metabolites':k,'nonzero_stoichiometry':int(s.nnz)}

out={'design':'Independent LP solver and direct S*v=0 verification of the previously reported MoCo biomass dependency',
     'model':str(SOURCE),'sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
     'objective_reaction':obj.id,'removed_biomass_coefficients':{x.id:v for x,v in cof.items()},
     'knockout_gene':'b0784 (moaD)','solver':'HiGHS '+highspy.Highs().version()}
out['wild_type']=solve(model)
with model:
 model.genes.get_by_id('b0784').knock_out();out['moaD_knockout_original']=solve(model)
with model:
 obj.subtract_metabolites(cof);out['biomass_edit_wild_type']=solve(model)
 model.genes.get_by_id('b0784').knock_out();out['biomass_edit_moaD_knockout']=solve(model)
assert abs(out['wild_type']['growth_per_h']-.9823718127269785)<1e-6
assert abs(out['moaD_knockout_original']['growth_per_h'])<1e-7
assert abs(out['biomass_edit_wild_type']['growth_per_h']-.9824851916537897)<1e-6
assert abs(out['biomass_edit_moaD_knockout']['growth_per_h']-out['biomass_edit_wild_type']['growth_per_h'])<1e-6
out['caveat']='Solver concordance validates the saved model computation, not a live-cell phenotype or universal biomass repair.'
Path('results/highs_moco_rescue.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
