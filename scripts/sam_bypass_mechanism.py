"""Exploratory (not pre-registered): mechanism of SAM rescues of purine genes in iML1515 (M9 glucose, SBML default medium).
For each SAM-rescued off-pathway gene, test whether the rescue survives removal of candidate route reactions."""
import json, warnings
import cobra
warnings.filterwarnings('ignore')
m = cobra.io.read_sbml_model('external/E_coli_GEM_validation/Models/iML1515.xml')
wt = m.slim_optimize()
GENES = ['b0523', 'b0522', 'b1131', 'b4177', 'b2476']  # purK purE purB purA purC (results/rescue_audit_iML1515_M9.csv)
ROUTES = {'none': [], 'AHCYSNS (SAH nucleosidase)': ['AHCYSNS'], 'HCYSMT (homocysteine S-methyltransferase)': ['HCYSMT'],
          'ADD (adenine deaminase)': ['ADD'], 'ADPT (adenine phosphoribosyltransferase)': ['ADPT'], 'ADD+ADPT': ['ADD', 'ADPT']}
out = {'wt_growth': wt, 'rescued_growth': {}}
for g in GENES:
    out['rescued_growth'][g] = {}
    for name, rxns in ROUTES.items():
        with m:
            m.genes.get_by_id(g).knock_out()
            for r in rxns: m.reactions.get_by_id(r).knock_out()
            m.add_boundary(m.metabolites.amet_c, type='sink', reaction_id='SUPPLY_amet', lb=-10, ub=0)
            out['rescued_growth'][g][name] = round(m.slim_optimize(error_value=0.0) or 0.0, 4)
out['reaction_genes'] = {r: m.reactions.get_by_id(r).gene_reaction_rule for r in ['AHCYSNS', 'HCYSMT', 'ADD', 'ADPT']}
json.dump(out, open('results/sam_bypass_mechanism.json', 'w'), indent=1); print(json.dumps(out, indent=1))
