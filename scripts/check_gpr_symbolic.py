"""Independent SymPy Boolean interpretation of every iJO1366 single-gene GPR knockout."""
import hashlib,json,re
from pathlib import Path
import cobra
from sympy import symbols
from sympy.logic.boolalg import And,Or,Not
from sympy.parsing.sympy_parser import parse_expr
src=Path('data/iJO1366.json');m=cobra.io.load_json_model(str(src));out=[];disagreements=[];n_tests=0;n_inactive=0
for r in m.reactions:
 rule=r.gene_reaction_rule
 if not rule:continue
 keys=sorted(g.id for g in r.genes)
 assert keys==sorted(g.id for g in r.genes),(r.id,keys)
 names={g:f'g{i}' for i,g in enumerate(keys)}
 expr=re.sub(r'(?<![A-Za-z0-9_.-])(?:'+'|'.join(re.escape(g) for g in sorted(keys,key=len,reverse=True))+r')(?![A-Za-z0-9_.-])',lambda x:names[x.group()],rule).replace(' and ',' & ').replace(' or ',' | ')
 syms={name:symbols(name, boolean=True) for name in names.values()}
 parsed=parse_expr(expr,local_dict=syms,global_dict={'And':And,'Or':Or,'Not':Not},evaluate=False)
 inactive=[]
 for g in keys:
  val=bool(parsed.subs({syms[names[k]]:k!=g for k in keys}))
  actual=bool(r.gpr.eval(knockouts={g}))
  n_tests+=1;n_inactive+=not val
  if val!=actual:disagreements.append({'reaction':r.id,'gene':g,'gpr':rule,'symbolic_active':val,'cobra_gpr_active':actual})
  if not val:inactive.append(g)
 out.append({'reaction':r.id,'genes':len(keys),'knockout_inactive':inactive})
moaD=sorted(x['reaction'] for x in out if 'b0784' in x['knockout_inactive'])
with m:
 m.genes.get_by_id('b0784').knock_out()
 disabled=sorted(x.id for x in m.genes.get_by_id('b0784').reactions if x.bounds==(0,0))
assert moaD==disabled,(moaD,disabled)
result={'design':'Independent symbolic interpretation of every non-empty iJO1366 GPR for all constituent single-gene knockouts','model_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'n_reactions_with_gpr':len(out),'n_reaction_gene_knockout_tests':n_tests,'n_inactive_calls':n_inactive,'n_disagreements':len(disagreements),'disagreements':disagreements,'moaD_symbolic_inactive_reactions':moaD,'moaD_actual_disabled_reactions':disabled,'per_reaction':out,'caveat':'Concordance tests model GPR rule semantics, not experimental essentiality or phenotype.'}
Path('results/gpr_symbolic.json').write_text(json.dumps(result,indent=2)+'\n');print({k:v for k,v in result.items() if k not in ('per_reaction','disagreements')})
