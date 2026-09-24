"""Essentiality-provenance auditor.

For every FBA-predicted-essential gene, decide WHY it is essential:
  - 'biomass_forced': the knockout only kills growth because the biomass
    objective requires compounds the cell can no longer produce. Removing
    exactly those compounds from the objective rescues growth.
  - 'network_forced': even with an unachievable biomass requirement removed,
    no flux routing restores growth - a genuine network lesion.

Each verdict is a re-runnable FBA experiment: falsifiable by construction.
"""
from __future__ import annotations

import cobra
import pandas as pd


def biomass_metabolites(model: cobra.Model) -> list[str]:
    """Consumed metabolites (negative coefficients) of the objective reaction."""
    from cobra.util.solver import linear_reaction_coefficients
    coeffs = linear_reaction_coefficients(model)
    rxn = next(iter(coeffs))
    return [m.id for m, c in rxn.metabolites.items() if c < 0], rxn.id


def unproducible_biomass_compounds(model: cobra.Model, biomass_mets: list[str]) -> list[str]:
    """Which biomass constituents can the (possibly knocked-out) model no longer make?
    Producibility = maximal flux through a temporary demand reaction > tiny."""
    out = []
    for mid in biomass_mets:
        with model:
            met = model.metabolites.get_by_id(mid)
            dem = cobra.Reaction(f"DM_audit_{mid}")
            dem.add_metabolites({met: -1.0})
            model.add_reactions([dem])
            model.objective = dem
            flux = model.slim_optimize(error_value=0.0)
        if flux is None or flux < 1e-6:
            out.append(mid)
    return out


def audit_gene(model: cobra.Model, gene_id: str, biomass_rxn_id: str,
               biomass_mets: list[str], growth_tol: float = 0.01) -> dict:
    """Classify one essential gene: biomass_forced vs network_forced, naming compounds."""
    wt = model.slim_optimize(error_value=0.0)
    with model:
        for rxn in model.genes.get_by_id(gene_id).reactions:
            rxn.knock_out()
        ko_growth = model.slim_optimize(error_value=0.0)
        if ko_growth is not None and ko_growth >= growth_tol * wt:
            return {"gene": gene_id, "verdict": "not_essential", "rescued_growth": ko_growth,
                    "unproducible": []}
        missing = unproducible_biomass_compounds(model, biomass_mets)
        if not missing:
            return {"gene": gene_id, "verdict": "network_forced_no_biomass_link",
                    "rescued_growth": ko_growth or 0.0, "unproducible": []}
        # rescue test: remove exactly the unproducible compounds from the objective
        bm = model.reactions.get_by_id(biomass_rxn_id)
        delta = {}
        for m in missing:
            met = model.metabolites.get_by_id(m)
            if met in bm.metabolites:
                delta[met] = -bm.metabolites[met]
        bm.add_metabolites(delta)
        rescued = model.slim_optimize(error_value=0.0)
    verdict = "biomass_forced" if (rescued is not None and rescued >= growth_tol * wt) else "network_forced"
    return {"gene": gene_id, "verdict": verdict, "rescued_growth": float(rescued or 0.0),
            "wt_growth": float(wt), "unproducible": missing}


def audit_model(model: cobra.Model, essential_genes: list[str]) -> pd.DataFrame:
    """Run the audit over a list of predicted-essential genes."""
    biomass_mets, bm_id = biomass_metabolites(model)
    rows = [audit_gene(model, g, bm_id, biomass_mets) for g in essential_genes]
    return pd.DataFrame(rows)
