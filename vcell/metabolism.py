"""Genome-scale constraint-based metabolism layer (FBA) of the virtual cell.

Uses COBRApy with real BiGG models (e_coli_core, iJO1366). Every number
produced here comes from an actual LP solve on a published reconstruction.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import cobra
import pandas as pd
from cobra.flux_analysis import single_gene_deletion

# Inorganic / non-carbon exchange metabolites that stay open in every condition.
_NON_CARBON_EXCHANGE = {"h2o", "h", "co2", "nh4", "pi", "o2", "so4", "k", "na1", "mg2", "fe2", "fe3"}


def load_model(path: str) -> cobra.Model:
    """Load a BiGG JSON model and snapshot its default medium."""
    model = cobra.io.load_json_model(path)
    model._vcell_default_medium = dict(model.medium)  # exchange_id -> lower bound
    return model


def reset_medium(model: cobra.Model) -> None:
    """Restore the medium snapshot taken at load time."""
    # cobra's model.medium stores uptake rates as POSITIVE numbers that map to
    # NEGATIVE lower bounds; flip the sign back when restoring.
    default = model._vcell_default_medium
    for rxn in model.exchanges:
        rxn.lower_bound = -default.get(rxn.id, 0.0)


def carbon_exchange_ids(model: cobra.Model) -> list[str]:
    """Exchange reactions whose metabolite formula contains carbon."""
    ids = []
    for rxn in model.exchanges:
        mets = list(rxn.metabolites)
        if len(mets) != 1:
            continue
        formula = (mets[0].formula or "")
        base = mets[0].id.rsplit("_", 1)[0]  # strip compartment suffix
        if "C" in formula and base not in _NON_CARBON_EXCHANGE:
            ids.append(rxn.id)
    return sorted(ids)


@dataclass
class Condition:
    """A growth condition: one carbon source + aerobiosis."""
    carbon_exchange: str
    aerobic: bool = True
    uptake_bound: float = -10.0
    name: str = field(default="")

    def __post_init__(self):
        if not self.name:
            self.name = f"{self.carbon_exchange}_{'aerobic' if self.aerobic else 'anaerobic'}"


def apply_condition(model: cobra.Model, cond: Condition, carbon_ids: list[str]) -> None:
    """Restore defaults, close every carbon exchange, open the chosen one, set O2."""
    reset_medium(model)
    for rid in carbon_ids:
        model.reactions.get_by_id(rid).lower_bound = 0.0
    model.reactions.get_by_id(cond.carbon_exchange).lower_bound = cond.uptake_bound
    o2 = model.reactions.get_by_id("EX_o2_e")
    o2.lower_bound = -1000.0 if cond.aerobic else 0.0


def growth_rate(model: cobra.Model, cond: Condition, carbon_ids: list[str]) -> float:
    """Maximise biomass under a condition; returns h^-1 (0 if infeasible)."""
    apply_condition(model, cond, carbon_ids)
    sol = model.optimize()
    return float(sol.objective_value) if sol.status == "optimal" else 0.0


def default_condition_panel(model: cobra.Model) -> list[Condition]:
    """Aerobic + anaerobic panel over every carbon source the model can import."""
    panel = []
    for rid in carbon_exchange_ids(model):
        panel.append(Condition(rid, aerobic=True))
        panel.append(Condition(rid, aerobic=False))
    return panel


def condition_matrix(model: cobra.Model, panel: list[Condition]) -> pd.DataFrame:
    """Growth rate for every condition in the panel."""
    carbon_ids = carbon_exchange_ids(model)
    rows = []
    for cond in panel:
        gr = growth_rate(model, cond, carbon_ids)
        rows.append({"condition": cond.name, "carbon": cond.carbon_exchange,
                     "aerobic": cond.aerobic, "growth_h": gr, "grows": gr > 1e-6})
    return pd.DataFrame(rows)


def gene_essentiality_scan(model: cobra.Model, fraction_of_optimum: float = 0.01,
                           processes: int = 2) -> pd.DataFrame:
    """In-silico single-gene deletion screen.

    A gene is predicted essential if its deletion drops maximal growth below
    `fraction_of_optimum` of the wild-type optimum (the standard FBA
    essentiality criterion used by the iJO1366/iML1515 papers).
    """
    reset_medium(model)
    wt = model.slim_optimize()
    res = single_gene_deletion(model, processes=processes)
    res = res.copy()
    res["growth"] = res["growth"].fillna(0.0)  # infeasible knockout = no growth
    res["gene"] = [next(iter(ids)) for ids in res["ids"]]
    res["growth_frac"] = res["growth"] / wt if wt and wt > 0 else 0.0
    res["predicted_essential"] = res["growth_frac"] < fraction_of_optimum
    res["wt_growth_h"] = wt
    return res[["gene", "growth", "growth_frac", "predicted_essential", "wt_growth_h", "status"]]
