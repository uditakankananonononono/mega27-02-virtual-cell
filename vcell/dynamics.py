"""Dynamic FBA (dFBA): time-course simulation of the virtual cell.

Static-optimization approach (Mahadevan et al. 2002, Biophys J 83:1331):
at each step, extracellular concentrations set Michaelis-Menten uptake
bounds, FBA returns fluxes, and concentrations + biomass are integrated
forward explicitly. Reproduces the classic E. coli diauxic shift
(glucose first, acetate second) on the core model - a qualitative
experimental benchmark every virtual cell must pass.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def dynamic_fba(model, substrates: dict, biomass0: float = 0.01,
                dt: float = 0.1, t_end: float = 12.0,
                vmax: float = 10.0, km: float = 0.015,
                o2_bound: float = -1000.0,
                regulators: list[tuple[str, str, float]] | None = None) -> pd.DataFrame:
    """Run dFBA. substrates: {exchange_rxn_id: initial conc (mmol/L)}.
    Returns a trajectory frame: time, biomass, growth, per-substrate conc."""
    from .metabolism import carbon_exchange_ids, reset_medium
    carbon_ids = carbon_exchange_ids(model)
    conc = dict(substrates)
    t = 0.0
    X = biomass0
    rows = []
    while t <= t_end + 1e-9 and X > 0:
        reset_medium(model)
        for rid in carbon_ids:
            model.reactions.get_by_id(rid).lower_bound = 0.0
        uptake = {}
        for rid in conc:
            s = max(conc[rid], 0.0)
            uptake[rid] = -vmax * s / (km + s) if s > 1e-12 else 0.0
            model.reactions.get_by_id(rid).lower_bound = uptake[rid]
        model.reactions.get_by_id("EX_o2_e").lower_bound = o2_bound
        # Boolean regulation: (sensor_exchange, target_exchange, cap).
        # If the sensor substrate is still importable, cap the target uptake -
        # the classic catabolite-repression motif (glucose represses acetate).
        for sensor, target, cap in (regulators or []):
            if uptake.get(sensor, 0.0) < -1e-6:
                model.reactions.get_by_id(target).lower_bound = -cap
        sol = model.optimize()
        if sol.status != "optimal":
            mu = 0.0
            fluxes = {rid: 0.0 for rid in conc}
        else:
            mu = float(sol.objective_value)
            fluxes = {rid: float(sol.fluxes[rid]) for rid in conc}
        rows.append({"time": t, "biomass": X, "growth": mu, **{f"conc_{k}": v for k, v in conc.items()}})
        # explicit Euler integration (documented in paper; dt small)
        X = X * np.exp(mu * dt) if mu > 0 else X
        for rid in conc:
            conc[rid] = max(conc[rid] + fluxes[rid] * X * dt, 0.0)
        t += dt
    return pd.DataFrame(rows)


def diauxie_metrics(traj: pd.DataFrame, primary: str, secondary: str) -> dict:
    """Quantify the diauxic shift: time when primary is exhausted and when
    secondary consumption starts; presence of two growth phases."""
    cp, cs = f"conc_{primary}", f"conc_{secondary}"
    t_primary_out = float(traj.loc[traj[cp] <= 1e-6, "time"].min()) if (traj[cp] <= 1e-6).any() else float("nan")
    sub_before = traj.loc[traj["time"] < t_primary_out - 1e-9, cs] if t_primary_out == t_primary_out else traj[cs]
    co_utilization = bool((sub_before.iloc[-1] < traj[cs].iloc[0] * 0.5)) if len(sub_before) else False
    final_biomass = float(traj["biomass"].iloc[-1])
    return {"t_primary_exhausted_h": t_primary_out,
            "co_utilization_primary_phase": co_utilization,
            "final_biomass_gL": final_biomass,
            "peak_growth_h": float(traj["growth"].max()),
            "n_timepoints": int(len(traj))}
