"""Rescue audit: flag supplement rescues whose gene is off the supplement's own pathway.

Finding behind this tool (paper H.7, results/pathway_concordance.json): on the
Bernstein 2023 E. coli RB-TnSeq benchmark, in-silico supplement rescues where the
knocked-out gene sits on the supplement's biosynthetic pathway have near-neutral
measured fitness (median -0.97), while off-pathway rescues are strongly deleterious
(median -4.24). Off-pathway rescues are therefore likely model artifacts: the
supplied metabolite lets FBA route around a gene the cell cannot actually bypass.

`rescue_audit` finds every gene whose knockout is lethal without a supplement and
viable with it, and labels each rescue `on_pathway` or `off_pathway` using a
gene -> pathway map (e.g. KEGG) and a supplement -> pathway map.
"""
from __future__ import annotations

import cobra


def _grow(model: cobra.Model) -> float:
    g = model.slim_optimize(error_value=0.0)
    return 0.0 if g is None or g != g else float(g)


def rescue_audit(model: cobra.Model, supplements: dict[str, set[str]],
                 gene_pathways: dict[str, set[str]], genes: list[str] | None = None,
                 tol: float = 0.01, supply_rate: float = 10.0) -> list[dict]:
    """Return one row per (gene, supplement) rescue.

    supplements: metabolite id (e.g. 'btn_c') -> set of pathway ids it belongs to.
    gene_pathways: gene id -> set of pathway ids.
    A rescue is a gene whose knockout grows < tol x wild type without the
    supplement and >= tol x wild type with an uptake sink for it.
    """
    wt = _grow(model)
    thr = tol * wt
    genes = genes if genes is not None else [g.id for g in model.genes]
    rows = []
    for gid in genes:
        with model:
            model.genes.get_by_id(gid).knock_out()
            base = _grow(model)
            if base >= thr:
                continue
            for met, pws in supplements.items():
                if met not in model.metabolites:
                    continue
                with model:
                    sink = model.add_boundary(model.metabolites.get_by_id(met), type="sink",
                                              reaction_id=f"SUPPLY_{met}", lb=-supply_rate, ub=0.0)
                    gr = _grow(model)
                if gr >= thr:
                    on = bool(gene_pathways.get(gid, set()) & set(pws))
                    rows.append({"gene": gid, "supplement": met, "growth_ko": base, "growth_rescued": gr,
                                 "label": "on_pathway" if on else "off_pathway",
                                 "flag": "" if on else "likely_artifact"})
    return rows


def load_kegg_links(path: str) -> dict[str, set[str]]:
    """Parse a KEGG REST `link/pathway/<org>` TSV into gene -> {pathway number}."""
    out: dict[str, set[str]] = {}
    for line in open(path):
        if "\t" not in line:
            continue
        g, p = line.rstrip("\n").split("\t")[:2]
        g = g.split(":", 1)[-1]
        p = p.split(":", 1)[-1]
        p = p[-5:] if p[-5:].isdigit() else p
        out.setdefault(g, set()).add(p)
    return out


def can_produce(model: cobra.Model, met_id: str, eps: float = 1e-6) -> bool:
    """True if the model (current bounds/medium) can make a positive net amount of met_id."""
    with model:
        dm = model.add_boundary(model.metabolites.get_by_id(met_id), type="demand", reaction_id=f"DMTEST_{met_id}")
        model.objective = dm
        v = model.slim_optimize(error_value=0.0)
    return v is not None and v == v and v > eps


def rescue_audit_production(model: cobra.Model, supplements: list[str], genes: list[str] | None = None,
                            tol: float = 0.01, supply_rate: float = 10.0) -> list[dict]:
    """KEGG-free variant: a rescue of gene g by supplement s is on_pathway when knocking out g
    abolishes de-novo production of s (g is needed to make s), else off_pathway (bypass)."""
    wt = _grow(model)
    thr = tol * wt
    genes = genes if genes is not None else [g.id for g in model.genes]
    rows = []
    for gid in genes:
        with model:
            model.genes.get_by_id(gid).knock_out()
            base = _grow(model)
            if base >= thr:
                continue
            for met in supplements:
                if met not in model.metabolites:
                    continue
                with model:
                    model.add_boundary(model.metabolites.get_by_id(met), type="sink",
                                       reaction_id=f"SUPPLY_{met}", lb=-supply_rate, ub=0.0)
                    gr = _grow(model)
                if gr >= thr:
                    on = not can_produce(model, met)
                    rows.append({"gene": gid, "supplement": met, "growth_ko": base, "growth_rescued": gr,
                                 "label": "on_pathway" if on else "off_pathway",
                                 "flag": "" if on else "likely_artifact"})
    return rows
