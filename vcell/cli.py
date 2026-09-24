"""Command-line entry: `python -m vcell audit-biomass MODEL [--genes ...] [--out CSV]`.

Works on any COBRA model file (BiGG JSON or SBML). If --genes is omitted, the
tool first runs a single-gene deletion screen and audits every gene whose
knockout drops growth below `--tol` x wild type.
"""
from __future__ import annotations

import argparse
import sys

import cobra

from vcell.audit import audit_model


def load_model(path: str) -> cobra.Model:
    if path.endswith(".json"):
        return cobra.io.load_json_model(path)
    if path.endswith((".xml", ".sbml")):
        return cobra.io.read_sbml_model(path)
    if path.endswith(".mat"):
        return cobra.io.load_matlab_model(path)
    raise ValueError(f"unsupported model format: {path}")


def predicted_essential(model: cobra.Model, tol: float) -> list[str]:
    wt = model.slim_optimize(error_value=0.0)
    out = []
    for g in model.genes:
        with model:
            g.knock_out()
            gr = model.slim_optimize(error_value=0.0)
        if gr is None or gr != gr or gr < tol * wt:
            out.append(g.id)
    return out


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="vcell")
    sub = p.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("audit-biomass", help="classify FBA-essential genes as biomass_forced vs network_forced")
    a.add_argument("model")
    a.add_argument("--genes", nargs="*", default=None)
    a.add_argument("--tol", type=float, default=0.01)
    a.add_argument("--out", default=None)
    b = sub.add_parser("bernstein-score", help="score a model on the Bernstein 2023 RB-TnSeq benchmark (their pipeline and PR-AUC)")
    b.add_argument("model", help="model name in the Bernstein repo (e.g. iML1515) or an SBML .xml path")
    b.add_argument("--bernstein-dir", default="external/E_coli_GEM_validation")
    b.add_argument("--supplement", nargs="*", default=[], help="intracellular metabolite ids to supply (e.g. btn thf)")
    b.add_argument("--out", default=None, help="write JSON summary here")
    r = sub.add_parser("rescue-audit", help="flag supplement rescues whose gene is off the supplement's pathway (likely FBA artifacts)")
    r.add_argument("model")
    r.add_argument("--supplement", nargs="+", required=True, help="MET=PW1,PW2 e.g. btn_c=00780")
    r.add_argument("--gene-pathways", required=True, help="KEGG link TSV: 'eco:b0001<TAB>path:eco00290'")
    r.add_argument("--genes", nargs="*", default=None)
    r.add_argument("--tol", type=float, default=0.01)
    r.add_argument("--out", default=None, help="write CSV here")
    args = p.parse_args(argv)
    if args.cmd == "rescue-audit":
        import csv
        from vcell.rescue import rescue_audit, load_kegg_links
        sup = {}
        for item in args.supplement:
            met, _, pws = item.partition("=")
            sup[met] = {x for x in pws.split(",") if x}
        rows = rescue_audit(load_model(args.model), sup, load_kegg_links(args.gene_pathways), args.genes, args.tol)
        if args.out:
            with open(args.out, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=["gene", "supplement", "growth_ko", "growth_rescued", "label", "flag"])
                w.writeheader(); w.writerows(rows)
        n_off = sum(r_["label"] == "off_pathway" for r_ in rows)
        print(f"{len(rows)} rescues: {len(rows) - n_off} on_pathway, {n_off} off_pathway (likely artifacts)")
        return 0
    if args.cmd == "bernstein-score":
        import json
        from vcell.bernstein import run_benchmark
        r = run_benchmark(args.model, args.bernstein_dir, args.supplement)
        summ = {k: v for k, v in r.items() if k not in ("sim", "fit", "genes", "carbon")}
        print(json.dumps(summ, indent=1))
        if args.out:
            json.dump(summ, open(args.out, "w"), indent=1)
        return 0
    model = load_model(args.model)
    genes = args.genes if args.genes else predicted_essential(model, args.tol)
    df = audit_model(model, genes)
    counts = df["verdict"].value_counts().to_dict()
    if args.out:
        df.to_csv(args.out, index=False)
    print(f"audited {len(df)} genes: {counts}")
    bf = df[df.verdict == "biomass_forced"]
    for _, r in bf.iterrows():
        print(f"  biomass_forced {r['gene']}: needs {','.join(r['unproducible'])}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
