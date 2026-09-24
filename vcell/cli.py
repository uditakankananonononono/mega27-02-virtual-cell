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
    args = p.parse_args(argv)
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
