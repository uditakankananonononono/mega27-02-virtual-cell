"""Loaders for experimental ground truth and genome data.

Ground truth: Gerdes et al. 2003 (J. Bacteriol. 185:5673) systematic
transposon-footprinting essentiality assertions for E. coli MG1655 -
4,294 protein-coding genes, each asserted E (essential) / N (nonessential)
/ X (conflicting) / ? (no data). We use E and N only; X and ? are excluded
from benchmarks (documented in the paper).
"""
from __future__ import annotations

import pandas as pd


def parse_gerdes_s1(path: str) -> pd.DataFrame:
    """Parse Gerdes 2003 Supplementary Table S1 into a tidy frame."""
    rows = []
    with open(path) as fh:
        for line in fh:
            f = [c.strip().strip('"') for c in line.rstrip("\n").split("\t")]
            if len(f) < 11 or not f[0].isdigit():
                continue  # headers / blank lines
            rows.append({
                "start": int(f[0]), "ergo_id": f[1], "gene_name": f[2],
                "assertion": f[3], "length_bp": int(f[4]) if f[4].isdigit() else None,
                "sp_id": f[7], "description": f[8], "bnumber": f[9],
            })
    df = pd.DataFrame(rows)
    df = df[df["bnumber"].str.match(r"^b\d{4}$")]
    return df


def essentiality_labels(gerdes: pd.DataFrame) -> pd.DataFrame:
    """Binary experimental labels: 1 = essential (E), 0 = nonessential (N)."""
    lab = gerdes[gerdes["assertion"].isin(["E", "N"])].copy()
    lab["essential"] = (lab["assertion"] == "E").astype(int)
    return lab[["bnumber", "gene_name", "essential"]]
