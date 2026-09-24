"""Uptake plausibility for supplement rescues.

A model rescue by an intracellular supplement is only biologically meaningful if the organism can import the
compound. `uptake_systems` returns the organism's TCDB transport systems whose curated substrates overlap a set of
ChEBI IDs (primary and secondary IDs should both be supplied, because TCDB uses secondary IDs).
"""
from __future__ import annotations
import re


def load_tcdb_substrates(path: str) -> dict[str, set[str]]:
    """TCDB getSubstrates table: 'TC<TAB>CHEBI:1;name|CHEBI:2;name'."""
    out = {}
    for line in open(path):
        tc, _, subs = line.rstrip("\n").partition("\t")
        if tc:
            out[tc] = set(re.findall(r"CHEBI:\d+", subs))
    return out


def load_organism_tc(path: str) -> dict[str, list[str]]:
    """UniProt TSV with columns Entry, Gene Names (primary), TCDB ('2.A.21.1.1;')."""
    out: dict[str, list[str]] = {}
    lines = open(path).read().splitlines()
    hdr = lines[0].split("\t")
    i_e, i_g, i_t = hdr.index("Entry"), hdr.index("Gene Names (primary)"), hdr.index("TCDB")
    for ln in lines[1:]:
        c = ln.split("\t") + [""] * 3
        for tc in [t for t in c[i_t].split(";") if t]:
            out.setdefault(tc, []).append(c[i_g] or c[i_e])
    return out


def uptake_systems(chebi_ids: set[str], tcdb: dict[str, set[str]], org_tc: dict[str, list[str]]) -> dict[str, list[str]]:
    return {tc: genes for tc, genes in org_tc.items() if tcdb.get(tc, set()) & set(chebi_ids)}
