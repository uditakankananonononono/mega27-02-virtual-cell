"""Biomass-constituent classes for iJO1366 (used to stratify the provenance audit).
cofactor_vitamin: trace cofactors/prosthetic groups/vitamins - candidates for carry-over or
cross-feeding artifacts (Bernstein et al. 2023 hypothesis class).
bulk: amino acids, nucleotides, lipids, cell wall, LPS - required in bulk; loss is plausibly lethal."""
COFACTOR_VITAMIN = {
    '10fthf_c', '2fe2s_c', '4fe4s_c', '2ohph_c', 'amet_c', 'bmocogdp_c', 'btn_c', 'coa_c', 'fad_c',
    'mlthf_c', 'mobd_c', 'nad_c', 'nadp_c', 'pheme_c', 'pydx5p_c', 'ribflv_c', 'sheme_c', 'thf_c',
    'thmpp_c', 'udcpdp_c'}
INORGANIC = {'ca2_c', 'cl_c', 'cobalt2_c', 'cu2_c', 'fe2_c', 'fe3_c', 'k_c', 'mg2_c', 'mn2_c', 'nh4_c',
             'ni2_c', 'so4_c', 'zn2_c', 'h2o_c'}


def compound_class(met_id: str) -> str:
    if met_id in COFACTOR_VITAMIN:
        return 'cofactor_vitamin'
    if met_id in INORGANIC:
        return 'inorganic'
    return 'bulk'


def gene_class(unproducible: list[str]) -> str:
    """cofactor_only if every unproducible compound is a cofactor/vitamin; else bulk_involved."""
    if not unproducible:
        return 'none'
    return 'cofactor_only' if all(compound_class(m) == 'cofactor_vitamin' for m in unproducible) else 'bulk_involved'
