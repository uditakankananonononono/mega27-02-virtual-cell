# Pre-registration (exploratory, descriptive): uptake plausibility of rescue-audit supplements
Written before TCDB/ChEBI matching is run (git commit time authoritative). n = 7 supplements, so no significance test.
Question: for each supplement that rescues knockouts in vcell rescue-audit (results/rescue_audit_iML1515_M9.csv), does
E. coli K-12 have a transporter in the Transporter Classification Database whose curated substrates include that compound?
Method: supplement ChEBI IDs from iJO1366 annotations, expanded through the ChEBI ontology (EBI OLS4 API: conjugate
acid/base, tautomer and enantiomer links, one hop); TCDB substrate table (TC -> ChEBI); E. coli K-12 TC assignments from
UniProt UP000000625 xref_tcdb. A supplement is "uptake-supported" if any K-12 TC system lists any expanded ChEBI ID.
Expectation written in advance (from the literature, e.g. ThiBPQ, PanF, biotin uptake): thm, pnto__R, btn supported;
amet (SAM), nad, thf, pydx5p not supported. If SAM, the main off-pathway rescuer, has no K-12 uptake system, that is
independent evidence the SAM rescues are artifacts of supplying an intracellular sink, and rescue-audit gains a
no_known_uptake flag. Any mismatch with the expectation is reported.

## Amendment 1 (after a first run showed an ID-matching defect; disclosed, first-run output kept in git history as none was committed)
TCDB substrate records use secondary ChEBI IDs (e.g. CHEBI:9530 for thiamine(1+), whose primary ID is CHEBI:18385), so the
first run missed ThiBPQ (TC 3.A.1.19.1). The expanded set now also includes every alt_id (secondary ID) of each expanded
term from OLS4. First-run result for the record: only pnto__R_c supported; thm_c and btn_c not.
