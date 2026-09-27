# Provisional accession reconciliation, 27 September 2026

The strict gate audit (`scripts/audit_strict_gates.py`, `results/datasets_accession_audit.csv`) left 88 manifest rows `provisional` because provider accession and usage were not individually verified. This reconciliation does NOT loosen that standard: a row is resolved only when (a) the local file exists and hashes to a recorded value, (b) a real provider identifier is present (PaxDb `#id:` header, BiGG survey accession, KEGG organism code), and (c) scored-use evidence exists in a saved result. The original audit CSV is untouched; output is an overlay (`results/provisional_accession_reconciliation.csv`).

Resolved 65 of 88: 57 PaxDb datasets across six scored arms (E. coli 19, M. tuberculosis 14, B. subtilis 5, P. aeruginosa 10, cross-species 9), 5 BiGG models with survey hashes, 2 KEGG maps (son, ppu) whose saved bytes exactly match live KEGG REST on 2026-09-27 and which have scored usage, and the TCDB substrate table (live row-set match; byte order differs). Seven other KEGG maps stay provisional: six have no scored usage evidence found, and `sme` drifted from live (5,843 saved vs 5,841 live rows). 23 rows remain provisional (label papers, EcoliWiki, SubtiWiki, Rousset table, GO/EcoCyc, COG, OMA, AlphaFold, PRECISE-1K, iMR1_799 XML, unused/drifted KEGG maps). Technical replicates within one study remain one study; this does not change the 123 model-inclusive gate floor or claim independent wet-lab datasets.

## Extension pass (2026-09-27, second overlay)

15 further rows resolved against local hashed files with in-repo scored/scripted use:
EcoCyc GAF (GO annotations), iMR1_799 GEM (cross-species MR1 arm), Rousset 2018 S12 (CRISPRi benchmark),
AlphaFold per-gene metrics, OMA HOG orthology, NCBI COG-2020, PRECISE-1K compendium (pinned commit),
DeJesus 2017 Table S3, SubtiWiki essential-gene page, Poulsen 2019 Datasets S5/S6 (incl. Turner/Lee/Skurnik/Liberati
label replications), Griffin 2011 Table 2.

80/88 provisional rows now reconciled. 8 remain provisional with stated reasons:
- EcoliWiki essential-gene list and GO go-basic.obo: no local hashed file in-repo.
- KEGG link/pathway codes ccs, evi, dsu, pact, sme, psb: enumerated in the audit but have no scored use in any
  result; resolving them would require fetching sources the analyses never consumed.

Resolution = provider identity + hashed local copy + cited in-repo use. It does not upgrade technical
replicates to independent studies and does not change the owner-approved 123 model-inclusive gate floor.
