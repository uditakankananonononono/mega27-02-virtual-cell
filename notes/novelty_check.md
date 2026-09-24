# Novelty check: essentiality-provenance audit (2026-09-24)

Prior art found (read, not recalled):
- Bernstein et al. 2023, Mol Syst Biol, "Evaluating E. coli genome-scale metabolic model accuracy with
  high-throughput mutant fitness data" (https://pmc.ncbi.nlm.nih.gov/articles/PMC10698504/ ;
  doi 10.15252/msb.202311566). iML1515 error analysis: 21 vitamin/cofactor biosynthesis genes
  (biotin, R-pantothenate, thiamin, THF, NAD+) gave false-essential calls; fixed by MANUALLY adding
  those cofactors via intracellular exchange; carry-over hypothesis. MoCo/molybdopterin not in their list.
- Joyce/Palsson-lab style error taxonomy: "Three factors underlying incorrect in silico predictions of
  essential metabolic genes" (https://pmc.ncbi.nlm.nih.gov/articles/PMC2248557/): names incomplete
  biomass definition as a cause of errors; no per-gene automated attribution.
- Xavier et al. 2017, "Integration of Biomass Formulations ... Universally Essential Cofactors in
  Prokaryotes" (https://pmc.ncbi.nlm.nih.gov/articles/PMC5249239/): curates which cofactors belong in biomass.
- MEMOTE (https://memote.readthedocs.io/en/latest/autoapi/memote/support/biomass/index.html):
  find_blocked_biomass_precursors / precursor-production tests run on the wild-type model; they do not
  attribute each knockout's lethality to specific biomass compounds.

Honest positioning:
- The phenomenon (biomass-cofactor-driven false essentiality) is KNOWN. Not claimed as a discovery.
- What is ours: an automated, per-gene, falsifiable rescue test (knock out -> find unproducible biomass
  compounds -> strip exactly those -> re-solve) that labels every FBA-essential call biomass_forced vs
  network_forced and names the compounds; packaged as `python -m vcell audit-biomass MODEL`.
  We are not aware of an existing tool that does this per-gene attribution; this is a hedge, not a proof.
- The MoCo cluster (bmocogdp_c, mobd_c) is a compound class outside Bernstein's five - an extension, to be
  quantified by the full audit.

## Supplement-bypass artifact (pathway concordance), checked 11:44 PM IST
Searches (web_search, 4 queries): GEM supplementation false rescue off-pathway; SAM supplementation FBA spurious rescue;
metabolite supplementation essentiality cross-feeding RB-TnSeq; Bernstein follow-ups.
Closest work: Bernstein et al. 2023 (vitamin/cofactor cross-feeding, on-pathway only; https://pmc.ncbi.nlm.nih.gov/articles/PMC10698504/);
SAM catabolism biology ("Excess S-adenosylmethionine inhibits methylation via catabolism to adenine",
https://www.nature.com/articles/s42003-022-03280-5) shows SAM is catabolised in vivo, consistent with the model using it as a nutrient;
"Model-driven analysis of mutant fitness experiments" (https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1008137).
No paper found that tests on- vs off-pathway rescue concordance against fitness. Not exhaustive: novelty is PROBABLE, not proven.
