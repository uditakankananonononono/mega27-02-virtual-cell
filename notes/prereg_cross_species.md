# Pre-registration: cross-species replication of the supplement-bypass finding

Written 2026-09-24 23:50 IST, before any fitness value of the second organism is
compared against rescue labels. Committed before the analysis script is run.

## Hypothesis (from H.7, E. coli, results/pathway_concordance.json)
In-silico supplement rescues where the knocked-out gene lies on the supplement's own
KEGG pathway have near-neutral measured fitness; off-pathway rescues are deleterious.

## Organism and data (independent of the E. coli discovery set)
- Shewanella oneidensis MR-1. GEM iMR1_799 (Pinchuk et al. 2010), PSAMM collection
  (github.com/zhanglab/psamm-model-collection, sbml/iMR1_799), converted with psamm sbmlexport.
- RB-TnSeq fitness: Price et al. 2018 Nature, genomics.lbl.gov/supplemental/bigfit/html/MR1/
  (fit_logratios_good.tab, fit_genes.tab, expsUsed).
- Gene -> pathway: KEGG REST link/pathway/son.
- Other organisms may be added only with this same protocol and must be reported
  whatever the result.

## Protocol
1. Medium: iMR1_799 on a defined minimal medium (lactate carbon source, O2, NH4, Pi, SO4,
   minerals), no vitamins or cofactors. If the model cannot grow on it, the minimal set of
   extra exchanges needed is reported and added.
2. Supplement panel, same as E. coli: biotin (00780), thiamine (00730), pantothenate
   (00770), THF (00790, 00670), NAD (00760), SAM (00270), PLP (00750), supplied as
   cytosolic sinks with `vcell rescue-audit` (tol 0.01, supply 10).
3. Gene label: on_pathway if any rescuing supplement shares a KEGG pathway with the gene,
   else off_pathway.
4. Outcome: per-gene median log2 fitness over all experiments in fit_logratios_good.tab.
   Genes absent from the fitness table (no usable insertions; often essential) are
   counted per group as a secondary outcome, not imputed.

## Criteria
- Primary: one-sided Mann-Whitney (on > off) on gene-level median fitness, p < 0.05 AND
  median(on) - median(off) > 0 -> REPLICATED.
- median difference <= 0 -> FALSIFIED in MR-1.
- fewer than 5 genes with fitness in either group -> UNDERPOWERED (reported as such, no claim).
- Secondary: fraction of genes absent from the fitness table, off vs on (prediction: off higher).

## Amendment 1 (2026-09-24 23:54 IST, after MR-1 was UNDERPOWERED; before any P. putida analysis)
MR-1 gave zero off-pathway rescues, so the test could not run. Add Pseudomonas putida
KT2440 under the same protocol and criteria, with these organism-specific inputs:
- GEM iJN1463 (BiGG, bigg.ucsd.edu/static/models/iJN1463.json.gz); medium: the model's
  own exchanges closed except glucose 10 (EX_glc__D_e), O2, NH4, Pi, SO4, water, protons,
  CO2 and inorganic ions; no vitamins.
- Fitness: Fitness Browser February 2024 release (figshare 10.6084/m9.figshare.25236931),
  db.StrainFitness.Putida. Gene fitness per experiment = mean strain fitness over strains
  with used == TRUE; gene outcome = median over experiments. (The gene-level table is only
  in the 2.3 GB feba.db, so strain data are aggregated.)
- Gene -> pathway: KEGG link/pathway/ppu.
- Supplements: same seven cofactors, BiGG ids btn_c, thm_c, pnto__R_c, thf_c, nad_c, amet_c, pydx5p_c.
Verdict rules unchanged; each organism is reported separately.
