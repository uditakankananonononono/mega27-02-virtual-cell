# PREREG VC2-R2 (draft): from descriptive census to verified cross-model finding
Date: 2026-09-26. Status: LOCKED at commit time; any amendment gets a new dated section.
Background: results/bigg_model_survey.csv - 108/108 BiGG v2 models fetched,
 parsed, hashed; 68/108 include a MoCo-named biomass constituent. Descriptive
 only; no growth scored. The E. coli iJO1366 finding (MoCo objective causes
 false moaD/moaC/moaE/mobA/moeB essentiality) is verified in ONE model.
Locked question: in how many of the 68 MoCo-positive models does removing the
 MoCo constituent(s) from the biomass objective flip at least one gene
 essentiality call under default medium? (structural falsifiability census)
Locked protocol:
 1. Load each hashed model JSON with COBRApy (versions pinned); verify SHA-256
    against bigg_model_survey.csv before use. No re-download.
 2. Default bounds as shipped; single-gene deletion via cobra.knock_out(),
    essential = growth < 1e-6 (same threshold as iJO1366 analysis).
 3. Compare (a) baseline objective vs (b) objective with MoCo-named
    constituents zeroed; count genes flipping essential->non-essential with
    retained growth >= 0.95 wild-type.
 4. Report per-model flip counts + the MoCo-gene-family labels; primary
    endpoint = fraction of the 68 models with >=1 flip (exact binomial CI).
 Falsifiability: any fraction in [0,1] is a valid result; the DISCOVERY claim
 requires >=10% of models with >=1 flip (pre-locked threshold) AND at least
 one independently re-verified model (second model family beyond iJO1366
 reproduced by manual inspection).
 Compute: 68 models x ~1-2k genes, COBRA LP - feasible on 2 cores (hours).
Negative handling: a null (<10%) is preserved as a documented negative and
 the claim reverts to the verified single-model finding; rule-6 pivot then
 targets the rescue-audit cross-species power gap instead.

## Amendment A1 (2026-09-26, before any scoring run)
Scorability: a model whose baseline (unmodified objective) FBA growth on its
stored medium bounds is < 1e-6 cannot produce an essentiality call and is
NON-SCORABLE. Primary fraction denominator = scorable MoCo-positive models;
non-scorable models are listed with their baseline growth value. iJO1366
serves as the positive control: the harness must reproduce the known
moaD-cluster rescue (essential at baseline, rescued to >=0.95 WT growth after
MoCo-constituent removal) before any other model result is accepted.
