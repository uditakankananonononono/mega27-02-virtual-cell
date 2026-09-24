# Pre-registration: protein abundance vs essentiality in M. tuberculosis H37Rv

Written 2026-09-25 04:11 IST, before any abundance-vs-essentiality number was computed for Mtb.

## Data
- PaxDb v5, taxid 83332, all 16 non-integrated per-dataset files (pax-db.org/downloads/latest/datasets/83332/).
- Essentiality: DeJesus et al. 2017 mBio Table S3 (HMM calls, saturating Tn-seq, H37Rv), mirrored in github.com/ajinich/mtb_tn_db data/annotations/DeJesus_mbio.xlsx.
  Essential = ES or ESD; non-essential = NE; GD, GA and Uncertain are excluded.
- Mapping: PaxDb IDs are STRING IDs 83332.<Rv number>, which match the ORF IDs directly.

## Primary
Per dataset AUROC of log10 abundance for essential vs non-essential. Prediction > 0.5.
One-sided sign test over datasets: p < 0.05 -> HOLDS IN MTB, otherwise NOT SHOWN.
Pooled inverse-variance AUROC with 95% CI reported.

## Duplicate rule (dataset counting)
Some files may be re-releases of the same data (e.g. PA_2013-7 and PA_201307). Two datasets with
Spearman > 0.99 over shared proteins and overlap > 90% count as one for the sign test and for the dataset gate.

## Secondary
GD genes should sit between ES and NE in median abundance (reported, not tested).
