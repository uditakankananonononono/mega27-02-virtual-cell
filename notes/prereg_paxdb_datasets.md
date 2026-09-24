# Pre-registration: is the protein-abundance signal robust across individual proteomics datasets? (PaxDb v5)
Written before the individual PaxDb datasets are downloaded or scored (git commit time authoritative).
Context: the integrated PaxDb abundance (pax_log_ppm) is one of the v3 features. Integrated values can hide
dataset-specific artefacts, so we score each individual E. coli K-12 (taxon 511145) dataset in PaxDb separately.
Data: every file in pax-db.org/downloads/latest/datasets/511145/ except WHOLE_ORGANISM-integrated (already used).
Mapping: PaxDb string ids 511145.bXXXX -> b-number. Genes: the 1249 model genes with Gerdes labels (results/aligned_predictions_all_models.csv).
Per dataset: coverage; AUROC of log10(abundance) for essentiality among covered genes (missing genes excluded).
PRIMARY: exact two-sided sign test over datasets of AUROC > 0.5 vs < 0.5; the signal is robust if p < 0.05 and the
majority exceeds 0.5. Secondary: inverse-variance (Hanley-McNeil SE) random-effects-free pooled AUROC; Spearman of
AUROC vs coverage (do low-coverage datasets differ?). All datasets reported.
