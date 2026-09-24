# Pre-registration: protein abundance vs essentiality outside E. coli (PaxDb x Fitness Browser)

Written 2026-09-25 04:08 IST, before any abundance-vs-fitness number was computed for these organisms.

## Question
H.19 found abundance predicts essentiality in all 19 E. coli PaxDb datasets. Does this hold in other bacteria?

## Data (all public, accession level)
PaxDb v5 per-dataset abundance files (pax-db.org/downloads/latest/datasets/<taxid>/):
- 160488 P. putida KT2440: PXD003826_Lidbury
- 1140 S. elongatus PCC 7942: PXD000510_Guerreiro
- 211586 S. oneidensis MR-1: GPM_201408
- 226186 B. thetaiotaomicron: GPM_201408, PXD014877_Mueller
- 882 D. vulgaris Hildenborough: 4 x zhang_2006 (Form/Lac x Exp/Stat)
Integrated (WHOLE_ORGANISM) files are excluded - they are derived.
Fitness Browser (Price et al. 2018) gene median fitness tables and FB_aaseqs protein sequences (already on disk).
STRING v12 protein sequences for the 5 taxa, used only to map PaxDb STRING IDs to FB locusIds by exact sequence match.

## Labels
Putative essential = protein in FB_aaseqs for that organism with no row in the gene median fitness table
(FB reports no fitness for genes with too few insertions, the standard RB-TnSeq essentiality proxy).
Nonessential = has a fitness row.

## Primary test
Per dataset: AUROC of log10 abundance for putative essential vs nonessential, over mapped genes detected in the dataset.
Prediction: AUROC > 0.5. Verdict across 9 datasets by one-sided sign test:
p < 0.05 -> GENERALISES; otherwise NOT SHOWN. Also report per-organism inverse-variance pooled AUROC (Hanley-McNeil SE).

## Secondary
1. Spearman(abundance, median fitness) among nonessential genes; prediction negative.
2. Length control: logistic regression essential ~ log abundance + log protein length; report abundance coefficient sign per dataset.
Caveat stated in advance: short genes also lack fitness rows for technical reasons, so the primary label is noisy.
