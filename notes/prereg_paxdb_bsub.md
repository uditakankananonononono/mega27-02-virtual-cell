# Pre-registration: protein abundance vs essentiality in B. subtilis 168

Written 2026-09-25 04:14 IST, before any abundance-vs-essentiality number was computed for B. subtilis.

- Abundance: PaxDb v5 taxid 224308, the 5 non-integrated datasets (PRIDE, GPM_201408, MSV000096603, PXD006444_control, PXD014877_Mueller).
- Essential: genes linked from the SubtiWiki "Essential genes" page (retrieved 2026-09-25; list based on Koo et al. 2017 deletion libraries).
  Matched by gene name (PaxDb column 1, case-insensitive). All other detected proteins with a gene name = non-essential (proxy).
- Primary: per-dataset AUROC of log10 abundance; one-sided sign test; p < 0.05 -> HOLDS IN B. SUBTILIS, else NOT SHOWN.
  With 5 datasets the minimum achievable p is 0.031, so all 5 must exceed 0.5.
- Counting: Mueller PXD014877 and MSV000096603 are multi-species studies already partly counted; under the strict rule
  Mueller collapses with the E. coli/B. theta entries, MSV000096603 counts once across species.
