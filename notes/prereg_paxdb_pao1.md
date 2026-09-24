# Pre-registration: protein abundance vs essentiality in P. aeruginosa PAO1

Written 2026-09-25 04:18 IST, before any abundance-vs-essentiality number was computed for P. aeruginosa.

- Abundance: PaxDb v5 taxid 208964, all 10 non-integrated datasets. Studies (strict rule): GPM_201408, MSV000096603
  (multi-species, already counted via B. subtilis), PXD004560, PXD009705 (3 files, one study), PXD021907, PXD025827,
  PeptideAtlas PA_2021-1, biomart_22024.
- Essentiality: Poulsen et al. 2019 PNAS (doi:10.1073/pnas.1900570116) Dataset S5 "Essential Genes" (528 genes, PA14 IDs).
  Essential = Essential Category "Core" (essential in all nine strains tested). Non-essential = named proteins not on the
  528-gene sheet at all. Genes on the sheet but not Core are excluded. PAO1 was not one of the nine strains, so
  mapping is by gene name (PaxDb column 1 vs sheet "Gene Name", case-insensitive); unnamed genes are dropped.
- Primary: per-dataset AUROC of log10 abundance; datasets from one study are averaged first, so the sign test runs over studies.
  One-sided sign test p < 0.05 -> HOLDS IN P. AERUGINOSA, else NOT SHOWN.
- Caveat stated in advance: name matching and the core-only label bias toward well-studied, named genes.
