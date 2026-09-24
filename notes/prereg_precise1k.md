# Pre-registration: transcriptome features from PRECISE-1K / iModulonDB (v9)
Written before the expression files are downloaded or scored (git commit time authoritative).
Hypothesis: essential metabolic genes are expressed at high, stable levels across conditions and sit in fewer
condition-specific regulatory modules, so expression adds signal beyond v2.
Data: SBRG PRECISE-1K (Lamoureux et al. 2023 Nat Commun; github.com/SBRG/precise1k, data/precise1k):
log_tpm_qc.csv (QC-passed samples), M.csv (gene weights per iModulon), imodulon_table.csv (per-iModulon threshold).
Features: mean log-TPM across samples; SD of log-TPM; number of iModulons with |M| > that iModulon's threshold.
Genes not in the compendium: training-fold median imputation plus a missing flag.
Model/folds identical to v2. v9_lr = v2 features + 3 features (+flag). PRIMARY: paired 2000-gene bootstrap
AUROC(v9_lr) - AUROC(v2_lr) on Gerdes labels; success iff CI lower bound > 0. Secondary: AUPRC; Rousset CRISPRi difference;
univariate AUROCs. Nulls reported as nulls.
