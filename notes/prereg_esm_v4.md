# Pre-registration: ESM-2 protein-language-model feature (v4), written 9:53 PM IST before embeddings exist

Feature: ESM-2 t6-8M (fair-esm) mean-pooled embedding (320-d) per U00096.3 CDS (scripts/esm_embed.py).
Primary model v4_lr: v2_lr feature set + PCA(16) of the ESM embedding, PCA and scaling fit inside each training fold
(ColumnTransformer), same LogisticRegression (C=0.5, balanced), same folds StratifiedKFold(3, shuffle, rs=7).
Comparator: v2_lr (pre-specified primary so far; AUROC 0.797 in results/ensemble_v2.json).
Test: paired bootstrap over genes (2,000, seed 11) of AUROC(v4_lr) - AUROC(v2_lr).
Claim a gain only if the 95% CI excludes 0. Secondary (reported, not claimed): esm_only_lr, v4 with PCA(32).
Negative result is preserved in the paper either way.
