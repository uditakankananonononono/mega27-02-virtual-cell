# Pre-registration: OMA 1:1 ortholog breadth as an essentiality feature (v7)
Written before any ortholog count is fetched or scored (git commit time authoritative).
Hypothesis: essential genes are conserved more widely, so the number of OMA 1:1 orthologs (and number of distinct
species/kingdom-level breadth) adds signal beyond frozen v2 LR.
Data: OMA Browser REST /api/protein/<UniProt>/orthologs/?rel_type=1:1 for model genes (b-number -> UniProt via
data/external/uniprot_ecoli.tsv). Features: log1p(n 1:1 orthologs); fraction of orthologs outside Enterobacterales is NOT
computed (taxonomy not returned uniformly); a gene absent from OMA gets 0.
Model/folds: identical to v2 (LR pipeline, StratifiedKFold(3, shuffle, rs=7)). v7_lr = v2 features + log1p(n_orth_1to1).
PRIMARY: paired 2000-gene bootstrap AUROC(v7_lr) - AUROC(v2_lr) on Gerdes labels; success iff CI lower bound > 0.
Secondary: univariate AUROC of the feature; same difference on Rousset 2018 CRISPRi labels. Nulls reported as such.
