# Pre-registration: COG functional category / phyletic spread and PDB structure coverage (v8)
Written before the COG file is downloaded or any feature is scored (git commit time authoritative).
Data: NCBI COG 2020 (ftp.ncbi.nlm.nih.gov/pub/COG/COG2020/data: cog-20.cog.csv, cog-20.def.tab), rows for assembly
GCF_000005845.2 (E. coli K-12 MG1655), matched to b-numbers by locus tag / RefSeq protein; UniProt xref_pdb for UP000000625.
Features: (a) one-hot COG functional category letters (J, K, L, C, E, F, G, H, I, P, Q, M, D, O, T, V, S, R, other);
(b) log1p(number of COG 2020 genomes in the gene's COG) [phyletic spread; max over COGs if several];
(c) log1p(number of PDB entries). Genes without COG get all-zero one-hot plus an explicit no_cog flag.
Model/folds identical to v2 (LR, StratifiedKFold(3, shuffle, rs=7)). v8_lr = v2 features + (a)+(b)+(c).
PRIMARY: paired 2000-gene bootstrap AUROC(v8_lr) - AUROC(v2_lr) on Gerdes labels; success iff CI lower bound > 0.
Secondary: v8 without (c) (PDB counts track how well studied a protein is, a known leakage route); Rousset CRISPRi difference.
Nulls reported as nulls.
