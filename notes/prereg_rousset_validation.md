# Pre-registration: external validation on genome-wide CRISPRi (Rousset et al. 2018 PLoS Genet)
Written before any model score is compared with these data (git commit time authoritative).
Data: Rousset et al. 2018, doi:10.1371/journal.pgen.1007749, S12 table (gene-level median log2 fold-change of
guides on the coding strand, `median_coding`), via Europe PMC supplementaryFiles for PMC6242692.
Genes mapped to b-numbers by UniProt primary gene name (data/external/uniprot_ecoli.tsv); unmapped genes dropped and counted.
Label: CRISPRi-essential = median_coding <= -5 (strong depletion); fixed before looking.
Scores (all out-of-fold, trained only on Gerdes 2003 labels): v2 LR (pre-specified primary model), FBA minimal (fba_min).
PRIMARY: AUROC of v2 LR vs the CRISPRi label; paired gene bootstrap (2000) of AUROC(v2 LR) - AUROC(fba_min).
v2 LR generalises if AUROC > 0.70 and the paired CI excludes 0.
Secondary: Spearman between v2 LR score and -median_coding; the same restricted to genes not in the Gerdes positives.
Caveat stated in advance: CRISPRi is polar within operons, so downstream-gene effects inflate apparent essentiality.
