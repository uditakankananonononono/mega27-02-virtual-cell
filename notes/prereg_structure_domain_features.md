# Pre-registration: protein structure confidence and Pfam domain features (v6)
Written before any feature is computed or scored (git commit time authoritative).
Hypothesis: essential E. coli proteins are more ordered (higher AlphaFold pLDDT) and carry Pfam families that are
essential in other members, so these features add signal beyond the frozen v2 LR feature set.
Data: AlphaFold DB API (/api/prediction/<UniProt>) fields globalMetricValue (mean pLDDT) and fractionPlddtVeryLow;
UniProt REST xref_pfam for proteome UP000000625. b-number to UniProt via data/external/uniprot_ecoli.tsv.
Features: plddt_mean, plddt_frac_vlow, n_pfam, pfam_te. pfam_te = max over a protein's Pfam families of the smoothed
essential rate (ess+1)/(n+2) of OTHER TRAINING-FOLD genes sharing the family, computed inside each fold only (no leakage);
0.5*prior if no family or no training-fold neighbours. Missing AlphaFold values imputed with training-fold median.
Model: same LR pipeline, same StratifiedKFold(3, shuffle, rs=7) as v2. v6_lr = v2 features + the 4 new features.
PRIMARY: paired 2000-gene bootstrap of AUROC(v6_lr) - AUROC(v2_lr) on Gerdes labels. Success iff 95% CI lower bound > 0.
Secondary: same on AUPRC; the same difference on the Rousset 2018 CRISPRi labels (notes/prereg_rousset_validation.md).
A null or negative result is reported as such.
