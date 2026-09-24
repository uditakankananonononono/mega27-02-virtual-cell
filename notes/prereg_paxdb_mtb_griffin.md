# Pre-registration: Mtb abundance test with an independent essentiality label (Griffin et al. 2011)

Written 2026-09-25 04:19 IST, before computing. Label replication of H.21.
- Labels: Griffin et al. 2011 PLoS Pathog Table 2 (Tn-seq, H37Rv, glycerol medium; mirror github.com/ajinich/mtb_tn_db
  data/SI_datasets/2011_Griffin_Sassetti/table_2.xlsx). Essential = p_val < 0.05 (the paper's criterion); non-essential = p_val >= 0.05.
- Abundance: the same 14 unique PaxDb Mtb datasets as H.21 (duplicate rule unchanged).
- Primary: per-dataset AUROC; one-sided sign test over the 14; p < 0.05 -> REPLICATES, else NOT SHOWN.
- Secondary: agreement of Griffin and DeJesus labels (Cohen's kappa) on genes called in both.
