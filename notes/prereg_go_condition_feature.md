# Pre-registration: fold-internal GO error-term feature (v5), written 11:42 PM IST before any v5 result

Motivation: results/go_error_enrichment.json (computed on all genes) found confident false positives enriched for
molybdopterin biosynthesis and ATP synthase. Using those terms directly as a feature would be circular.
Procedure (nested, no test-fold information): for each outer fold of StratifiedKFold(3, shuffle, rs=7):
 1. inner 3-fold OOF v2_lr scores on the training genes only;
 2. training confident FPs = non-essential training genes in the top 10% of inner OOF scores;
 3. one-sided Fisher test per propagated GO term (terms with >=3 training genes), BH q<0.05, keep up to 10 terms;
 4. feature go_err = 1 if a gene carries any kept term (propagated), else 0;
 5. fit v2 features + go_err (same LR, C=0.5, balanced) on training genes; score the test fold.
Primary: v5_lr AUROC vs v2_lr AUROC, paired bootstrap over genes (2,000, seed 11). Gain claimed only if 95% CI excludes 0.
Secondary: AUPRC; list of terms kept per fold. Negative preserved either way.
