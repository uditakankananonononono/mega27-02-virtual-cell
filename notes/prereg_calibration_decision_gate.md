# PREREG VC2-R1 (draft, lock before scoring): ensemble decision-quality gate
Date: 2026-09-26. Author: revival builder. Status: LOCKED at commit time; any amendment gets a new dated section.
Question: can the stacked essentiality ensemble beat the iJO1366 FBA minimal
baseline on BOTH ranking and hard calls on the frozen Gerdes 2003 overlap?
Locked gates:
 G1 (ranking): OOF AUROC of final predictor > 0.666 (iJO1366 FBA minimal AUROC
   on same 125-gene overlap, repo-measured). [already held by 0.7225 ensemble]
 G2 (hard calls): McNemar discordant pairs vs FBA minimal must NOT be a
   significant loss (two-sided p >= 0.05), AND net discordant wins >= 0.
   Current state: 27:50 loss, p=0.012 -> FAIL. This is the gap to close.
Permitted levers (frozen): Platt/isotonic calibration of OOF scores, threshold
 selection on TRAIN folds only (never the eval fold), class-weight re-fit of
 the v2 LR / oversampled-LR head, probability averaging with FBA-flux sign
 features. Forbidden: threshold tuning on the evaluation fold, re-scoring the
 Gerdes labels, dropping genes from the overlap.
Evaluation: same frozen 3-fold split as results/ v2 (hashes recorded in repo);
 report AUROC, F1, discordant counts + exact McNemar p, all seed-fixed.
Negative handling: if G2 cannot pass after the lever set is exhausted, record
 the negative and pivot per rule 6 (ChatGPT redirection), e.g. to a
 selective-classification gate (coverage-curve beat) with a fresh prereg.

Amendment A1 (2026-09-26, locked before any gate evaluation under this protocol)
Baseline reproduction (harness check, not a new result): the locked "current
state" numbers reproduce exactly from repo artifacts - v1 ensemble OOF
(results/ensemble_oof.csv) thresholded at its best-F1 point over a 99-point
quantile grid vs fba_min thresholded at its own best-F1 point gives b=27,
c=50, F1s 0.4561 / 0.4119, matching results/ensemble_results.json. Both of
those thresholds were tuned on ALL genes (eval-inclusive) - the flaw this
gate exists to remove.
Locked protocol:
 1. Folds: reconstruct the frozen v2 split via StratifiedKFold(3, shuffle,
    rs=7) over the aligned row order (identical construction to
    scripts/run_ensemble_v2.py; split depends only on label vector + order).
 2. Comparator: FBA minimal hard call with threshold selected per eval fold
    on TRAIN folds only, same MCC rule as the candidates (symmetric
    protocol). The repo's fixed structural rule fba_min > 0.5 is reported as
    a sensitivity, not the gate comparator.
 3. Threshold selection rule (all candidates): per eval fold k, threshold =
    argmax MCC over the other two folds' OOF scores on a 199-point quantile
    grid of the train scores; ties break toward the threshold nearest 0.5.
    Eval-fold labels/scores never enter selection.
 4. Candidate ladder (pre-ordered; stop at first PASS; all rungs reported):
    R1 v2_lr OOF + train-only MCC threshold (primary).
    R2 isotonic calibration of v2_lr (fit on train folds only) + same rule.
    R3 Platt (logistic) calibration of v2_lr (fit on train folds only).
    R4 v2-feature LR re-fit with random oversampling of the minority class
       inside each training fold (seed 7; class_weight off once oversampled).
    R5 probability average (v2_lr + v2_gbm)/2 + train-only MCC threshold.
 5. Gate statistics (pooled over all 1249 genes): discordant counts b, c vs
    the comparator, exact two-sided binomial McNemar p; G2 PASSes iff
    p >= 0.05 AND b - c >= 0. AUROC of the rung's score (unchanged by
    monotone calibration) must remain > 0.666 (G1). F1 and MCC reported.
 6. Multiplicity: the ladder order is the testing order; the first passing
    rung is the reported predictor. If no rung passes, the negative is
    recorded and the pivot clause of the base prereg fires (rule 6).
