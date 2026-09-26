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
