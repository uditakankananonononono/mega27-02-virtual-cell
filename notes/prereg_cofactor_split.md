# Pre-registration: data-selected cofactor supplements on the Bernstein benchmark (written 9:33 PM IST, before any attribution result)

Base: Bernstein all-corrections iML1515 (published SOTA, PR-AUC 0.764 in results/bernstein/iML1515_bernstein_allcorr.json).
Candidates: the 19 audit-derived cofactors (results/auto_cofactor_set.json).
Split: carbon sources sorted by name; even positions = selection half, odd positions = held-out half.
Selection: greedy forward addition of single cofactors maximising PR-AUC (Bernstein metric) on the selection half, stop when no gain.
Test: the selected set is re-run with the FULL pipeline (run_bernstein_benchmark.py) and scored on the held-out half and on all 25.
Break criterion: held-out-half PR-AUC above allcorr's held-out-half PR-AUC, with paired bootstrap over genes (2,000 resamples) 95% CI of the difference excluding 0.
If the selected set is empty or the CI includes 0: NEGATIVE, preserved.

## Addendum (9:56 PM IST, before any result on the new base)
Comparator correction: Bernstein et al. report their fully corrected iML1515 analysis (all-corrections model + 5 vitamins btn, pnto__R, thm, thf, nad) at PR-AUC 0.843 (0.844 with TALA=0). The 0.764 base above is the all-corrections model WITHOUT vitamins and is not the published best.
New base: all-corrections + 5 vitamins (tag iML1515_allcorr_vit5), reproduced with their pipeline.
Same procedure: single-supplement rescue attribution of the remaining 19-set compounds on that base, greedy selection on the same selection carbons, full-run verification, paired gene bootstrap on held-out carbons.
Break requires held-out PR-AUC above the reproduced allcorr_vit5 held-out PR-AUC with 95% CI excluding 0, AND all-25 PR-AUC above 0.843.
The first-round result (btn, thmpp, coa over allcorr without vitamins) is reported as a secondary finding only.
