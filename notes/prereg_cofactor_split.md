# Pre-registration: data-selected cofactor supplements on the Bernstein benchmark (written 9:33 PM IST, before any attribution result)

Base: Bernstein all-corrections iML1515 (published SOTA, PR-AUC 0.764 in results/bernstein/iML1515_bernstein_allcorr.json).
Candidates: the 19 audit-derived cofactors (results/auto_cofactor_set.json).
Split: carbon sources sorted by name; even positions = selection half, odd positions = held-out half.
Selection: greedy forward addition of single cofactors maximising PR-AUC (Bernstein metric) on the selection half, stop when no gain.
Test: the selected set is re-run with the FULL pipeline (run_bernstein_benchmark.py) and scored on the held-out half and on all 25.
Break criterion: held-out-half PR-AUC above allcorr's held-out-half PR-AUC, with paired bootstrap over genes (2,000 resamples) 95% CI of the difference excluding 0.
If the selected set is empty or the CI includes 0: NEGATIVE, preserved.
