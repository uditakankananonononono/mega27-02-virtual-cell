# Pre-registration: pathway concordance of supplement rescue (written 11:43 PM IST, before the pathway split was computed)

Observation that motivated it (already seen): pairs rescued by the 5 Bernstein vitamins have median RB-TnSeq fitness -0.69,
pairs rescued by SAM/PLP/folates on the corrected all-corrections model have median -3.5 to -4.1.
Hypothesis (H-PC, "supplement-bypass artifact"): a model rescue is biologically real (high fitness; cross-feeding or carry-over)
when the knocked-out gene lies in the supplement's own biosynthesis pathway, and spurious (low fitness) when the supplement
rescues an off-pathway gene by being metabolised as a nutrient/one-carbon source.
Pathway map (KEGG eco): btn 00780; thm 00730; pnto__R 00770; thf 00790+00670; nad 00760; amet 00270; pydx5p 00750; 10fthf/mlthf 00670+00790.
Rescued pairs: iML1515 base -> +5 vitamins (their Part 6), and single-supplement rescues on the corrected all-corrections base
(results/bernstein/attrib_allcorr_fixed.npz). A vitamin-rescued pair is on-pathway if the gene is in any of the 5 vitamin pathways.
Primary test: one-sided Mann-Whitney U, on-pathway fitness > off-pathway fitness, over pairs; effect = difference of medians with
gene-cluster bootstrap 95% CI (2,000 resamples of genes). Falsified if p >= 0.05 or the CI includes 0.
Secondary: same test within the vitamin set only (removes the supplement-identity confound); fraction of pairs with fitness > -2.
