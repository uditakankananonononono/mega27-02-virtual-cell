"""Main-body 2026 calibration and structural-discovery updates, preserving the locked failure."""
import json,os
from build_paper import P,H1,H2,BODY,EQ,tbl,fig,ROOT,PageBreak

def story_calibration(story,R):
    bg=R('bigg_model_survey.json'); census=R('bigg_moco_growth_summary.json');c=R('calibration_gate.json');a=R('calibration_lineage_audit.json'); second=R('bigg_second_family_verification.json'); rv=R('rousset_validation.json'); gd=R('rousset_genomic_dependence_audit.json'); mn=R('moco_specificity_null_ijn1463.json')
    assert census['scorable']==68 and census['models_with_at_least_one_flip']==63
    assert c['verdict']['first_passing_rung'] is None and len(c['ladder'])==5
    story += [PageBreak(),P('7. Later locked studies and audit of the calibration claim',H1),
              P('7.1 BiGG census: distinguishing model structure from biology',H2),
              P(f"The BiGG v2 listing supplied {bg['listed']} unique model accessions, each downloaded and hashed. "
                f"Of these, {len(bg['moco_biomass_models'])} have a positive-weight biomass objective consuming a MoCo-named "
                "compound; four snapshots without an annotated positive objective remain in the full denominator rather than being called negative. "
                "In the stored-medium knockout counterfactual, 63 of 68 scorable MoCo-positive models flip at least one gene call "
                "after the named MoCo biomass coefficients are removed; the reported binomial interval is 0.8367-0.9757. "
                "The total is 511 flipped model-gene calls, not 511 independent genes or wet-lab rescue experiments. "
                "The five zero-flip model IDs are retained in the machine-readable census. The positive-control iJO1366 restores "
                "the moa/mob cluster, including moaD b0784, under the exact counterfactual described in Section 4.3. "
                "See results/bigg_moco_growth_summary.json and per-model results/bigg_moco_growth/*.json.",BODY),
              P('This is a structural falsifiability result: the biomass objective forces a particular knockout prediction. '
                'A flip does not establish that the organism survives the knockout in its real medium. Closely related E. coli '
                'strain models share reaction templates, so neither a simple binomial confidence interval nor the 63/68 fraction '
                'can be read as 68 biologically independent replications. The result supports an audit of biomass-composition '
                'assumptions across reconstruction families, not the assertion that every objective is wrong. A separate '
                'manual second-family verification was subsequently completed on Pseudomonas putida iJN1463 (next paragraph).',BODY),
              P(f"Independent-code second-family check: in the downloaded Pseudomonas putida KT2440 iJN1463 snapshot, three negative MoCo objective coefficients "
                f"({', '.join(second['objective_negative_moco_coefficients'])}) were removed explicitly. The baseline wild-type optimum was {second['baseline_wt_growth']:.6f}/h; "
                f"the patched optimum was {second['patched_wt_growth']:.6f}/h. A fresh complete gene-deletion run found {second['baseline_essential_genes']} baseline-essential genes, "
                f"of which {len(second['flips_from_fresh_manual_counterfactual'])} flipped to growth at least 95% of the original wild type. "
                f"The independently calculated gene IDs are {', '.join(second['flips_from_fresh_manual_counterfactual'])}; each was separately zero-growth under its original model. "
                'The raw model snapshot SHA-256 was checked against the provider-index ledger before use; results/bigg_second_family_verification.json records coefficients, growth and gene IDs. '
                'This meets the narrow second-model-family structural criterion, but Pseudomonas viability in a lab was not tested; related model templates may still limit independence.',BODY),
              P('7.2 Locked decision gate and candidate ladder',H2),
              P('The follow-up gate was written down before scoring in notes/prereg_calibration_decision_gate.md. It asked '
                'whether the stacked gene-essentiality score crossed two tests on the 1,249 aligned genes: G1 AUROC >0.666 '
                'against iJO1366 minimal FBA, and G2 non-significant two-sided McNemar discordance (p >=0.05) with at least '
                'as many model wins as FBA wins. A1 locked a symmetric MCC-based threshold selector on the other two folds for '
                'both candidate and comparator, plus an ordered five-rung ladder. The earlier 27:50 discordant loss came from '
                'eval-inclusive F1-tuned thresholds, so A1 retains it as a harness reproduction, not as a clean baseline. '
                'The gate text stays fixed despite seeing a favorable-direction p-value.',BODY)]
    data=[['rung','AUROC','F1','model-only correct','FBA-only correct','McNemar 2-sided p','literal G2']]
    for r in c['ladder']:
        m=r['mcnemar_vs_comparator'];data.append([r['rung'],f"{r['auroc']:.4f}",f"{r['f1']:.4f}",str(m['b_ensemble_wins']),str(m['c_fba_wins']),f"{m['p_exact_2sided']:.5f}",'PASS' if m['g2_pass'] else 'FAIL'])
    story += tbl(data,'Table 21. Frozen G2 ladder. The significance clause p >=0.05 fails even when the sign favors the ensemble. G1 holds for every rung; no rung passes both clauses.')
    story += [P('R1 has 93 model-only correct versus 67 FBA-only correct calls, with exact p=0.04777. '
                'R2 has 98 versus 58, p=0.00170. Every rung gains more than it loses, but every rung has p<0.05, '
                'so the literal precommitted G2 FAILS. Changing it to a one-sided loss test after seeing results would '
                'misrepresent the preregistration. The observation can motivate a *fresh* study; it cannot retrospectively '
                'repair this one. The discrimination gain is real as a numerical OOF summary, but its inferential claims '
                'also need the lineage audit in the next section.',BODY),
              P('7.3 Information-flow audit: labels behind out-of-fold scores',H2),
              P(a['leak_path'],BODY), P(a['additional_leak_path'],BODY),
              P('The distinction is between the labels used by the threshold optimizer and those used upstream to train '
                'its score inputs. An OOF score for fold j was produced by a model trained on all genes outside j. '
                'When the outer evaluation fold k differs from j, those training genes include k. Thus the outer-fold '
                'labels can influence the score surface on which a threshold for k is selected. The base CNN, GNN and '
                'k-mer OOF predictions may themselves use label-trained complementary folds. Calling the full process '
                'strictly train-only is false, even though the final threshold optimizer takes only other-fold label vectors. '
                'A properly nested outer split must refit every supervised base model without k, generate inner-fold scores '
                'from only k-excluded training genes, choose threshold inside that pool, and score k once. A fresh untouched '
                'external cohort is preferable for generalization.',BODY)]
    rows=[['fold','train genes','eval genes','FBA threshold','calls vs fixed >0.5','MCC (fold)']]
    for row in a['fold_rows']:
        rows.append([str(row['fold']),str(row['n_train']),str(row['n_eval']),f"{row['threshold']:.6f}",str(row['threshold_vs_fixed_calls_differ']),f"{row['MCC_eval_fold_threshold']:.3f}"])
    story += tbl(rows,'Table 22. Comparator threshold audit on the same 1,249 genes. Numerical threshold values differ, but all pooled FBA calls match the structural >0.5 rule. This does not remove the candidate-score leakage.')
    story += [P('The fold-2 FBA threshold 0.110729 looks unstable beside 0.590300 on folds 0 and 1. Yet there are zero '
                'differences in comparator hard calls versus >0.5 over all three eval folds. This is an outcome of discrete '
                'model score gaps, not a license to choose a favorable threshold post hoc. The FBA comparator is physically '
                'defined by model growth fraction; the supervised stack and its thresholds are the principal leakage concern. '
                'We retain the literal saved gate result as a descriptive audit, mark the nested estimate as missing, and '
                'withhold a definitive hard-call benchmark-beat claim until its nested refit or independent cohort is complete.',BODY),
              P('7.4 Judge-informed research direction: where a virtual cell adds information',H2),
              P('The first recorded ChatGPT judge round was asked to critique the calibration paradox and propose a '
                'testable novelty analysis. Its critique identified fold provenance, leakage, comparator-threshold '
                'stability and multiplicity; its scientific suggestion was to ask whether the ensemble advantage '
                'over FBA concentrates in genes with strong transcriptional regulation but weak direct metabolic '
                'constraint. We wrote notes/prereg_discordance_regimes.md as an explicit *future* mechanism test. '
                'It requires a versioned E. coli transcription-factor-target map, a metabolic-constraint definition '
                'verified against code, matched/negative control partitions and an independent essentiality screen. '
                'STRING association degree is not a substitute for TF regulation. Existing Gerdes OOF outcomes '
                'cannot be recycled as confirmation after they inspired the hypothesis. This new protocol is a '
                'methodological novelty, not an observed biological discovery.',BODY),
              P('7.5 Independent assay as a different question',H2),
              P(f"The published Rousset et al. 2018 genome-wide CRISPRi assay supplies a biologically distinct but non-identical essentiality yardstick (https://doi.org/10.1371/journal.pgen.1007749). A separate pre-registration in notes/prereg_rousset_validation.md fixed median coding-strand depletion <= -5 as the binary outcome, UniProt gene-name mapping, and a ranking comparison against FBA. Of {rv['n_rousset_genes']} source rows, {rv['n_mapped']} map to locus IDs and {rv['n_evaluated']} overlap scored genes, including {rv['n_crispri_essential']} CRISPRi positives. The saved Gerdes-trained v2 score has AUROC {rv['auroc_v2_lr']:.4f} versus FBA {rv['auroc_fba_min']:.4f}; paired gene-bootstrap difference {rv['paired_v2_minus_fba']['mean']:+.4f} [95% CI {rv['paired_v2_minus_fba']['ci95'][0]:+.4f}, {rv['paired_v2_minus_fba']['ci95'][1]:+.4f}]. The study's locked numerical ranking bar was met, as detailed in results/rousset_validation.json and scripts/run_rousset_validation.py.",BODY),
              P('This is genuinely different assay labels, not a new organism, gene universe, or controlled independent training environment. CRISPRi knockdown can exert polar effects in bacterial operons, whereas transposon footprinting and FBA knockouts measure different perturbations. The external label assessment does NOT retroactively validate the defective Gerdes hard-call gate: it evaluates ranking, not a newly selected threshold or matched hard-call specificity. Some CRISPRi measurements may overlap the biology or published annotations used in external features; this needs source provenance review before a strong independence claim. Seven positives in the Gerdes-nonessential subset are too few to support a broad rescue statement despite its attractive AUROC. Report the distinct assay alongside, not in place of, the failed G2 and the MoCo mechanism.',BODY),
              P('7.6 Assay mapping and bias pathways',H2),
              P('The external assay was not used to tune the saved v2 score in the reported test, but mapping one primary gene name to one b-number can undercount aliases or collapse isoforms. Source rows without clean mapping are excluded by design, which changes both prevalence and gene set. A robust follow-up must show sensitivity to ambiguous-name mappings, operon-level grouped bootstrap, and the frozen CRISPRi median depletion threshold, not choose a more favorable cutoff after viewing AUROC. The confidence interval above resamples genes as if independent: adjacent genes in an operon and shared model features weaken that assumption. An operon-cluster bootstrap and separate Keio loss-of-function screen would better support generalization; they have not been run under this calibration gate. The biological test remains whether a virtual-cell output predicts condition-specific viability outside the training label source, not whether a ranking threshold can be adjusted until a gate passes.',BODY),
              P('7.7 Post-result genomic-neighborhood uncertainty sensitivity',H2),
              P('A gene-wise bootstrap treats nearby genes as independent draws despite operons and shared perturbation effects. To stress-test that arithmetic without redefining the CRISPRi gate, a separate post-result script groups the evaluated genes into fixed genomic b-number windows of width 5, 10, 20 or 50 IDs and resamples the occupied windows. This is a dependence sensitivity, not a true operon-cluster bootstrap: locus-number adjacency is an imperfect proxy, and windows can split operons or combine unrelated genes. Every window definition and random seed is retained in scripts/audit_rousset_genomic_dependence.py; the same 1,214 genes and fixed scores are used throughout. We did not refit a model, optimize a threshold, or use this descriptive result as a new registered discovery.',BODY),

              P(f"The paired AUROC difference on these {gd['n_genes']} mapped genes is {gd['score_difference']:+.4f}. Three thousand resamples at each genomic window width give intervals: 5 IDs [{gd['window_bootstrap'][0]['diff_ci95'][0]:+.4f}, {gd['window_bootstrap'][0]['diff_ci95'][1]:+.4f}], 10 IDs [{gd['window_bootstrap'][1]['diff_ci95'][0]:+.4f}, {gd['window_bootstrap'][1]['diff_ci95'][1]:+.4f}], 20 IDs [{gd['window_bootstrap'][2]['diff_ci95'][0]:+.4f}, {gd['window_bootstrap'][2]['diff_ci95'][1]:+.4f}], and 50 IDs [{gd['window_bootstrap'][3]['diff_ci95'][0]:+.4f}, {gd['window_bootstrap'][3]['diff_ci95'][1]:+.4f}]. All are positive in this chosen neighborhood surrogate, but none addresses annotation-source leakage, assay differences, batch effects, mapped-gene selection, or new-organism transfer. They strengthen only the limited claim that the saved ranking advantage is not driven solely by gene-wise resampling assumptions. The originally saved literal G2 failure and upstream score-lineage flaw stand unchanged (results/rousset_genomic_dependence_audit.json).",BODY),
              P('7.8 Is MoCo unusually disruptive? A one-model negative control',H2),
              P(f"A second ChatGPT role-play critique identified a missing specificity control: would the same seven-gene objective-dependent rescue appear after deleting other biomass ingredients? The exploratory post-result protocol in notes/prereg_moco_specificity_null.md was written before controls were computed but was not committed until afterward; it selected every one of {len(mn['controls'])} non-MoCo negative objective coefficients within a factor of two of at least one MoCo coefficient in the saved Pseudomonas iJN1463 model. The source JSON SHA-256 was checked again. Starting from {mn['baseline_essential_genes']} original zero-growth gene deletions, we removed each control coefficient one at a time and retested growth at >=95% original wild type under the same model and medium. The three-term MoCo bundle rescues {mn['positive_bundle']['n_baseline_essential_rescued']} knockout calls; {mn['n_controls_with_rescue_ge_moco_bundle']} of {len(mn['controls'])} single-term controls rescue at least as many. These are adenosylcobalamin (21 calls), biotinyl 5-AMP (8), pyoverdine (14) and thiamine diphosphate (8). Thirty-five controls rescue none. Every ID/count is retained in results/moco_specificity_null_ijn1463.json and the executable script.",BODY),
              P('This is a meaningful exploratory negative result for uniqueness, not a refutation of the MoCo-specific mechanistic gene IDs. Deleting even one alternative objective constituent can cause an equal or larger number of predicted essentiality flips. But comparison of one control coefficient with a three-ingredient MoCo bundle does not match perturbation size; coefficient proximity does not match molecule class, network connectivity or model frequency. Consequently we do not quote an exchangeability-based p-value or rank MoCo as a uniquely causal metabolite from this scout. The result redirects the claimed computational discovery from "MoCo uniquely breaks virtual-cell essentiality" to "biomass composition can make distinct gene calls sensitive to multiple small required terms." Multi-term category- and connectivity-matched controls across model families and condition-matched biological knockout evidence would be needed for a stronger claim. The old Gerdes G2 gate remains failed.',BODY),
              P('7.9 Revision of the headline and outstanding work',H2),

              P('The historical abstract reports the 0.7225 v1 AUROC and 27:50 hard-call loss, which is accurate for '
                'its original eval-inclusive threshold exercise but not the end of the project. The later MoCo census '
                'adds cross-reconstruction structural scope, and the v2 OOF AUROC is approximately 0.7967, yet a '
                'repeated-label score-lineage audit limits the calibration gate. Therefore the defensible headline '
                'today is: a modular model made testable predictions and exposed a biomass-objective artifact; the '
                'later hard-call improvement is suggestive, not a clean held-out win. An independent phenotypic '
                'screen and condition-matched carbon-source evaluation remain needed for the strongest claim.',BODY)]
    return story
