"""Paper content: Appendix H - head-to-head on the Bernstein et al. 2023 RB-TnSeq benchmark."""
import os, sys, json, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_paper import P, H1, H2, BODY, EQ, tbl, fig, ROOT, PageBreak
from reportlab.lib.units import inch

DESC = {'iML1515_base': 'iML1515 as released (their pipeline)',
        'iML1515_bernstein_vitamins': 'iML1515 + 5 hand-picked vitamins (Bernstein Part 6)',
        'iML1515_auto_cofactors': 'iML1515 + 19 audit-derived cofactors (ours, label-free)',
        'iML1515_bernstein_allcorr': 'INVALID: all-corrections file with its vitamin exchanges closed (our loader bug)',
        'iML1515_allcorr_plus_auto': 'all-corrections + audit-derived cofactors (ours)',
        'iML1515_bernstein_allcorr_fixed': 'Bernstein all-corrections model, their saved vitamins kept (published SOTA, reproduced)',
        'iML1515_allcorr_selected': 'INVALID: same bug; btn/thmpp/coa re-added the removed vitamins',
        'iJO1366_base': 'iJO1366 (their pipeline)', 'iJR904_base': 'iJR904 (their pipeline)', 'iAF1260_base': 'iAF1260 (their pipeline)'}


def story_bern(story, R):
    rows = [json.load(open(f)) for f in sorted(glob.glob(os.path.join(ROOT, 'results/bernstein/*.json'))) if not f.endswith('_ids.json')]
    rows.sort(key=lambda r: -r['pr_auc_bernstein_metric'])
    story += [PageBreak(), P('Appendix H. Head-to-head on the Bernstein et al. (2023) RB-TnSeq benchmark', H1),
              P('H.1 Why this benchmark', H2),
              P('Appendix E showed that Gerdes 2003 (rich LB medium) cannot test whether cofactor requirements in the '
                'biomass objective create false essentiality calls, because the model and the experiment use different '
                'media. Bernstein et al. (2023) built the condition-matched benchmark for exactly this question. It '
                'compares FBA growth calls with RB-TnSeq fitness (Price et al. 2018, E. coli BW25113 Keio_ML9 library) '
                'for every model gene on 25 defined carbon sources. They report that hand-supplying five vitamins '
                '(biotin, R-pantothenate, thiamin, tetrahydrofolate, NAD+) fixes many false negatives in iML1515, and '
                'they release an all-corrections model. Their fully corrected analysis (all-corrections model plus the five vitamins) is reported at PR-AUC 0.843 (0.844 with transaldolase constrained); that is the published state of the art for this task. The all-corrections model file alone, without vitamins, scores lower in their pipeline.', BODY),
              P('H.2 Protocol', H2),
              P('To rule out any protocol drift, we do not re-implement their benchmark. scripts/run_bernstein_benchmark.py '
                'executes the function cells of their released notebook (MIT licence, vendored under '
                'external/E_coli_GEM_validation) for model loading, gene and carbon-source matching, the BW25113 strain '
                'correction, removal of non-conditionally essential genes, carbon-source mapping, FBA simulation '
                '(carbon uptake -10, other media -1000 mmol/gDW/h, growth threshold 0.001) and the metric. We changed '
                'two things: their model_adjustments function reads a notebook-global variable, which we bind '
                'explicitly; and when a saved SBML model is loaded, exchanges that the file leaves open are re-opened after '
                'model_adjustments (which switches all exchanges off), exactly as their own all-corrections cell does. The metric is their precision-recall AUC, where labels are the model growth calls, the '
                'score is negative fitness and the positive class is no-growth:', BODY),
              P('AUC<sub>PR</sub> = &int; P(r) dr, &nbsp; P(t) = |{ f &lt; -t } &cap; {model no-growth}| / |{ f &lt; -t }|, &nbsp; R(t) = |{ f &lt; -t } &cap; {no-growth}| / |{no-growth}|    (H1)', EQ),
              P('Supplements are added as intracellular exchange reactions, the same construction as their Part 6. Our '
                'automated set is every cofactor/vitamin-class biomass compound that the iJO1366 provenance audit found '
                'unproducible in at least one lethal knockout. It was derived without any fitness data.', BODY)]
    t = [['variant', 'PR-AUC', 'ROC-AUC', 'no-growth calls', 'genes x carbons']]
    for r in rows:
        t.append([DESC.get(r['variant'], r['variant']), f"{r['pr_auc_bernstein_metric']:.3f}", f"{r['roc_auc']:.3f}",
                  str(r['n_model_nogrowth']), f"{r['n_genes']} x {r['n_carbon']}"])
    story += tbl(t, 'Table 19. Bernstein benchmark results (their code, their metric).',
                 widths=[3.0 * inch, 0.7 * inch, 0.7 * inch, 1.0 * inch, 1.1 * inch])
    fp = os.path.join(ROOT, 'figures', 'fig11_bernstein_benchmark.png')
    if os.path.exists(fp):
        story += fig(fp, 5.6 * inch, 'Figure 11. PR-AUC on the Bernstein benchmark by model variant.')
    story += [P('H.3 Reading', H2), P(R('bernstein_verdict.json')['text'], BODY)]
    story += [P('H.4 A reproduction error, kept on the record', H2),
              P('Our first run of the released all-corrections model scored 0.764. The cause was our loader: it closed every '
                'exchange in the saved file, and their model_adjustments also switches exchanges off, so the three vitamin '
                'uptakes the authors saved open (EX_btn_e, EX_thm_e, EX_pnto__R_e) were silently removed. Biotin, thiamin and '
                'pantothenate synthesis knockouts (b0775, b0776, b3990, b0133, b0134) then showed zero growth on every carbon '
                'source. A split-half search over our 19 audit-derived cofactors selected biotin, thiamin diphosphate and CoA '
                'and appeared to beat the model by +0.071 PR-AUC on held-out carbons (95% CI 0.037-0.111). That gain only '
                'restored the removed vitamins: after the fix the all-corrections model scores 0.839 on its own, identical to '
                'the "selected" run. The rows marked INVALID in Table 19 are kept so the error stays auditable. Lesson: an '
                'apparent benchmark win that re-discovers the benchmark authors\' own fix is a symptom of a broken baseline.', BODY)]
    v4 = R('ensemble_v4_esm.json')
    story += [P('H.5 Protein language model feature (pre-registered negative)', H2),
              P(f"ESM-2 (t6, 8M parameters; fair-esm) mean-pooled embeddings were computed for all 4,300 CDS of U00096.3. "
                f"Pre-registered test (notes/prereg_esm_v4.md): v2 features plus a 16-component PCA of the embedding fit inside "
                f"each training fold. AUROC {v4['v4_lr']['auroc']:.3f} vs {v4['v2_lr']['auroc']:.3f} for v2; paired bootstrap "
                f"difference {v4['paired_v4lr_vs_v2lr']['mean']:+.3f} (95% CI {v4['paired_v4lr_vs_v2lr']['ci95'][0]:+.3f} to "
                f"{v4['paired_v4lr_vs_v2lr']['ci95'][1]:+.3f}). The embedding alone reaches AUROC {v4['esm_only_lr']['auroc']:.3f}: "
                f"real signal, but redundant with the existing features. Verdict: negative.", BODY)]
    go = R('go_error_enrichment.json')
    story += [P('H.6 Where the ensemble is confidently wrong: GO enrichment', H2),
              P('Hard errors of the v2 ensemble were tested for Gene Ontology enrichment (goatools, propagated annotations '
                'from the EcoCyc GAF, universe = annotated benchmark genes) with the one-sided hypergeometric test and '
                'Benjamini-Hochberg correction:', BODY),
              P('p = &Sigma;<sub>i&ge;k</sub> C(K,i) C(N-K,n-i) / C(N,n), &nbsp; q<sub>(j)</sub> = min<sub>l&ge;j</sub> ( m p<sub>(l)</sub> / l )    (H2)', EQ)]
    t2 = [['error set', 'GO term', 'study', 'universe', 'BH q']]
    for x in go['false_positives'][:8]:
        t2.append(['confident FP', x['name'][:62], x['study'], x['pop'], f"{x['p_fdr_bh']:.1e}"])
    story += tbl(t2, f"Table 20. GO terms enriched among confident false positives (non-essential genes in the top 10% of v2 scores; n={go['n_fp']}). False negatives (n={go['n_fn']}): no term at q<0.05.",
                 widths=[1.0 * inch, 3.0 * inch, 0.7 * inch, 0.8 * inch, 0.8 * inch])
    story += [P('Molybdopterin cofactor biosynthesis and the ATP synthase are required only anaerobically or are bypassed by '
                'fermentation in rich medium; the ensemble inherits a condition mismatch from its FBA inputs, the same class '
                'of error Bernstein et al. fixed with vitamins. This is a hypothesis for condition-aware features, not a claim.', BODY)]
    v5 = R('ensemble_v5_go.json')
    story += [P(f"A nested test turned this into a feature: GO terms enriched among training-fold confident false positives (selected "
                f"inside each training fold only; notes/prereg_go_condition_feature.md) became a binary flag. AUROC {v5['v5_lr']['auroc']:.3f} vs "
                f"{v5['v2_lr']['auroc']:.3f} (difference {v5['paired_v5_vs_v2']['mean']:+.3f}, 95% CI {v5['paired_v5_vs_v2']['ci95'][0]:+.4f} to "
                f"{v5['paired_v5_vs_v2']['ci95'][1]:+.4f}): a small, significant loss. Negative.", BODY)]
    pc = R('pathway_concordance.json'); a = pc['primary_all_rescued']; w = pc['secondary_within_vitamins']; n = pc['secondary_within_nonvitamin']
    story += [P('H.7 Finding: the supplement-bypass artifact (pre-registered)', H2),
              P('Every model no-growth call that a supplement turns into growth is a candidate correction. We asked which of these '
                'rescues the experiment agrees with. Hypothesis (notes/prereg_pathway_concordance.md, committed before the split was '
                'computed): a rescue is real when the knocked-out gene lies in the supplement\'s own KEGG biosynthesis pathway, and '
                'spurious when the model rescues an off-pathway gene by metabolising the supplement as a nutrient or one-carbon source. '
                'With f the RB-TnSeq fitness of a rescued gene-carbon pair and on(g,s) = 1 if gene g is in a KEGG pathway of supplement s:', BODY),
              P('&Delta; = median{ f : on = 1 } - median{ f : on = 0 }, &nbsp; H<sub>0</sub>: &Delta; &le; 0, &nbsp; CI by resampling genes    (H3)', EQ)]
    t3 = [['set', 'on: pairs (genes)', 'off: pairs (genes)', 'median f', 'f > -2', 'Delta [95% CI]', 'p']]
    for lab, x in (('all rescues (primary)', a), ('vitamins only', w), ('SAM/PLP/folates only', n)):
        t3.append([lab, f"{x['n_on_pairs']} ({x['n_on_genes']})", f"{x['n_off_pairs']} ({x['n_off_genes']})",
                   f"{x['median_on']:.2f} / {x['median_off']:.2f}", f"{x['frac_gt_-2_on']:.2f} / {x['frac_gt_-2_off']:.2f}",
                   f"{x['median_diff']:.2f} [{x['cluster_boot_ci95'][0]:.2f}, {x['cluster_boot_ci95'][1]:.2f}]", f"{x['mwu_p_one_sided']:.0e}"])
    story += tbl(t3, 'Table 21. Pathway concordance of supplement rescue on the Bernstein benchmark. median f and f > -2 are given as on / off; one-sided Mann-Whitney; gene-cluster bootstrap CI.',
                 widths=[1.3 * inch, 1.0 * inch, 1.0 * inch, 0.9 * inch, 0.8 * inch, 1.0 * inch, 0.5 * inch])
    story += [P('The primary criterion is met: on-pathway rescues sit near neutral fitness, off-pathway rescues are almost all truly '
                'deleterious. Caveats: the primary contrast is confounded with supplement identity (vitamin rescues are almost all '
                'on-pathway; SAM and folate rescues almost all off-pathway, e.g. succinate dehydrogenase sdhCDAB and galactose genes '
                'rescued by SAM). The within-group contrasts agree in direction but rest on 1 off-pathway vitamin gene and 5 on-pathway '
                'non-vitamin genes. One organism, one benchmark. A literature search found no prior test of this concordance '
                '(notes/novelty_check.md); novelty is probable, not proven. Falsifier: a supplement whose off-pathway rescues have '
                'near-neutral fitness in an independent fitness dataset. Practical consequence: a supplement correction should be '
                'accepted only for on-pathway rescues; this rule would have rejected our own 19-compound blanket set (PR-AUC 0.548).', BODY)]

    mr = R('cross_species_mr1.json'); pu = R('cross_species_putida.json')
    story += [P('H.8 Cross-species replication (pre-registered; underpowered / not significant)', H2),
              P('The protocol was fixed in notes/prereg_cross_species.md (git 00e7a9d, before any MR-1 fitness value was compared; '
                'amendment 1 adding P. putida in git 999b51c, before that analysis ran). Each organism uses its own genome-scale model on '
                'a vitamin-free minimal medium, the same seven cofactor supplements supplied as cytosolic sinks with vcell rescue-audit, '
                'KEGG gene-pathway links for that organism, and RB-TnSeq fitness from the Fitness Browser (Price et al. 2018). A gene is '
                'on-pathway if any supplement that rescues it shares a KEGG pathway with it. The outcome is the median log2 fitness of the '
                'gene over all experiments; criteria as in (H3) at gene level, with fewer than 5 genes per group declared underpowered.', BODY)]
    t4 = [['organism', 'WT', 'pairs', 'genes', 'w/ fitness', 'median f', 'p', 'verdict']]
    for name, x, gk in (('MR-1, iMR1_799, lac', mr, 'wt_growth_lactate'), ('KT2440, iJN1463, glc', pu, 'wt_growth_glucose')):
        mo = f"{x['median_on']:.2f} / {x['median_off']:.2f}" if 'median_on' in x else 'n/a'
        pv = f"{x['mwu_p_one_sided']:.2f}" if 'mwu_p_one_sided' in x else 'n/a'
        t4.append([name, f"{x[gk]:.2f}", str(x['n_rescue_pairs']), f"{x['genes_on']} / {x['genes_off']}",
                   f"{x['with_fitness_on']} / {x['with_fitness_off']}", mo, pv, x['verdict'].split(' (')[0]])
    cm = R('cross_species_carveme.json')
    NAMES = {'Caulo': 'Caulobacter', 'Cola': 'Echinicola', 'PS': 'Dechlorosoma', 'Ponti': 'Pontibacter',
             'Smeli': 'S. meliloti', 'SyringaeB728a': 'P. syringae'}
    for o, x in cm['per_organism'].items():
        mo = f"{x['median_on']:.2f} / {x['median_off']:.2f}" if 'median_on' in x else 'n/a'
        pv = f"{x['mwu_p']:.2f}" if 'mwu_p' in x else 'n/a'
        t4.append([NAMES[o] + ', CM, glc', f"{x['wt']:.2f}", str(x['n_rescue_pairs']), f"{x['genes_on']} / {x['genes_off']}",
                   f"{x['with_fitness_on']} / {x['with_fitness_off']}", mo, pv, x['verdict'].split(' (')[0]])
    pp, ns = cm['primary_pooled_percentile'], cm['secondary_pooled_excluding_SAM']
    for lab, x in (('pooled CM (primary)', pp), ('pooled CM, no SAM', ns)):
        t4.append([lab, '', '', '', f"{x['n_on']} / {x['n_off']}", f"pct {x['median_pct_on']:.2f} / {x['median_pct_off']:.2f}", f"{x['mwu_p']:.2f}", x['verdict'].split(' (')[0]])
    story += tbl(t4, 'Table 22. Cross-species test of the supplement-bypass finding. genes, w/ fitness and median f are on / off; WT = wild-type growth (1/h); lac/glc = lactate/glucose medium; CM = CarveMe-built model; pct = within-organism fitness percentile; one-sided Mann-Whitney on gene medians (results/cross_species_mr1.json, results/cross_species_putida.json, results/cross_species_carveme.json).',
                 widths=[1.45 * inch, 0.4 * inch, 0.45 * inch, 0.6 * inch, 0.7 * inch, 0.95 * inch, 0.4 * inch, 1.35 * inch])
    story += [P(f"In MR-1 every rescue was on-pathway, so the test could not run (underpowered, no claim). The PSAMM export of "
                "iMR1_799 also left three multi-compound pseudo-exchanges (casamino acids, gelatin, Tween 20) open in both directions, "
                "which let the model grow with no carbon source; we closed them and report it as a model defect. In P. putida the "
                f"difference points the predicted way (+{pu['median_diff']:.2f} log2 units; {100*pu['absent_frac_off']:.0f}% of off-pathway "
                f"genes lack any fitness value, often a mark of essentiality, versus {100*pu['absent_frac_on']:.0f}% on-pathway) but p = "
                f"{pu['mwu_p_one_sided']:.2f} fails the pre-registered bar. All 12 off-pathway P. putida rescues come from SAM, the same "
                "supplement behind most E. coli off-pathway rescues, so the supplement-identity confound is not broken by the second "
                "species. iJN1463 has no free thiamine metabolite, so thiamine was not tested there. "
                f"Amendments 3-4 widened the test to six Fitness Browser organisms with models built by CarveMe (DIAMOND, SCIP). The "
                f"inclusion rule for KEGG identifier matching was mis-specified as written (it excluded every organism by construction) and was "
                f"corrected, with disclosure, before any model was built. The pooled pre-registered primary on within-organism fitness "
                f"percentiles is not significant (difference {pp['diff']:+.2f}, p = {pp['mwu_p']:.2f}), nor is the test without SAM-only "
                f"rescues (p = {ns['mwu_p']:.2f}); an exploratory count goes against the hypothesis ("
                f"{100*cm['exploratory_absent_fraction']['off'][2]:.0f}% of off-pathway vs {100*cm['exploratory_absent_fraction']['on'][2]:.0f}% of "
                f"on-pathway genes lack fitness values). Honest status: strongly supported in E. coli, not replicated in eight other "
                f"bacteria (same direction, not significant). We therefore state the finding as E. coli-specific, most likely tied to how "
                f"iML1515 metabolises SAM, and not as a general law.", BODY)]
    story += [P('H.9 Tool: vcell rescue-audit', H2),
              P('The finding is packaged as a command that runs on any COBRA model: python -m vcell rescue-audit MODEL --supplement '
                'btn_c=00780 thf_c=00790,00670 ... --gene-pathways kegg_links.tsv --out rescues.csv. It knocks out each gene, supplies '
                'each supplement through a sink, and labels every rescue on_pathway or off_pathway (flag likely_artifact). It is covered '
                'by a hermetic test on the E. coli core model (icd rescued by 2-oxoglutarate is on-pathway; enolase rescued by pyruvate '
                'is off-pathway). On iML1515 with glucose minimal medium it reports 44 rescues, 7 off-pathway: SAM rescuing five purine '
                'biosynthesis genes (purK, purE, purB, purA, purC), biotin rescuing fabH and thiamine rescuing b4407 '
                '(results/rescue_audit_iML1515_M9.csv). Neither COBRApy nor MEMOTE, the tools used here, labels '
                'supplement rescues by pathway concordance.', BODY)]
    sm = R('sam_bypass_mechanism.json'); rg = sm['rescued_growth']['b0523']
    story += [P('H.10 Mechanism of the SAM bypass (exploratory, not pre-registered)', H2),
              P(f"SAM is the supplement behind most off-pathway rescues in E. coli and all of them in P. putida. In iML1515 on glucose "
                f"minimal medium SAM rescues knockouts of five purine biosynthesis genes (purK, purE, purB, purA, purC). Removing single "
                f"reactions shows the route (results/sam_bypass_mechanism.json): every one of the five rescues needs SAH nucleosidase "
                f"(AHCYSNS, mtnN/b0159; rescued growth {rg['none']:.2f} to {rg['AHCYSNS (SAH nucleosidase)']:.2f} per h). SAM is "
                f"methyl-transferred to S-adenosylhomocysteine, the nucleosidase releases adenine (Rhea RHEA:17805; UniProt P0AF12), and "
                f"adenine enters purine salvage. The rescued mutant grows at {rg['none']:.2f} per h against {sm['wt_growth']:.2f} for the "
                f"wild type: with homocysteine S-methyltransferase (HCYSMT) the model uses the supplied SAM as a bulk carbon and purine "
                f"source; without it the rescue remains at {rg['HCYSMT (homocysteine S-methyltransferase)']:.2f} per h. Two practical "
                "rules follow and are built into vcell rescue-audit: flag off-pathway rescues, and flag any rescue whose growth exceeds the "
                "unmodified model (column exceeds_wt), because a cofactor supplied at trace need cannot raise growth above wild type.", BODY)]
    rv = R('rousset_validation.json')
    story += [P('H.11 External validation on genome-wide CRISPRi (pre-registered)', H2),
              P(f"To test whether the v2 model generalises beyond the Gerdes 2003 transposon labels it was trained on, we scored its "
                f"out-of-fold predictions against an independent assay: the genome-wide CRISPRi screen of Rousset et al. 2018 (S12 table, "
                f"gene-level median log2 fold-change of coding-strand guides). The test was written down before scoring "
                f"(notes/prereg_rousset_validation.md): label CRISPRi-essential if median log2FC <= -5; primary metric AUROC of v2 LR with a "
                f"paired 2000-sample gene bootstrap against FBA. Of {rv['n_rousset_genes']} CRISPRi genes, {rv['n_mapped']} mapped to "
                f"b-numbers and {rv['n_evaluated']} overlap the model gene set ({rv['n_crispri_essential']} CRISPRi-essential). v2 LR reaches "
                f"AUROC {rv['auroc_v2_lr']:.3f} against {rv['auroc_fba_min']:.3f} for FBA (paired difference "
                f"{rv['paired_v2_minus_fba']['mean']:.3f}, 95% CI {rv['paired_v2_minus_fba']['ci95'][0]:.3f} to "
                f"{rv['paired_v2_minus_fba']['ci95'][1]:.3f}); the pre-registered verdict is {rv['verdict']} "
                f"(results/rousset_validation.json). Spearman between score and depletion is {rv['spearman_v2_vs_depletion']:.2f}. "
                f"Limits: the labels overlap heavily with Gerdes (only {rv['gerdes_nonessential_subset']['n_crispri_essential']} "
                f"CRISPRi-essential genes are Gerdes-nonessential, so the AUROC of {rv['gerdes_nonessential_subset']['auroc_v2_lr']:.3f} on that "
                "subset rests on very few positives), and CRISPRi is polar within operons. This is agreement across assays, not a new "
                "state of the art.", BODY)]
    v6 = R('v6_structure_domain.json')
    story += [P('H.12 Structure confidence and Pfam families add nothing (pre-registered negative)', H2),
              P(f"We tested whether predicted protein order and domain family membership carry essentiality signal beyond v2 "
                f"(notes/prereg_structure_domain_features.md, committed before scoring). Features: AlphaFold DB mean pLDDT and fraction "
                f"of very-low pLDDT residues ({v6['af_coverage']} of {v6['n_genes']} genes), number of Pfam families, and a fold-internal "
                f"Pfam essentiality rate computed only from training-fold genes. On their own these features reach AUROC "
                f"{v6['auroc']['new_features_only']:.3f} (univariate pLDDT {v6['univariate_auroc']['plddt_mean']:.3f}). Added to v2 on the "
                f"same folds, AUROC moves from {v6['auroc']['v2_lr']:.3f} to {v6['auroc']['v6_lr']:.3f} (paired difference "
                f"{v6['primary_auroc_v6_minus_v2']['diff']:+.4f}, 95% CI {v6['primary_auroc_v6_minus_v2']['ci95'][0]:+.4f} to "
                f"{v6['primary_auroc_v6_minus_v2']['ci95'][1]:+.4f}); AUPRC and the Rousset CRISPRi check are also flat. Verdict: "
                f"{v6['verdict']} (results/v6_structure_domain.json). Nearly all E. coli proteins in the model are well folded "
                f"(median pLDDT near 95), so order does not separate essential from non-essential enzymes.", BODY)]
    v7 = R('v7_oma.json')
    story += [P('H.13 Evolutionary conservation adds nothing either (pre-registered negative)', H2),
              P(f"Conserved genes are often said to be essential. We tested this against v2 using the depth of each gene's hierarchical "
                f"orthologous group in the OMA Browser (number of taxonomic levels, from E. coli up to LUCA; notes/prereg_oma_conservation.md, "
                f"with amendment 1 made before scoring because the 1:1-ortholog endpoint was too slow). {v7['oma_nonzero']} of "
                f"{v7['n_genes']} genes have a HOG; most of the model's genes are ancient (356 reach LUCA, 702 Bacteria). Alone the "
                f"feature reaches AUROC {v7['univariate_auroc']['log_orth']:.3f}; added to v2 the AUROC changes by "
                f"{v7['primary_auroc_v7_minus_v2']['diff']:+.4f} (95% CI {v7['primary_auroc_v7_minus_v2']['ci95'][0]:+.4f} to "
                f"{v7['primary_auroc_v7_minus_v2']['ci95'][1]:+.4f}). Verdict: {v7['verdict']} (results/v7_oma.json). Within metabolic "
                "genes, conservation carries weak signal that v2 already captures through its network and sequence features.", BODY)]
    v8 = R('v8_cog_pdb.json')
    story += [P('H.14 COG categories, phyletic spread and PDB coverage (pre-registered null with a positive trend)', H2),
              P(f"The third feature test added NCBI COG 2020 functional-category indicators, the number of COG genomes carrying each "
                f"gene's COG, and the number of PDB entries (notes/prereg_cog_pdb.md). {v8['cog_coverage']} of {v8['n_genes']} genes have a "
                f"COG and {v8['pdb_nonzero']} have a PDB entry. These features are informative on their own (AUROC "
                f"{v8['auroc']['new_features_only']:.3f}; COG spread alone {v8['univariate_auroc']['log_cog_spread']:.3f}), unlike pLDDT or "
                f"OMA depth. Added to v2, AUROC rises from {v8['auroc']['v2_lr']:.3f} to {v8['auroc']['v8_lr']:.3f}, but the paired 95% CI "
                f"({v8['primary_auroc_v8_minus_v2']['ci95'][0]:+.4f} to {v8['primary_auroc_v8_minus_v2']['ci95'][1]:+.4f}) includes zero, so the "
                f"pre-registered verdict is {v8['verdict']} (results/v8_cog_pdb.json). Dropping the PDB count, the leakage-prone term, "
                f"leaves the gain unchanged ({v8['secondary_v8_noPDB_minus_v2']['diff']:+.4f}). The v2 model already absorbs most of what "
                "conservation and function categories say about metabolic-gene essentiality.", BODY)]
    v9 = R('v9_precise1k.json')
    story += [P('H.15 Transcriptome features: primary null, positive secondaries (pre-registered)', H2),
              P(f"The fourth feature test used the PRECISE-1K E. coli RNA-seq compendium ({v9['n_samples']} QC-passed samples) and its "
                f"iModulon decomposition (notes/prereg_precise1k.md): mean log-TPM, SD of log-TPM across conditions and the number of "
                f"iModulons a gene belongs to. Essential genes are highly expressed (univariate AUROC {v9['univariate_auroc']['expr_mean']:.3f}), "
                f"stable across conditions (SD AUROC {v9['univariate_auroc']['expr_sd']:.3f}, i.e. lower variability predicts essentiality) "
                f"and belong to fewer iModulons ({v9['univariate_auroc']['n_imod']:.3f}). Together the three features reach "
                f"AUROC {v9['auroc']['new_features_only']:.3f} alone, the strongest external block tested. Added to v2, the pre-registered primary "
                f"AUROC gain is {v9['primary_auroc_v9_minus_v2']['diff']:+.4f} (95% CI {v9['primary_auroc_v9_minus_v2']['ci95'][0]:+.4f} to "
                f"{v9['primary_auroc_v9_minus_v2']['ci95'][1]:+.4f}): {v9['verdict']}. Two pre-registered secondaries exclude zero: AUPRC "
                f"{v9['secondary_auprc_v9_minus_v2']['diff']:+.4f} ({v9['secondary_auprc_v9_minus_v2']['ci95'][0]:+.4f} to "
                f"{v9['secondary_auprc_v9_minus_v2']['ci95'][1]:+.4f}) and AUROC on the independent Rousset CRISPRi labels "
                f"{v9['secondary_rousset_auroc_v9_minus_v2']['diff']:+.4f} ({v9['secondary_rousset_auroc_v9_minus_v2']['ci95'][0]:+.4f} to "
                f"{v9['secondary_rousset_auroc_v9_minus_v2']['ci95'][1]:+.4f}). These are secondary endpoints without multiplicity correction, "
                "so we report them as supporting evidence that expression carries signal v2 lacks, not as a confirmed improvement. "
                "One iModulon name differed between files (Superoxide in M.csv, SoxS in the thresholds) and was matched by hand "
                "(results/v9_precise1k.json).", BODY)]
    up = R('uptake_plausibility.json'); ps = up['per_supplement']
    story += [P('H.16 Can E. coli take the rescuing supplements up? (exploratory)', H2),
              P(f"A supplement rescue in the model is only meaningful if the cell can import the compound. We asked, for the seven "
                f"supplements in the iML1515 rescue audit, whether E. coli K-12 has a transporter in the Transporter Classification "
                f"Database whose curated substrates include the compound (notes/prereg_uptake_plausibility.md; ChEBI IDs expanded "
                f"through the ChEBI ontology). Amendment 1, disclosed in the note, added secondary ChEBI IDs after the first run missed "
                f"ThiBPQ because TCDB uses CHEBI:9530 for thiamine. Result (results/uptake_plausibility.json): only pantothenate (PanF, "
                f"TC 2.A.21.1.1) and thiamine (ThiBPQ, TC 3.A.1.19.1) have a K-12 uptake system. SAM, NAD, tetrahydrofolate and "
                f"pyridoxal 5'-phosphate have none, as expected in advance; biotin was expected to be supported but has no K-12 entry in "
                f"TCDB (its E. coli uptake route is not characterised there). {up['off_pathway_rescues_without_uptake']} of "
                f"{up['off_pathway_rescues_total']} off-pathway rescues use a supplement with no known K-12 uptake system, including all "
                f"{ps['amet_c']['rescues'].get('off_pathway', 0)} SAM rescues. This is independent support that the SAM rescues are artifacts "
                "of supplying an intracellular metabolite. With seven supplements it is descriptive, not a test, and absence from TCDB is "
                "weaker evidence than presence.", BODY)]
    x5 = R('cross_species_xs5.json'); ec = R('xs5_ecoli_consistency.json')
    story += [P('H.17 Forty-four organisms: the model-internal label is degenerate (pre-registered, UNDERPOWERED)', H2),
              P(f"Amendment 5 replaced the KEGG pathway label with a model-internal one (a rescue is on-pathway if the knockout abolishes "
                f"production of the supplement) so that every Fitness Browser bacterium could be tested. We built CarveMe models for all "
                f"{x5['n_organisms_analysed']} eligible organisms (none excluded; results/xs5/). The label did not separate anything: all "
                f"{x5['primary_pooled_percentile']['n_on']} rescued genes with fitness data were labelled on-pathway and none off-pathway, in "
                f"every organism, so the pre-registered primary and the no-SAM secondary are both {x5['primary_pooled_percentile']['verdict']} by "
                f"their own rule (results/cross_species_xs5.json). The E. coli consistency check shows why ({ec['n_rows_labelled']} of "
                f"{ec['n_rows_H7']} H.7 rows labelled, all on-pathway; results/xs5_ecoli_consistency.json): a lethal knockout that a supplement "
                "rescues nearly always also blocks de-novo synthesis of that supplement, because the supplement shares precursors with the "
                "blocked product (SAM and NAD contain adenosine, so every purine gene is needed to make them). We documented this in "
                "amendment 6 before computing the primary. The negative is informative for method design: production-blocked labels cannot "
                "stand in for pathway membership, and the E. coli finding (H.7) remains tested only with the KEGG label in eight organisms (H.8).", BODY)]
    x7 = R('cross_species_xs7.json'); pa = x7['primary_pooled_all']; ns = x7['secondary_pooled_all_excluding_SAM']
    story += [P('H.18 KEGG-label test extended (amendment 7): not significant, and falsified without SAM', H2),
              P(f"Because the model-internal label failed, amendment 7 extended the original KEGG-label test to every amendment-5 organism "
                f"whose KEGG gene identifiers match Fitness Browser locus ids (results/xs7_kegg_discovery.json). Only P. putida qualified "
                f"(76% of KEGG pathway genes match); the other organisms use numeric or RefSeq locus ids that KEGG does not index, and "
                f"Synechococcus fell below the 50% rule. Pooled with the six amendment-3 organisms ({len(x7['organisms_in_pool'])} organisms, "
                f"{pa['n_on']} on- vs {pa['n_off']} off-pathway genes), on-pathway genes sit at the {100*pa['median_pct_on']:.0f}th fitness "
                f"percentile against {100*pa['median_pct_off']:.0f}th for off-pathway genes (p = {pa['mwu_p']:.2f}): {pa['verdict']}. Without "
                f"SAM rescues the difference is {ns['diff']:+.3f} (p = {ns['mwu_p']:.2f}): {ns['verdict']} (results/cross_species_xs7.json). "
                "Outside E. coli the supplement-bypass signal is carried, if at all, by SAM alone.", BODY)]
    px = R('paxdb_datasets.json')
    story += [P('H.19 The abundance signal holds in every proteomics dataset (pre-registered)', H2),
              P(f"Integrated PaxDb abundance is a v3 feature. To check it is not an artefact of integration we scored each of the "
                f"{px['n_datasets']} individual E. coli K-12 datasets in PaxDb v5 separately (notes/prereg_paxdb_datasets.md). Highly "
                f"abundant proteins are more often essential in all {px['n_auroc_gt_0.5']} of {px['n_scored']} datasets (sign test "
                f"p = {px['sign_test_p']:.1e}; AUROC range {px['auroc_range'][0]:.3f} to {px['auroc_range'][1]:.3f}; inverse-variance pooled "
                f"AUROC {px['pooled_auroc_ivw']:.3f}, 95% CI {px['pooled_ci95'][0]:.3f} to {px['pooled_ci95'][1]:.3f}): "
                f"{px['verdict']} (results/paxdb_datasets.json); AUROC rises with coverage (Spearman {px['spearman_auroc_vs_coverage']:.2f}). "
                "A parsing error (a fourth column in some files) first left 11 datasets empty; it was fixed before interpretation.", BODY)]
    pc = R('paxdb_cross.json'); pp = pc['primary']; po = pc['per_organism']
    story += [P('H.20 Abundance predicts essentiality in five other bacteria (pre-registered)', H2),
              P(f"H.19 is E. coli only. We matched {pp['n_datasets']} PaxDb v5 datasets from five Fitness Browser organisms "
                "(P. putida, S. elongatus, S. oneidensis, B. thetaiotaomicron, D. vulgaris) to RB-TnSeq data by exact protein sequence "
                "via STRING v12 (notes/prereg_paxdb_cross.md). Genes with no fitness row were labelled putatively essential. "
                f"Abundant proteins were more often essential in {pp['n_auroc_above_half']} of {pp['n_datasets']} datasets "
                f"(sign test p = {pp['sign_test_p']:.3f}): {pp['verdict']}. Pooled AUROC per organism ranged from "
                f"{min(v['pooled_auroc'] for v in po.values()):.2f} (S. elongatus) to {max(v['pooled_auroc'] for v in po.values()):.2f} "
                f"(B. thetaiotaomicron). The abundance effect stayed positive after adjusting for protein length in "
                f"{pc['secondary']['n_abundance_coef_positive_after_length']} of {pc['secondary']['n']} datasets. Negative: among non-essential "
                f"genes, abundance was negatively correlated with median fitness in only {pc['secondary']['n_spearman_negative']} of "
                f"{pc['secondary']['n']} datasets, so the signal separates essential from non-essential genes but does not grade fitness "
                "(results/paxdb_cross.json). The label is a proxy: short genes can lack fitness rows for technical reasons.", BODY)]
    pm = R('paxdb_mtb.json'); q = pm['primary']
    story += [P('H.21 The abundance signal holds in M. tuberculosis (pre-registered)', H2),
              P(f"We scored all 16 PaxDb v5 M. tuberculosis H37Rv datasets against the DeJesus et al. (2017) saturating Tn-seq calls "
                f"({pm['label_counts']['essential']} essential, {pm['label_counts']['nonessential']} non-essential; notes/prereg_paxdb_mtb.md). "
                f"Two files were re-releases of the PeptideAtlas build and were merged by the pre-registered duplicate rule, leaving {q['n_unique']}. "
                f"Abundant proteins were more often essential in {q['n_auroc_above_half']} of {q['n_unique']} (sign test p = {q['sign_test_p']:.4f}; "
                f"pooled AUROC {q['pooled_auroc']:.3f}, 95% CI {q['ci95'][0]:.3f} to {q['ci95'][1]:.3f}): {q['verdict']} (results/paxdb_mtb.json). "
                f"The one exception is a 201-gene dataset (AUROC {q['auroc_range'][0]:.2f}). Growth-defect genes sat between essential and "
                f"non-essential genes in median abundance in {pm['secondary_gd_between']} of {q['n_unique']} datasets. Together with H.19 and H.20 "
                "the abundance-essentiality link holds across seven bacterial species from four phyla, which makes it a general feature "
                "rather than an E. coli quirk; it is, however, a known correlate, and what is new here is its per-dataset robustness, not its existence.", BODY)]
    pb = R('paxdb_bsub.json'); qb = pb['primary']; sc = R('datasets_strict_count.json')
    story += [P('H.22 B. subtilis, and how datasets are counted (pre-registered)', H2),
              P(f"Against the {pb['n_essential_list']} SubtiWiki essential genes, abundant proteins were more often essential in "
                f"{qb['n_auroc_above_half']} of {qb['n']} B. subtilis PaxDb datasets (sign test p = {qb['sign_test_p']:.3f}, the smallest "
                f"possible with five; pooled AUROC {qb['pooled_auroc']:.3f}, 95% CI {qb['ci95'][0]:.3f} to {qb['ci95'][1]:.3f}): {qb['verdict']} "
                "(results/paxdb_bsub.json; notes/prereg_paxdb_bsub.md). "
                f"Dataset counting: the manifest lists {sc['n_manifest_entries']} accession-level entries. Counting one per source study "
                "(re-deposits, quantification variants and fractions of one study collapse; derived tables not counted) gives "
                f"{sc['n_strict']}; that strict number is the one we report. The rule and the collapsed groups are in "
                "results/datasets_strict_count.json.", BODY)]
    pa = R('paxdb_pao1.json'); qa = pa['primary']
    story += [P('H.23 P. aeruginosa (pre-registered)', H2),
              P(f"Against the {pa['n_core_named']} named core essential genes of Poulsen et al. (2019), matched to PAO1 by gene name, "
                f"abundant proteins were more often essential in {qa['n_auroc_above_half']} of {qa['n_studies']} PaxDb studies (files from one "
                f"study averaged first; sign test p = {qa['sign_test_p']:.4f}; study AUROC {qa['auroc_range'][0]:.2f} to {qa['auroc_range'][1]:.2f}): "
                f"{qa['verdict']} (results/paxdb_pao1.json; notes/prereg_paxdb_pao1.md). PAO1 was not among the nine strains Poulsen screened "
                "and unnamed genes were dropped, so this label favours well-studied genes. Nine species now show the link.", BODY)]
    pg = R('paxdb_mtb_griffin.json'); qg = pg['primary']
    story += [P('H.24 The Mtb result survives an independent essentiality label (pre-registered)', H2),
              P(f"Re-scoring the {qg['n']} Mtb datasets of H.21 against Griffin et al. (2011) Tn-seq calls ({pg['n_griffin_essential']} essential of "
                f"{pg['n_griffin_labelled']}; agreement with DeJesus kappa = {pg['label_agreement']['cohen_kappa']:.2f}) gives AUROC > 0.5 in "
                f"{qg['n_auroc_above_half']} of {qg['n']} (sign test p = {qg['sign_test_p']:.4f}, median AUROC {qg['median_auroc']:.3f}): {qg['verdict']} "
                "(results/paxdb_mtb_griffin.json). AUROCs are lower than with DeJesus labels, as expected from the older, less saturated screen; "
                "the same 201-gene dataset is again the exception.", BODY)]
    return story
