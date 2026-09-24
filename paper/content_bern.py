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
    return story
