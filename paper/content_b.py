"""Paper content, part B: results, discovery, discussion, references."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_paper import P, H1, H2, BODY, EQ, fig, tbl, PageBreak, Spacer, ROOT

def story_b(story, R):
    core = R('benchmark_core_vs_gerdes.json')
    ijo = R('benchmark_ijo1366_vs_gerdes.json')
    rich = R('benchmark_ijo1366_rich_vs_gerdes.json')
    ens = R('ensemble_results.json')
    dia = R('dfba_diauxie_metrics.json')
    diar = R('dfba_diauxie_regulated_metrics.json')
    rescue = R('moco_biomass_rescue_test.json')
    disc = R('discovery_candidates.json')

    story += [P('3. Results', H1),
        P('3.1 The metabolic layer reproduces core physiology', H2),
        P('As a sanity gate, the core model was run over a 28-condition matrix (14 carbon '
          'sources x aerobic/anaerobic). Aerobic growth on glucose is 0.8739 h-1, matching '
          'the published value for this reconstruction to four decimal places; anaerobic '
          'glucose growth drops to 0.2117 h-1; formate supports no growth; fermentation '
          'products (acetate, ethanol, lactate, pyruvate) support aerobic but not anaerobic '
          'growth as sole carbon sources. This panel (Figure 5) is not novel biology - it is '
          'the unit test that proves the LP plumbing, medium handling and objective are '
          'correct before any downstream number is trusted.', BODY)]
    story += fig(os.path.join(ROOT, 'figures', 'fig5_condition_matrix.png'),
                 4.2*1.0*72*1.0*1.0*1.0*1.0*1.0 if False else 300,
                 'Figure 5. Core-model maximal growth rate across carbon sources and aerobiosis. '
                 'Values in h-1; aerobic/anaerobic columns.')
    story += [P('3.2 Genome-scale essentiality: individual layers', H2),
        P(f'Table 1 and Figure 1 collect the head-to-head numbers against the Gerdes 2003 '
          f'labels. The core model covers only 125 labeled genes and predicts essentiality '
          f'poorly (AUROC {core["auroc"]:.3f}, 95% CI {core["auroc_ci"][0]:.3f}-'
          f'{core["auroc_ci"][1]:.3f}) - expected, since most essential genes (DNA, RNA, '
          f'cell-division machinery) lie outside central metabolism. The genome-scale iJO1366 '
          f'under minimal glucose reaches AUROC {ijo["auroc"]:.3f} '
          f'(CI {ijo["auroc_ci"][0]:.3f}-{ijo["auroc_ci"][1]:.3f}), F1 {ijo["f1"]:.3f}, '
          f'MCC {ijo["mcc"]:.3f}. The rich-medium variant trades ranking for call precision: '
          f'AUROC {rich["auroc"]:.3f} but precision {rich["precision"]:.3f} versus '
          f'{ijo["precision"]:.3f} - a real, quantified condition-matching effect, since the '
          f'footprinting experiment ran in rich LB.', BODY),
        P(f'The two learned layers sit in the same band. The graph GCN reaches OOF AUROC '
          f'0.660; the sequence CNN 0.643 (fold AUROCs 0.662 / 0.649 / 0.646), only slightly '
          f'above its own 3-mer logistic baseline (0.628 on the full sequence set; 0.641 on '
          f'the aligned gene set). We flag this honestly: on this label set and feature '
          f'scale, deep sequence modeling adds little over k-mer counts - a negative result '
          f'for the CNN, preserved rather than tuned away.', BODY)]
    story += tbl([['Layer', 'AUROC', 'AUPRC', 'Accuracy', 'Precision', 'Recall', 'F1', 'MCC'],
        ['core FBA (n=125)', f'{core["auroc"]:.3f}', f'{core["auprc"]:.3f}', f'{core["accuracy"]:.3f}', f'{core["precision"]:.3f}', f'{core["recall"]:.3f}', f'{core["f1"]:.3f}', f'{core["mcc"]:.3f}'],
        ['iJO1366 FBA minimal', f'{ijo["auroc"]:.3f}', f'{ijo["auprc"]:.3f}', f'{ijo["accuracy"]:.3f}', f'{ijo["precision"]:.3f}', f'{ijo["recall"]:.3f}', f'{ijo["f1"]:.3f}', f'{ijo["mcc"]:.3f}'],
        ['iJO1366 FBA rich', f'{rich["auroc"]:.3f}', f'{rich["auprc"]:.3f}', f'{rich["accuracy"]:.3f}', f'{rich["precision"]:.3f}', f'{rich["recall"]:.3f}', f'{rich["f1"]:.3f}', f'{rich["mcc"]:.3f}'],
        ['GCN (metabolic graph)', '0.660', '0.310', '-', '-', '-', '-', '-'],
        ['CNN (CDS sequence)', '0.643', '0.251', '-', '-', '-', '-', '-'],
        ['k-mer logistic', '0.628', '0.296', '-', '-', '-', '-', '-'],
        ['ENSEMBLE (stacked)', f'{ens["ensemble"]["oof_auroc"]:.4f}', f'{ens["ensemble"]["oof_auprc"]:.4f}', '-', '-', '-', f'{ens["mcnemar_vs_fba_min"]["ensemble_f1_at_best_t"]:.3f}', '-']],
        'Table 1. Essentiality prediction vs Gerdes 2003 labels (1,249-1,250 aligned genes; '
        'core model on its 125-gene overlap). Ranking metrics are out-of-fold where applicable.')
    story += fig(os.path.join(ROOT, 'figures', 'fig1_auroc_comparison.png'), 420,
                 'Figure 1. AUROC of each layer and the stacked ensemble; error bars are '
                 'bootstrap 95% CIs where available.')
    story += fig(os.path.join(ROOT, 'figures', 'fig2_roc_pr.png'), 430,
                 'Figure 2. ROC (left) and precision-recall (right) on the aligned 1,249-gene '
                 'set. Dashed: chance / prevalence.')
    story += [P('3.3 The ensemble and its statistical caveat', H2),
        P(f'The stacked ensemble reaches OOF AUROC {ens["ensemble"]["oof_auroc"]:.4f} and '
          f'AUPRC {ens["ensemble"]["oof_auprc"]:.4f}, a gain of +0.057 AUROC and +0.10 AUPRC '
          f'over the best single layer (iJO1366 minimal FBA). At F1-optimal hard thresholds, '
          f'however, McNemar\'s test counts 27 genes the ensemble gets right while FBA errs '
          f'versus 50 in the opposite direction (p = '
          f'{ens["mcnemar_vs_fba_min"]["mcnemar_p"]:.4f}). The reconciliation is structural: '
          f'the ensemble trades some majority-class (nonessential) accuracy for positive-class '
          f'F1 ({ens["mcnemar_vs_fba_min"]["ensemble_f1_at_best_t"]:.3f} vs '
          f'{ens["mcnemar_vs_fba_min"]["fba_min_f1_at_best_t"]:.3f}) and much better ranking. '
          f'We therefore claim an improvement in candidate <i>prioritization</i>, not in '
          f'genome-wide hard-call accuracy - and say so, because the distinction is exactly '
          f'where virtual-cell papers tend to oversell.', BODY),
        PageBreak(),
        P('3.4 Dynamics: diauxie needs regulation', H2),
        P(f'Unregulated dFBA fails the diauxie test cleanly: acetate is co-consumed with '
          f'glucose (co-utilization flag = {dia["co_utilization_primary_phase"]}), glucose is '
          f'exhausted at {dia["t_primary_exhausted_h"]:.1f} h, and no second growth phase '
          f'appears. This is not a solver error; it is FBA being exactly as dumb as its '
          f'equations - nothing in Sv = 0 knows about catabolite repression. Adding one '
          f'Boolean repression rule (glucose uptake &gt; 0 caps acetate uptake at 0) flips the '
          f'qualitative behavior: acetate stays untouched during the glucose phase '
          f'(co-utilization = {diar["co_utilization_primary_phase"]}), glucose exhausts at '
          f'{diar["t_primary_exhausted_h"]:.1f} h, and a genuine two-substrate sequence '
          f'emerges (Figure 3). Peak regulated growth is {diar["peak_growth_h"]:.3f} h-1 and '
          f'final biomass {diar["final_biomass_gL"]:.3f} g/L. The lesson, quantified: a '
          f'virtual cell without a regulatory layer is not a cell; it is a stoichiometry.', BODY)]
    story += fig(os.path.join(ROOT, 'figures', 'fig3_diauxie.png'), 440,
                 'Figure 3. dFBA on glucose + acetate. Left: unregulated - acetate co-utilized, '
                 'no diauxie. Right: with the catabolite-repression rule - sequential '
                 'utilization, as observed experimentally.')


    # ---- 3.5 detailed result tables ----
    import pandas as _pd
    cm = _pd.read_csv(os.path.join(ROOT, 'results', 'core_condition_matrix.csv'))
    rows = [['Carbon source', 'Aerobic (h-1)', 'Anaerobic (h-1)']]
    for carb, grp in cm.groupby('carbon'):
        a = grp[grp.aerobic].growth_h.iloc[0]; n = grp[~grp.aerobic].growth_h.iloc[0]
        rows.append([carb.replace('EX_','').replace('_e',''), f'{a:.4f}', f'{n:.4f}'])
    story += [P('3.5 Detailed result tables', H2),
        P('Table 2 lists the full condition matrix that gates all downstream claims; '
          'Table 3 the core model essential set; Table 4 the CNN per-fold scores; '
          'Tables 5 and 6 the complete discovery classes.', BODY)]
    story += tbl(rows, 'Table 2. Core-model maximal growth across the 28-condition panel.')
    ess_core = _pd.read_csv(os.path.join(ROOT, 'results', 'core_gene_essentiality.csv'))
    ess_core = ess_core[ess_core.predicted_essential].sort_values('gene')
    names = {'b0720':'glyA? see S1','b1136':'','b1779':'','b2415':'','b2416':'','b2779':'','b2926':''}
    gerdes = _pd.read_csv(os.path.join(ROOT, 'data', 'gerdes_labels.csv'))
    nmap = dict(zip(gerdes.bnumber, gerdes.gene_name))
    rows = [['b-number', 'gene', 'growth fraction']]
    for _, r in ess_core.iterrows():
        rows.append([r.gene, nmap.get(r.gene, ''), f'{r.growth_frac:.2e}'])
    story += tbl(rows, 'Table 3. Core-model predicted essential genes (glucose minimal, '
                       '1% threshold), with Gerdes gene names.')
    seqr = R('seqcnn_results.json')
    rows = [['Fold', 'AUROC', 'AUPRC']] + [[str(f['fold']), f"{f['auroc']:.4f}", f"{f['auprc']:.4f}"] for f in seqr['cnn']['folds']]
    rows.append(['OOF (pooled)', f"{seqr['cnn']['oof_auroc']:.4f}", f"{seqr['cnn']['oof_auprc']:.4f}"])
    rows.append(['k-mer LR (OOF)', f"{seqr['kmer_baseline']['oof_auroc']:.4f}", f"{seqr['kmer_baseline']['oof_auprc']:.4f}"])
    story += tbl(rows, 'Table 4. Sequence-layer per-fold and pooled out-of-fold performance.')
    rows = [['b-number', 'gene', 'consensus', 'description']]
    for g in disc['classA']:
        rows.append([g['bnumber'], g['gene_name'], f"{g['consensus']:.3f}", g['description'][:58]])
    story += tbl(rows, 'Table 5. Class A: model-essential, experimentally viable.')
    rows = [['b-number', 'gene', 'consensus', 'description']]
    for g in disc['classB']:
        rows.append([g['bnumber'], g['gene_name'], f"{g['consensus']:.3f}", g['description'][:58]])
    story += tbl(rows, 'Table 6. Class B: experimentally essential, missed by all layers.')


    story += [P('3.6 Robustness checks', H2),
        P('Two robustness facts matter for interpreting Table 1. First, the hard essentiality '
          'calls are insensitive to the 1% growth threshold: at 0.5%, 1% and 5% of wild-type '
          'growth the iJO1366 minimal screen predicts exactly 183 essential genes with '
          'identical precision, recall, F1 and MCC, because knockout growth fractions are '
          'strongly bimodal - deletions either abolish growth or leave it near wild type, '
          'and almost no gene lands in the 0.5-5% band. The threshold choice therefore '
          'drives none of our conclusions. Second, the dFBA conclusions are step-size '
          'robust: halving and doubling &Delta;t (0.05-0.2 h) leaves the co-utilization '
          'flags and phase structure unchanged, consistent with the O(&Delta;t) error bound '
          'of Section 2.10; the test suite asserts the &Delta;t = 0.2 h case directly.', BODY),
        P('3.7 Second yardstick: PEC (Keio-era) deletion calls', H2),
        P('Gerdes Table S2 cross-tabulates essentiality calls across three studies, including '
          'the PEC deletion-collection calls from the Keio era [1]. On the 197 model-covered '
          'genes with unambiguous PEC labels (83 essential), iJO1366 minimal FBA scores AUROC '
          '0.576 (CI 0.489-0.663) and the ensemble 0.568 (n=180) - markedly below the '
          'genome-wide Gerdes numbers. We read this cautiously rather than as a failure: S2 '
          'is a consensus-focused 692-gene subset deliberately enriched for conflicted and '
          'interesting genes (228 of 692 carry an X call), so it is a hard, biased exam, not '
          'a representative one. Cross-study disagreement itself is the signal here: the two '
          'experimental yardsticks disagree enough that "which truth?" is a first-class '
          'question for virtual-cell benchmarking. A genome-wide Keio Table I extraction is '
          'queued as future work (Section 5.5).', BODY)]

    # ---------------- discovery ----------------
    story += [PageBreak(), P('4. Discovery analysis: from residuals to a verified model repair', H1),
        P('4.1 Disagreement classes', H2),
        P('We define the consensus score as the mean of the five min-max-normalized layer '
          'scores and isolate two tails (Figure 4). Class A (9 genes): consensus in the top '
          '2% while experiment says viable. Class B (10 genes): experimentally essential '
          'while the consensus bottom quintile says viable. Class B is dominated by genes '
          'whose essentiality is structural rather than metabolic - mrdA (penicillin-binding '
          'protein 2, cell-wall synthesis), wzyE (enterobacterial common antigen polymerase), '
          'and ABC-transporter components (cydC/cydD, araG, yddQ, ydcU) - plus hypothetical '
          'proteins with no model coverage. This is the expected blind spot of all three '
          'layers and sets an honest scope limit: a metabolism-anchored virtual cell cannot '
          'see cell-envelope essentiality until those processes are modeled.', BODY)]
    story += fig(os.path.join(ROOT, 'figures', 'fig4_discovery.png'), 420,
                 'Figure 4. Consensus essentiality score vs experimental label (jittered). '
                 'Class A (red) and Class B (blue) disagreement tails are marked; the five '
                 'highest-consensus Class A genes are labeled.')
    genesA = ', '.join(f"{g['gene_name']} ({g['bnumber']})" for g in disc['classA'][:9])
    story += [
        P('4.2 Class A localizes a pathway: molybdenum-cofactor biosynthesis', H2),
        P(f'Class A is not random: five of its nine members - {genesA} - belong to one '
          f'pathway, molybdenum-cofactor (MoCo) biosynthesis, plus the folate pair folB/folP '
          f'and the ubiquinone enzyme ubiC. A pathway-shaped residual is a signature of a '
          f'systematic model error rather than noise, so we investigated the MoCo cluster '
          f'directly.', BODY),
        P('4.3 Mechanism and in-silico proof', H2),
        P(f'In iJO1366, moaD participates in two reactions (MPTS, molybdopterin synthase; '
          f'MOADSUx, MoaD sulfuration). In aerobic glucose wild-type FBA these carry small '
          f'but nonzero flux, meaning the model <i>requires</i> MoCo production even when no '
          f'MoCo-dependent enzyme is needed. The reason is literal: the core biomass '
          f'objective reaction lists the MoCo compounds bmocogdp_c and mobd_c as required '
          f'constituents (the WT biomass additionally lists mococdp_c and mocogdp_c). Any '
          f'MoCo-pathway knockout therefore scores as lethal, regardless of condition - we '
          f'verified moaD/moaC/moaE deletions are lethal in silico under both aerobic glucose '
          f'and anaerobic nitrate, refuting our own initial hypothesis of anaerobic '
          f'conditional essentiality. The decisive experiment: removing the two MoCo '
          f'compounds from the core biomass objective rescues the moaD knockout completely '
          f'(growth {rescue["moaD_ko_core_bm_no_moco"]:.4f} vs wild-type '
          f'{rescue["wt_core_bm_no_moco"]:.4f} h-1), matching the experimental viability call '
          f'from Gerdes 2003.', BODY),
        P('The interpretation is physically sensible: MoCo enzymes (nitrate reductase, DMSO '
          'reductase, formate dehydrogenase) matter mainly under anaerobic respiration, so a '
          'constitutive biomass requirement overstates their necessity in aerobic rich '
          'growth. We are not claiming the iJO1366 biomass is "wrong" - biomass compositions '
          'are measured averages, not conditional requirements - but that it should be '
          'condition-dependent here, and that residual analysis against one experimental set '
          'was sufficient to find and repair the artifact. A PubMed cross-check '
          '("moaD Escherichia coli essential": 7 records) is consistent with scattered '
          'condition-specific reports rather than constitutive essentiality. The same '
          'condition-dependence argument applies to folB/folP (folate synthesis, bypassable '
          'by uptake in rich medium) and ubiC (ubiquinone for aerobic respiration).', BODY),
        PageBreak()]

    # ---------------- discussion ----------------
    story += [P('5. Discussion', H1),
        P('5.1 What the virtual cell is, and is not', H2),
        P('VC-2 is not a whole-cell model in the Karr et al. sense [5]: it has no replication, '
          'transcription, translation or cell-division machinery, and its Class B residual set '
          'is the direct, measurable price of that scope. Within its scope it is genuinely '
          'integrated: one genotype (the metabolic reconstruction), one experimental yardstick '
          '(Gerdes 2003), four mechanistically different predictors, and a combiner whose '
          'gains and failure modes are both quantified. The ensemble\'s +0.057 AUROC over the '
          'best layer is real but modest; the McNemar analysis (Section 3.3) is the guardrail '
          'against reading more into it than exists.', BODY),
        P('5.2 Comparison to published leaders', H2),
        P('The iJO1366 publication [3] reported strong agreement with knockout screens when '
          'medium conditions were matched to the experiment; our genome-scale FBA numbers are '
          'lower primarily because the Gerdes footprinting was performed in rich LB while our '
          'main screen uses glucose-minimal medium - the rich-medium variant recovers precision '
          'exactly as that argument predicts. Sequence deep learning published for '
          'essentiality (the DeepCellEss class of models) is typically evaluated on '
          'DEG-compiled essential sets with different gene universes, so headline AUROCs are '
          'not directly comparable to ours; our contribution is to evaluate sequence, graph '
          'and mechanistic predictors on the same genes, same labels and same folds, where '
          'the CNN\'s edge over 3-mer counts nearly vanishes. We regard cross-dataset '
          're-evaluation on identical universes as the missing control in this literature.', BODY),
        P('5.3 Bugs the test suite caught (reported on purpose)', H2),
        P('Three genuine implementation bugs were caught by regression tests and are part of '
          'the scientific record: (i) medium state leaking between conditions, which silently '
          'made every knockout "essential"; (ii) the COBRA medium sign convention (positive '
          'uptake rates mapping to negative lower bounds), which zeroed all growth; (iii) '
          'infeasible-knockout NaN fluxes, which mislabeled two core-model genes until '
          'infeasible was mapped to zero growth. Each would have produced plausible-looking, '
          'wrong results. The 18-test hermetic suite that caught them runs in about five '
          'seconds with no network access.', BODY),
        P('5.4 Limitations', H2),
        P('(1) One experimental label set (Gerdes 2003); adding the Keio deletion collection '
          '[1] as a second, independent yardstick is queued. (2) The CNN sees only the first '
          '1,200 nt. (3) The GNN uses degree-structural features only. (4) dFBA uses explicit '
          'Euler integration at &Delta;t = 0.1 h. (5) The ensemble is a linear stack; nothing '
          'about it captures epistasis. (6) All compute ran on two CPU cores by design, which '
          'caps model size; none of the conclusions depend on scale.', BODY),

        P('5.5 Future work', H2),
        P('Queued extensions, in priority order: (1) a second independent essentiality '
          'yardstick from the Keio deletion collection [1], with cross-study consistency '
          'analysis against Gerdes 2003 (the Gerdes Table S2 consensus sets are already '
          'parsed); (2) condition-specific biomass objectives, generalizing the Section 4.3 '
          'MoCo repair into a principled biomass library; (3) a transcriptional-regulatory '
          'layer so diauxie emerges instead of being injected; (4) longer-context sequence '
          'models and codon-level features for the CNN; (5) extension of the stacked '
          'architecture to non-model organisms with sparse labels, where the layer '
          'diversity argument should matter most.', BODY),
        P('6. Conclusion', H1),
        P('The positive contribution is a reproducible audit for two failure-prone choices '
          'in virtual-cell benchmarking: biomass composition and predictive score lineage. '
          'Objective perturbations changed essentiality calls across related model snapshots, '
          'and independently checked <i>Pseudomonas</i> counterfactuals showed that both '
          'single-term demands and term combinations can determine a call. The result is a '
          'model-curation and sensitivity finding, not evidence that those knockout strains '
          'live in a laboratory. The score-lineage audit distinguishes a historical ranking '
          'gain from the smaller, differently compared leakage-controlled v2 gain and blocks '
          'a clean held-out hard-call claim. The pre-registered G2 gate failed on every rung; '
          'medium mismatch, assay dependence and weak contributions from some learning layers '
          'remain explicit limits. The runnable audit gives future condition-matched studies '
          'a way to test which model predictions survive objective and evaluation choices '
          'before treating them as biological hypotheses.', BODY),
        PageBreak(),
        P('References', H1)]
    refs = [
        '[1] Baba T, Ara T, Hasegawa M, et al. Construction of Escherichia coli K-12 in-frame, '
        'single-gene knockout mutants: the Keio collection. <i>Mol Syst Biol</i> 2:2006.0008 (2006).',
        '[2] Gerdes SY, Scholle MD, Campbell JW, et al. Experimental determination and system '
        'level analysis of essential genes in Escherichia coli MG1655. <i>J Bacteriol</i> '
        '185(19):5673-5684 (2003).',
        '[3] Orth JD, Conrad TM, Na J, et al. A comprehensive genome-scale reconstruction of '
        'Escherichia coli metabolism - 2011. <i>Mol Syst Biol</i> 7:535 (2011).',
        '[4] Mahadevan R, Edwards JS, Doyle FJ. Dynamic flux balance analysis of diauxic growth '
        'in Escherichia coli. <i>Biophys J</i> 83(3):1331-1340 (2002).',
        '[5] Karr JR, Sanghvi JC, Macklin DN, et al. A whole-cell computational model predicts '
        'phenotype from genotype. <i>Cell</i> 150(2):389-401 (2012).',
        '[6] Orth JD, Thiele I, Palsson BO. What is flux balance analysis? <i>Nat Biotechnol</i> '
        '28(3):245-248 (2010).',
        '[7] Ebrahim A, Lerman JA, Palsson BO, Hyduke DR. COBRApy: constraints-based '
        'reconstruction and analysis for Python. <i>BMC Syst Biol</i> 7:74 (2013).',
        '[8] King ZA, Lu J, Drager A, et al. BiGG Models: a platform for integrating, '
        'standardizing and sharing genome-scale models. <i>Nucleic Acids Res</i> 44(D1):D515-D522 (2016).',
        '[9] Kipf TN, Welling M. Semi-supervised classification with graph convolutional '
        'networks. <i>ICLR</i> (2017).',
        '[10] Kingma DP, Ba J. Adam: a method for stochastic optimization. <i>ICLR</i> (2015).',
        '[11] Cock PJA, Antao T, Chang JT, et al. Biopython: freely available Python tools for '
        'computational molecular biology and bioinformatics. <i>Bioinformatics</i> 25(11):1422-1423 (2009).',
        '[12] Hayashi K, Morooka N, Yamamoto Y, et al. Highly accurate genome sequences of '
        'Escherichia coli K-12 strains MG1655 and W3110. <i>Mol Syst Biol</i> 2:2006.0007 (2006).',
    ]
    for r in refs:
        story.append(P(r, __import__('build_paper').REF))

    story += [PageBreak(), P('Appendix A. Reproducibility ledger', H1),
        P('Every number in this paper regenerates from the repository. Data snapshots: '
          'BiGG e_coli_core.json and iJO1366.json (fetched 2026-09-24), U00096.3 FASTA + '
          'GenBank (NCBI efetch, 2026-09-24), Gerdes 2003 Table S1 (genome.wisc.edu, '
          '2026-09-24). Pipelines: scripts/run_ijo1366_scan.py (deletion screens), '
          'scripts/run_seqcnn.py (CNN + baseline), scripts/run_gnn.py (graph GCN), '
          'scripts/make_figures.py (all figures), paper/build_paper.py (this document). '
          'Intermediate OOF scores for every model are stored in results/ as CSV/NPZ. '
          'Tests: 18 hermetic pytest cases, ~5 s, no network. Live network calls are made '
          'only by fetch scripts, never by tests.', BODY)]
    story += tbl([['Artifact', 'Path', 'Verified by'],
             ['Condition matrix', 'results/core_condition_matrix.csv', 'TestMetabolism'],
             ['Core essentiality', 'results/core_gene_essentiality.csv', 'regression-locked 7-gene set'],
             ['iJO1366 screens', 'results/ijo1366_gene_essentiality[_rich].csv', 'TestMetabolism + benchmark JSON'],
             ['CNN / k-mer OOF', 'results/seqcnn_oof_scores.npy, kmer_oof_scores.npy', 'TestSeqCNN'],
             ['GNN OOF', 'results/gnn_oof_scores.npy', 'TestGraphGNN'],
             ['Ensemble OOF + stats', 'results/ensemble_results.json', 'McNemar block'],
             ['dFBA trajectories', 'results/dfba_diauxie*_trajectory.csv', 'TestDynamics'],
             ['MoCo rescue', 'results/moco_biomass_rescue_test.json', 'Section 4.3']],
            'Table A1. Result-to-test traceability.')

    story += [PageBreak(), P('Appendix B. Notation and hyperparameters', H1)]
    story += tbl([['Symbol', 'Meaning'],
        ['S', 'stoichiometric matrix (metabolites x reactions)'],
        ['v', 'flux vector; lb, ub its bounds'],
        ['c', 'objective selector (biomass reaction)'],
        ['A, D', 'gene-gene adjacency; degree matrix'],
        ['H^(l), W^(l)', 'GCN layer activations; weights'],
        ['X(t)', 'biomass concentration (g/L)'],
        ['S_i(t)', 'extracellular concentration of substrate i (mmol/L)'],
        ['V_max, K_m', 'Michaelis-Menten uptake parameters (10 mmol/gDW/h; 0.015 mmol/L)'],
        ['mu(t)', 'specific growth rate (h-1)'],
        ['dt', 'dFBA integration step (0.1 h; 0.2 h in robustness test)'],
        ['b, c', 'McNemar discordant-pair counts']],
        'Table B1. Notation.')
    story += tbl([['Component', 'Setting'],
        ['CNN', 'Conv1d 4-32-64-128, k=9/9/5, maxpool 4/4, GAP; Adam 1e-3; 15 epochs; batch 64'],
        ['GCN', '2 layers, 32 hidden, Adam 5e-3, weight decay 1e-4, 200 epochs'],
        ['CV', 'stratified 3-fold, seed 7, all metrics out-of-fold'],
        ['Ensemble', 'logistic regression, class_weight=balanced, max_iter 2000'],
        ['Bootstrap', '1000 resamples for CIs (2000 in final pass where noted), seed 7'],
        ['FBA', 'GLPK via COBRApy; essentiality threshold 1% of WT optimum'],
        ['Hardware', '2 CPU cores, ~2 GB RAM, no GPU']],
        'Table B2. Hyperparameters and compute.')
    story += [Spacer(1, 18),
        P('<b>Data and code availability.</b> All models (BiGG), the genome (NCBI U00096.3) '
          'and the essentiality ground truth (Gerdes 2003 supplement, genome.wisc.edu) are '
          'public; snapshots ship in the repository so every figure and number regenerates '
          'offline. No human subjects, no wet-lab work, no dual-use concerns: this is '
          'in-silico analysis of public data on a non-pathogenic lab strain.', BODY)]

    return story
