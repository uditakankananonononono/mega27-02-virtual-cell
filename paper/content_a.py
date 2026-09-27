"""Paper content, part A: title, abstract, introduction, methods."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_paper import P, H1, H2, BODY, EQ, TITLE, PageBreak, Spacer

def story_a(story, R):
    # ---------------- title page ----------------
    story += [Spacer(1, 120),
        P('VC-2: Auditing Biomass-Objective Sensitivity and Score Leakage in an '
          '<i>Escherichia coli</i> Virtual-Cell Benchmark', TITLE),
        Spacer(1, 24),
        P('MEGA-PROGRAM-27, Item 2 - Virtual Cell', H1),
        P('Author: Udita Phookan', BODY),
        P('24 September 2026', BODY),
        Spacer(1, 60),
        P('<b>Repository:</b> mega27-02-virtual-cell &nbsp;&nbsp;|&nbsp;&nbsp; '
          '<b>Code:</b> Python 3.10, COBRApy, PyTorch 2.14 (CPU), scikit-learn &nbsp;&nbsp;|&nbsp;&nbsp; '
          '<b>Tests:</b> 27/27 tests passing at the current working release', BODY),
        PageBreak()]

    # ---------------- abstract ----------------
    story += [P('Abstract', H1),
        P('A virtual-cell benchmark is only as reliable as its biomass objective and its evaluation '
          'lineage. We built an executable <i>Escherichia coli</i> model with flux balance analysis, '
          'sequence and graph learners, a stacked score, and dynamic simulation, then asked two '
          'audit questions: which essentiality calls change when individual or combined biomass '
          'demands are perturbed, and which predictive comparisons survive scrutiny of training '
          'folds, thresholds, and experimental conditions? The contribution is an auditable '
          'model-sensitivity and leakage-control workflow, not a new biological viability claim. '
          'In 68 scorable MoCo-positive BiGG model snapshots, 63 changed at least one predicted '
          'gene-essentiality call after the named coefficients were removed. Those snapshots '
          'are not 68 independent organisms or phenotype replications. A separately coded '
          '<i>Pseudomonas</i> iJN1463 counterfactual rescued seven originally zero-growth '
          'knockouts after its three MoCo coefficients were removed together; its 102-term '
          'single-removal audit found 23 terms rescuing at least one call, while no individual '
          'MoCo term did. Pairwise controls further limit any claim of MoCo uniqueness. These '
          'results expose combinatorial objective sensitivity inside reconstructions, not '
          'experimental survival.', BODY),
        P('The original Gerdes 2003 hard-call G2 gate failed on all five locked rungs. A later '
          'score-lineage audit found that threshold selection reused out-of-fold scores whose '
          'training complements could contain evaluation labels; the hard-call counts are '
          'descriptive, not clean held-out evidence. Gerdes transposon labels also come from '
          'a different medium than the glucose-minimal simulation. The historical stacked '
          'ranking gain of +0.057 AUROC over FBA is not a leakage-controlled gain. In the '
          'separate v2 analysis, the leakage-controlled feature set gained +0.024 AUROC over '
          'the reproduced, rescaled v1 score, a different comparator. We retain these '
          'comparisons with their distinct denominators and limitations rather than claim '
          'a benchmark-beating hard-call predictor. The model and audits run on two CPU '
          'cores; condition-matched phenotype, fully nested refits and independent '
          'validation remain open.', BODY),
        P('<b>Keywords:</b> virtual cell, flux balance analysis, gene essentiality, graph '
          'convolutional network, convolutional neural network, dynamic FBA, model ensemble, '
          '<i>Escherichia coli</i>, biomass objective function', BODY),
        PageBreak()]

    # ---------------- introduction ----------------
    story += [P('1. Introduction', H1),
        P('A constraint-based virtual cell turns reaction stoichiometry and a biomass objective '
          'into growth predictions. Those predictions may shift when the objective demands '
          'different metabolites, even though no measured phenotype changes. This project '
          'tests that dependence explicitly, with per-accession model snapshots, knockout '
          'counterfactuals, individual-term and pairwise perturbations, and independent '
          'linear-program checks. The mechanistic output is a map of which calls are forced '
          'by objective choices and which are robust to those changes.', BODY),
        P('The second problem is evaluation integrity. We implemented metabolic, sequence, '
          'graph and dynamic layers and compared essentiality scores with published Gerdes '
          '2003 labels and a distinct CRISPRi assay. Different media and assay mechanisms '
          'make those yardsticks useful but not interchangeable. Rescaling a stack can improve '
          'performance without adding information; upstream label reuse can make threshold '
          'selection appear independent when it is not. The score-lineage and comparator '
          'audits therefore accompany every ranking or hard-call result. The failed locked '
          'hard-call gate remains failed. Sections 2-4 describe the model, Section 7 the '
          'cross-reconstruction objective sensitivity and calibration audit, and the '
          'appendices retain full positive and negative checks.', BODY),
        P('1.1 Why essentiality is the right first benchmark', H2),
        P('A gene-deletion call gives each model disagreement a traceable address: a reaction '
          'rule, a biomass ingredient, a medium assumption or a score threshold. It does not '
          'give a condition-independent truth. Gerdes transposon labels were measured in a '
          'rich-medium setting, whereas the principal genome-scale FBA screen used a '
          'glucose-minimal simulation. We use the labels to locate discrepancies and compare '
          'ranking methods, but do not treat their disagreement as proof that a particular '
          'model biomass term is erroneous. Objective perturbations identify model-dependent '
          'calls; phenotype experiments would have to decide which predictions hold in vivo.', BODY),
        P('2. Methods', H1),
        P('2.1 Models and experimental ground truth', H2),
        P('Two published constraint-based reconstructions were used exactly as distributed by '
          'the BiGG Models knowledgebase [8]: e_coli_core (95 reactions, 137 genes; the '
          'central-metabolism teaching model) and iJO1366 (2,583 reactions, 1,805 metabolites, '
          '1,367 genes; the 2011 genome-scale reconstruction of E. coli K-12 MG1655 [3]). '
          'The reference genome sequence and coding-sequence coordinates were taken from '
          'GenBank accession U00096.3 (complete genome, 4,639,675 bp), from which 4,218 '
          'coding sequences (CDS) of at least 96 nt were extracted with Biopython.', BODY),
        P('Experimental essentiality labels come from Gerdes et al. 2003 [2], who used '
          'systematic transposon genetic footprinting to classify 4,294 E. coli MG1655 '
          'protein-coding genes. Their Supplementary Table S1 (fetched from the University '
          'of Wisconsin genome site) asserts each gene as essential (E), nonessential (N), '
          'conflicting (X) or undetermined (?). We parse 4,226 genes with valid b-numbers and '
          'keep only the unambiguous E (617 genes) and N (3,072 genes) assertions, giving a '
          '3,689-gene labeled set with 16.7% prevalence. Genes asserted X or ? are excluded '
          'from every metric. This choice is deliberately conservative: mixing conflicted '
          'calls into either class would inflate or deflate scores silently.', BODY),
        P('2.2 Constraint-based layer: flux balance analysis', H2),
        P('FBA represents metabolism as a stoichiometric matrix S (metabolites x reactions) '
          'and solves the linear program', BODY),
        P('maximize &nbsp; c<sup>T</sup>v &nbsp;&nbsp; subject to &nbsp;&nbsp; Sv = 0, '
          '&nbsp; lb &le; v &le; ub,', EQ),
        P('where v is the flux vector and c selects the biomass objective reaction. Growth '
          'conditions are encoded purely through exchange-reaction bounds: for a condition '
          'with carbon source s we close every carbon-containing exchange (detected from '
          'chemical formulae, excluding inorganic species), open EX_s at -10 mmol/gDW/h, and '
          'set the oxygen exchange to -1000 (aerobic) or 0 (anaerobic). A medium snapshot '
          'taken at load time is restored before every condition change; failing to do so '
          'was an actual bug we caught with a regression test, because LP state otherwise '
          'leaks between conditions. All LP solves use GLPK through COBRApy [7].', BODY),
        P('2.3 In-silico gene essentiality screen', H2),
        P('Single-gene deletions were computed with COBRApy\'s single_gene_deletion, which '
          'applies each gene\'s gene-protein-reaction (GPR) rules and zeroes the affected '
          'reactions. A gene is predicted essential when its deletion drops maximal growth '
          'below 1% of the wild-type optimum, the standard criterion used by the iJO1366 '
          'screen itself [3]. Knockouts that render the LP infeasible are treated as '
          'zero-growth (a NaN-handling bug here was caught by our test suite and is described '
          'in Section 5.3). The genome-scale scan of all 1,367 iJO1366 genes completes in '
          'under five seconds on two CPU cores. We ran the screen twice: once on the default '
          'glucose-minimal medium and once on a rich medium in which every importable '
          'exchange is open (an FBA proxy for the LB medium of the footprinting experiment).', BODY),
        P('2.4 Sequence layer: 1D-CNN on coding sequence', H2),
        P('Each CDS is one-hot encoded (4 channels, first 1,200 nt, zero-padded) and passed '
          'through a compact DeepCellEss-style stack: Conv1d(4-32, k=9) - ReLU - MaxPool(4) - '
          'Conv1d(32-64, k=9) - ReLU - MaxPool(4) - Conv1d(64-128, k=5) - ReLU - global average '
          'pool - linear head, trained with BCEWithLogitsLoss and a positive-class weight of '
          'N<sub>neg</sub>/N<sub>pos</sub> to handle the 1:5 class imbalance. Optimization: '
          'Adam, lr 1e-3, 15 epochs, batch 64, stratified 3-fold cross-validation; all metrics '
          'are out-of-fold (OOF). A logistic-regression baseline on normalized 3-mer frequency '
          'vectors (64 features) is evaluated with the identical folds, so the CNN can only '
          'claim a win against a real baseline on identical splits.', BODY),
        P('2.5 Graph layer: GCN over the metabolic gene graph', H2),
        P('From the iJO1366 stoichiometry we build a gene-gene graph: genes are nodes, and '
          'two genes share an edge when reactions they encode (via GPR rules) share at least '
          'one metabolite, excluding 16 currency metabolites (H+, H2O, ATP/ADP/AMP, '
          'NAD(H), NADP(H), CoA, phosphate, CO2, O2, NH4+) and hyper-connected hubs '
          '(metabolites touching more than 25 reactions) to avoid clique explosions. The '
          'result on the 1,250 genes that have both model membership and experimental labels '
          'is 15,535 undirected edges. Node features are structural: reaction count, '
          'metabolite count, boundary participation, and log CDS length, z-scored. The model '
          'is a two-layer graph convolutional network [9] with the symmetric normalization', BODY),
        P('H<sup>(l+1)</sup> = &sigma;( D<sup>-1/2</sup> (A + I) D<sup>-1/2</sup> '
          'H<sup>(l)</sup> W<sup>(l)</sup> ),', EQ),
        P('implemented in pure PyTorch (dense adjacency, 32 hidden units, 200 epochs, '
          'Adam 5e-3, weight decay 1e-4) so the test suite needs no compiled graph library. '
          'The same 3-fold stratified protocol is used.', BODY),
        P('2.6 Dynamics layer: regulated dynamic FBA', H2),
        P('Dynamic FBA follows the static-optimization approach (SOA) of Mahadevan et al. '
          '2002 [4]. At each step of size &Delta;t, extracellular substrate concentrations '
          'S<sub>i</sub> set uptake bounds through a Michaelis-Menten rule', BODY),
        P('lb<sub>i</sub>(t) = -V<sub>max</sub> S<sub>i</sub>(t) / (K<sub>m</sub> + S<sub>i</sub>(t)),', EQ),
        P('FBA returns the growth rate &mu;(t) and exchange fluxes v<sub>i</sub>(t), and the '
          'state is integrated explicitly:', BODY),
        P('X(t+&Delta;t) = X(t) e<sup>&mu;(t)&Delta;t</sup>, &nbsp;&nbsp; '
          'S<sub>i</sub>(t+&Delta;t) = S<sub>i</sub>(t) + v<sub>i</sub>(t) X(t) &Delta;t.', EQ),
        P('The test system is the classic diauxie: E. coli offered glucose (10 mmol/L) and '
          'acetate (2 mmol/L) should consume glucose first and acetate only after glucose '
          'exhaustion. Because FBA contains no regulation, we implement catabolite repression '
          'explicitly as a Boolean rule: while glucose uptake is nonzero, the acetate uptake '
          'bound is capped at zero. Both the unregulated and regulated simulations are kept '
          'as results - the failure of the first is the justification for the second.', BODY),
        P('2.7 Ensemble and statistics', H2),
        P('Five signals per gene - FBA minimal score (1 - growth fraction), FBA rich score, '
          'GNN OOF score, CNN OOF score, k-mer baseline OOF score - are stacked with a '
          'logistic-regression combiner, itself evaluated out-of-fold under 3-fold '
          'stratification on the 1,249 genes present in every signal. Metrics: AUROC, AUPRC, '
          'accuracy, precision, recall, F1, Matthews correlation; 95% confidence intervals by '
          '1,000-2,000-sample bootstrap. Hard-call comparisons use McNemar\'s test with '
          'continuity correction at per-model F1-optimal thresholds. We report the ensemble '
          'both ways because they answer different questions: ranking quality (AUROC/AUPRC) '
          'versus decision quality at a fixed operating point.', BODY),
        PageBreak()]
    return story
