"""Paper content, part A: title, abstract, introduction, methods."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_paper import P, H1, H2, BODY, EQ, TITLE, PageBreak, Spacer

def story_a(story, R):
    # ---------------- title page ----------------
    story += [Spacer(1, 120),
        P('VC-2: A Modular Virtual Cell for <i>Escherichia coli</i> K-12 Integrating '
          'Constraint-Based Metabolism, Sequence and Graph Deep Learning, and '
          'Dynamic Simulation, Benchmarked Against Experimental Gene Essentiality', TITLE),
        Spacer(1, 24),
        P('MEGA-PROGRAM-27, Item 2 - Virtual Cell', H1),
        P('Author: Udita Phookan (computational work executed with Instinct)', BODY),
        P('24 September 2026', BODY),
        Spacer(1, 60),
        P('<b>Repository:</b> mega27-02-virtual-cell &nbsp;&nbsp;|&nbsp;&nbsp; '
          '<b>Code:</b> Python 3.10, COBRApy, PyTorch 2.14 (CPU), scikit-learn &nbsp;&nbsp;|&nbsp;&nbsp; '
          '<b>Tests:</b> 23/23 hermetic tests passing', BODY),
        PageBreak()]

    # ---------------- abstract ----------------
    story += [P('Abstract', H1),
        P('A "virtual cell" must do two things at once: reproduce what experiments already '
          'show, and fail in ways that teach us something. We built VC-2, a modular in-silico '
          '<i>Escherichia coli</i> K-12 MG1655 cell with four interacting layers: (i) a '
          'constraint-based metabolic layer running flux balance analysis (FBA) on the published '
          'e_coli_core and iJO1366 reconstructions, including a full 1,367-gene single-deletion '
          'screen under minimal and rich media; (ii) a sequence layer, a one-dimensional '
          'convolutional network trained on 4,218 real coding sequences extracted from the '
          'U00096.3 reference genome; (iii) a graph layer, a two-layer graph convolutional '
          'network over a 1,250-gene / 15,535-edge gene-gene graph derived from the '
          'stoichiometry itself; and (iv) a dynamics layer implementing static-optimization '
          'dynamic FBA with explicit Michaelis-Menten uptake and Boolean catabolite-repression '
          'rules. All layers were scored against the Gerdes et al. 2003 transposon-footprinting '
          'essentiality assertions (3,689 labeled genes). Individually the layers are mediocre '
          'and differently wrong (AUROC 0.63-0.67). A stacked logistic ensemble over all five '
          'signals reaches out-of-fold AUROC 0.7225 and AUPRC 0.414, a +0.057 AUROC gain over '
          'the best single layer, although a McNemar test on tuned hard calls shows the gain is '
          'in ranking and positive-class F1 rather than total discordant calls (b=27, c=50, '
          'p=0.012). Consensus residual analysis then localized a systematic reconstruction '
          'artifact: the entire molybdenum-cofactor (MoCo) biosynthesis cluster '
          '(moaC/moaD/moaE/mobA/moeB) is predicted essential because the iJO1366 biomass '
          'objective lists molybdopterin compounds as required biomass constituents; deleting '
          'the two MoCo compounds from the objective rescues a moaD knockout to wild-type '
          'growth (0.9825 vs 0.9824 h-1), matching the experimental viability call. Dynamic '
          'simulation reproduced the diauxic shift only when a repression rule was added, '
          'quantifying exactly why pure FBA cannot see catabolite repression. All negative '
          'results are preserved. Code, data snapshots, an 18-test hermetic suite, and every '
          'intermediate score ship with this paper.', BODY),
        P('<b>Keywords:</b> virtual cell, flux balance analysis, gene essentiality, graph '
          'convolutional network, convolutional neural network, dynamic FBA, model ensemble, '
          '<i>Escherichia coli</i>, biomass objective function', BODY),
        PageBreak()]

    # ---------------- introduction ----------------
    story += [P('1. Introduction', H1),
        P('The dream of a whole-cell computational model - a program that takes a genotype and '
          'an environment and returns the phenotype of a living cell - is decades old and still '
          'largely unmet. The strongest existence proof remains the Karr et al. whole-cell model '
          'of <i>Mycoplasma genitalium</i> [5], which coupled 28 submodels to simulate one cell '
          'cycle. For <i>Escherichia coli</i>, the best-developed layer is metabolism: the '
          'iJO1366 reconstruction [3] accounts for 1,366 genes, 2,251 metabolic reactions and '
          '1,136 unique metabolites, and flux balance analysis on it makes phenotypic '
          'predictions that agree well with knockout and growth screens under matched '
          'conditions. Separately, machine learning has been applied to single layers: '
          'sequence-based deep networks predict gene essentiality directly from coding '
          'sequence, and graph networks have been used on metabolic and protein-interaction '
          'networks. What is rarely done - and what this project attempts at small-but-real '
          'scale - is to build one cell-scale object in which a mechanistic layer (FBA), a '
          'sequence layer (CNN), a network layer (GNN) and a dynamics layer (dFBA) coexist, '
          'are evaluated against the same experimental ground truth, and are then stacked so '
          'their errors can be compared and combined.', BODY),
        P('Three questions drive the work. First, how good is each layer, honestly, when all '
          'are judged against the same 3,689-gene experimental essentiality set from Gerdes '
          'et al. 2003 [2]? Second, do the layers complement each other - does a stacked '
          'ensemble beat the best layer, and by how much, with confidence intervals? Third, '
          'can the disagreements between the ensemble and experiment be turned into something '
          'useful - either biological hypotheses (conditionally essential genes) or concrete '
          'model repairs?', BODY),
        P('The answers are, respectively: individually mediocre (AUROC 0.63-0.67); yes for '
          'ranking (ensemble AUROC 0.7225, +0.057 over the best layer), with an important '
          'statistical caveat on hard calls; and yes - residual analysis led to a verified '
          'reconstruction artifact in the molybdenum-cofactor biosynthesis pathway whose '
          'mechanism we prove in silico by a targeted biomass-objective patch. Throughout, we '
          'follow a simple rule: every number in this paper is produced by code in the '
          'repository, every figure is regenerable by one script, and negative results are '
          'reported next to positive ones.', BODY),
        P('1.1 Why essentiality is the right first benchmark', H2),
        P('Gene essentiality is the rare phenotype that is simultaneously binary, '
          'genome-wide, mechanistically interpretable and experimentally measured at scale. '
          'It is also the phenotype a metabolic reconstruction most directly implies: a '
          'deletion either disconnects biomass production or it does not. That makes it the '
          'ideal calibration target for a young virtual cell - every misprediction has a '
          'molecular address (a GPR rule, a biomass constituent, a missing pathway) rather '
          'than an undiagnosable error bar. Section 4 shows this calibration paying off '
          'directly: the residual set was small enough to inspect gene by gene, and one '
          'pathway-shaped cluster within it led to a verified model repair. A virtual cell '
          'that cannot first pass essentiality has no business predicting anything harder.', BODY),
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
