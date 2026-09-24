"""Paper content: Appendix A - mathematical derivations and proofs (expansion to 50 pages)."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_paper import P, H1, H2, BODY, EQ, PageBreak


def eq(t):
    return P(t, EQ)


def story_math(story, R):
    s = story
    s += [PageBreak(), P('Appendix A. Mathematical derivations and proofs', H1),
          P('This appendix collects every formula used by the VC-2 pipeline, each with a derivation '
            'or proof sketch, so that every number in the paper can be traced from definition to code. '
            'Equations are numbered (A1)-(A24). Notation: S is the m x n stoichiometric matrix, v the '
            'flux vector, lb and ub the bounds, c the objective selector, and G the gene set with '
            'gene-protein-reaction (GPR) map R(g), the set of reactions whose Boolean GPR evaluates to '
            'false when gene g is deleted.', BODY)]

    s += [P('A.1 Monotonicity of the knockout optimum', H2),
          P('Let F = {v : Sv = 0, lb &le; v &le; ub} and z* = max{c<sup>T</sup>v : v in F}. Deleting gene '
            'g imposes v<sub>j</sub> = 0 for every j in R(g), producing F<sub>g</sub> = F &cap; '
            '{v<sub>j</sub> = 0, j in R(g)}.', BODY),
          eq('z*<sub>g</sub> = max { c<sup>T</sup>v : v in F<sub>g</sub> } &le; z*    (A1)'),
          P('Proof. F<sub>g</sub> &sube; F, and the maximum of a linear function over a subset cannot '
            'exceed the maximum over the superset. Corollary: the growth-fraction score '
            'r<sub>g</sub> = z*<sub>g</sub>/z* lies in [0, 1] up to solver tolerance, which justifies '
            'the essentiality score 1 - r<sub>g</sub> used throughout. Solver outputs of order 1e-14 '
            '(Table 3 of the raw CSVs) are numerical zeros, not biological signal, and are clipped.', BODY),
          P('A.2 Double deletions are submodular-free', H2),
          P('For genes g and h, F<sub>gh</sub> &sube; F<sub>g</sub> &cap; F<sub>h</sub>, hence', BODY),
          eq('z*<sub>gh</sub> &le; min(z*<sub>g</sub>, z*<sub>h</sub>)    (A2)'),
          P('but no lower bound in terms of single deletions exists: synthetic lethality '
            '(z*<sub>g</sub> = z*<sub>h</sub> = z*, z*<sub>gh</sub> = 0) is exactly the case of two '
            'isozymes or parallel routes. This is why single-deletion FBA cannot see redundancy loss, '
            'and why isozyme GPR errors (Bernstein et al., 2023) dominate FBA false negatives.', BODY)]

    s += [P('A.3 Producibility of a biomass constituent', H2),
          P('For metabolite i, add a demand column d<sub>i</sub> = -e<sub>i</sub> to S and solve', BODY),
          eq('&pi;<sub>i</sub>(F<sub>g</sub>) = max { u : S v + d<sub>i</sub> u = 0, v in bounds, v<sub>j</sub> = 0 for j in R(g), u &ge; 0 }    (A3)'),
          P('Metabolite i is producible in the knockout iff &pi;<sub>i</sub> &gt; &epsilon; '
            '(we use &epsilon; = 1e-6). The unproducible set is U<sub>g</sub> = {i in B : '
            '&pi;<sub>i</sub>(F<sub>g</sub>) &le; &epsilon;}, where B is the set of consumed biomass '
            'metabolites.', BODY),
          P('A.4 The biomass rescue test and its soundness', H2),
          P('Write the biomass column as b = -&Sigma;<sub>i in B</sub> &beta;<sub>i</sub> e<sub>i</sub> + '
            'p (p collects products such as ADP, Pi, H+). The rescued column strips exactly '
            'U<sub>g</sub>:', BODY),
          eq('b<sup>(g)</sup> = b + &Sigma;<sub>i in U<sub>g</sub></sub> &beta;<sub>i</sub> e<sub>i</sub>    (A4)'),
          eq('verdict(g) = biomass_forced  iff  z*<sub>g</sub>(b<sup>(g)</sup>) &ge; &tau; z*,   else network_forced    (A5)'),
          P('Proposition 1 (soundness of biomass_forced). If verdict(g) = biomass_forced, then the '
            'lethality of g in the original model is entirely attributable to the requirement for '
            'U<sub>g</sub>: the knockout network can carry biomass flux at &ge; &tau; of wild type for '
            'every other constituent simultaneously. Proof: the rescued LP is feasible with value '
            '&ge; &tau; z*, and its optimal v is a flux distribution in F<sub>g</sub> that '
            'produces each i in B \\ U<sub>g</sub> at rate &beta;<sub>i</sub> z*<sub>g</sub>(b<sup>(g)</sup>). '
            'Since the original knockout optimum is zero while the only difference between the two LPs '
            'is the U<sub>g</sub> coefficients, U<sub>g</sub> is a sufficient cause of lethality.', BODY),
          P('Proposition 2 (network_forced is conservative). If verdict(g) = network_forced, then '
            'even an objective that demands none of the unproducible compounds cannot be met; the '
            'lesion is in shared precursor or energy supply. Proof: by (A1) applied to the rescued '
            'objective. Note that network_forced does not claim the gene is truly essential in vivo - '
            'only that no biomass relabelling rescues it.', BODY),
          P('Proposition 3 (minimality). U<sub>g</sub> is minimal among sets defined by single-compound '
            'producibility: removing any i in U<sub>g</sub> from the strip set leaves a compound with '
            '&pi;<sub>i</sub> = 0 in the objective, and a biomass reaction consuming an unproducible '
            'compound at positive coefficient forces z* = 0. Hence every element of U<sub>g</sub> is '
            'necessary for rescue.', BODY),
          P('Falsifiable prediction. If biomass_forced labels reflect model artifacts rather than '
            'biology, genes so labelled should be enriched for experimentally non-essential genes relative '
            'to network_forced genes. This is tested with the one-sided Fisher exact test (A18) in '
            'Section 4.4.', BODY)]

    s += [P('A.5 AUROC as a Mann-Whitney statistic', H2),
          P('For scores s, positives P (essential, n<sub>1</sub>) and negatives N (n<sub>0</sub>):', BODY),
          eq('AUROC = (1/(n<sub>1</sub> n<sub>0</sub>)) &Sigma;<sub>i in P</sub> &Sigma;<sub>j in N</sub> [ 1(s<sub>i</sub> &gt; s<sub>j</sub>) + &frac12; 1(s<sub>i</sub> = s<sub>j</sub>) ]    (A6)'),
          P('Derivation: the ROC curve is (FPR(t), TPR(t)) as t sweeps; integrating TPR dFPR over '
            'the step function gives, for each negative j, the fraction of positives scoring above it, '
            'which summed over j and normalized is (A6). Ties receive half credit because a tied pair '
            'lies on a diagonal ROC segment. Consequence for FBA: binary FBA scores have massive ties, '
            'so AUROC collapses to (TPR + TNR)/2 at the single threshold, which explains why pure FBA '
            'AUROCs cluster near 0.63-0.67.', BODY),
          eq('AUROC<sub>binary</sub> = &frac12;(TPR + 1 - FPR)    (A7)'),
          P('A.6 Average precision', H2),
          eq('AP = &Sigma;<sub>k</sub> (R<sub>k</sub> - R<sub>k-1</sub>) P<sub>k</sub>    (A8)'),
          P('with P<sub>k</sub>, R<sub>k</sub> precision and recall at the k-th threshold. The '
            'baseline for a random ranker equals the prevalence &pi; (0.172 here), so AUPRC must '
            'always be read against &pi;. The lift AP/&pi; is reported alongside.', BODY),
          P('A.7 Matthews correlation coefficient', H2),
          eq('MCC = (TP TN - FP FN) / sqrt((TP+FP)(TP+FN)(TN+FP)(TN+FN))    (A9)'),
          P('MCC equals the Pearson correlation of the two binary indicator vectors (proof: expand '
            'the covariance of Bernoulli variables), so it is symmetric in classes and insensitive to '
            'prevalence in the sense that a constant predictor scores 0.', BODY)]

    s += [P('A.8 Percentile bootstrap confidence intervals', H2),
          P('For statistic T and B resamples of genes with replacement, T*<sub>b</sub> = T(X*<sub>b</sub>). '
            'The 95% interval is', BODY),
          eq('CI<sub>95</sub> = [ Q<sub>0.025</sub>(T*), Q<sub>0.975</sub>(T*) ]    (A10)'),
          P('Validity rests on the bootstrap distribution of T* - T approximating that of T - &theta;; '
            'for U-statistics such as AUROC this holds asymptotically (Hall, 1992). Resamples with a '
            'single class are discarded (the code reports n_boot_used). Genes are the resampling unit, '
            'so the interval reflects gene-set uncertainty, not experimental replicate noise.', BODY),
          P('A.9 Paired comparison of two classifiers: McNemar', H2),
          P('On paired hard calls, let b = number of genes model A gets right and B wrong, c the reverse. '
            'Under H<sub>0</sub> (equal error rates) each discordant gene is a fair coin, b ~ Bin(b+c, &frac12;):', BODY),
          eq('p = 2 &Sigma;<sub>k=0</sub><sup>min(b,c)</sup> C(b+c, k) 2<sup>-(b+c)</sup>    (A11)'),
          eq('&chi;<sup>2</sup> = (|b - c| - 1)<sup>2</sup> / (b + c)    (A12)'),
          P('(A12) is the continuity-corrected normal approximation, derived from the variance '
            '(b+c)/4 of Bin(b+c, &frac12;). With b = 27, c = 50 the exact test gives p = 0.012 against '
            'the ensemble at the 0.5 cut, which is why the paper reports the ranking gain and the '
            'hard-call loss together.', BODY),
          P('A.10 DeLong-free AUROC difference by paired bootstrap', H2),
          eq('&Delta; = AUROC<sub>A</sub> - AUROC<sub>B</sub>, &nbsp; CI from paired resamples (same gene indices for A and B)    (A13)'),
          P('Pairing preserves the correlation between models scored on identical genes and is '
            'strictly tighter than comparing two independent intervals.', BODY)]

    s += [P('A.11 Stacked generalization without leakage', H2),
          P('Base scores are produced out-of-fold: for fold k with held-out set I<sub>k</sub>, model '
            'f<sup>(-k)</sup> is trained on the complement and scores I<sub>k</sub>. The meta-learner is', BODY),
          eq('h(x) = &sigma;( w<sub>0</sub> + &Sigma;<sub>m</sub> w<sub>m</sub> f<sub>m</sub><sup>(-k(x))</sup>(x) )    (A14)'),
          P('Because each base prediction for gene x was made without seeing x, the meta-features are '
            'distributed as test-time features, and the meta-learner\'s own OOF evaluation is unbiased '
            'up to fold-to-fold variance (Wolpert, 1992). Non-learned FBA scores need no folding.', BODY),
          P('A.12 Logistic loss and its gradient', H2),
          eq('L(w) = -&Sigma;<sub>i</sub> [ y<sub>i</sub> log p<sub>i</sub> + (1-y<sub>i</sub>) log(1-p<sub>i</sub>) ] + &lambda;||w||<sup>2</sup>/2, &nbsp; &nabla;L = X<sup>T</sup>(p - y) + &lambda;w    (A15)'),
          P('Derivation: d/dz log &sigma;(z) = 1 - &sigma;(z). The Hessian X<sup>T</sup>diag(p(1-p))X + '
            '&lambda;I is positive definite, so the meta-learner has a unique optimum.', BODY),
          P('A.13 Class-weighted loss for the CNN and GCN', H2),
          eq('L<sub>w</sub> = -&Sigma;<sub>i</sub> [ &omega;<sub>1</sub> y<sub>i</sub> log p<sub>i</sub> + &omega;<sub>0</sub>(1-y<sub>i</sub>) log(1-p<sub>i</sub>) ], &nbsp; &omega;<sub>c</sub> = n/(2 n<sub>c</sub>)    (A16)'),
          P('With &omega;<sub>c</sub> = n/(2n<sub>c</sub>) each class contributes equal total weight, which '
            'is the minimizer of expected balanced error for a calibrated model.', BODY)]

    s += [P('A.14 Receptive field of the sequence CNN', H2),
          eq('RF<sub>L</sub> = 1 + &Sigma;<sub>l=1</sub><sup>L</sup> (k<sub>l</sub> - 1) &Pi;<sub>j&lt;l</sub> s<sub>j</sub>    (A17)'),
          P('For kernel sizes k<sub>l</sub> and strides s<sub>j</sub> (pooling counted as stride). The '
            'receptive field bounds the longest motif the network can represent before global pooling; '
            'global max pooling then makes the gene score invariant to motif position.', BODY),
          P('A.15 Fisher exact test for provenance enrichment', H2),
          P('Let the 2 x 2 table be [[a, b], [c, d]] with rows biomass_forced / network_forced and columns '
            'non-essential / essential in vivo. Conditioning on margins, a follows the hypergeometric law', BODY),
          eq('P(A = a) = C(a+b, a) C(c+d, c) / C(n, a+c)    (A18)'),
          eq('p<sub>one-sided</sub> = &Sigma;<sub>k &ge; a</sub> P(A = k), &nbsp; OR = ad / bc    (A19)'),
          P('The one-sided alternative encodes the prediction of A.4: biomass_forced calls are '
            'more often wrong (non-essential in vivo) than network_forced calls.', BODY),
          P('A.16 Graph convolution as normalized propagation', H2),
          eq('H<sup>(l+1)</sup> = &sigma;( D&#771;<sup>-1/2</sup> A&#771; D&#771;<sup>-1/2</sup> H<sup>(l)</sup> W<sup>(l)</sup> ), &nbsp; A&#771; = A + I    (A20)'),
          P('The symmetric normalization has spectrum in (-1, 1], so repeated propagation is '
            'contractive and does not explode; two layers aggregate the 2-hop gene neighbourhood via '
            'shared metabolites.', BODY),
          P('A.17 Dynamic FBA', H2),
          eq('dX/dt = &mu;(t) X, &nbsp; dS<sub>k</sub>/dt = v<sub>k</sub>(t) X, &nbsp; &mu;(t) = max c<sup>T</sup>v s.t. Sv = 0, v<sub>k</sub> &ge; -V<sub>max,k</sub> S<sub>k</sub>/(K<sub>m,k</sub> + S<sub>k</sub>)    (A21)'),
          eq('X<sub>t+&Delta;t</sub> = X<sub>t</sub> e<sup>&mu;&Delta;t</sup>, &nbsp; S<sub>t+&Delta;t</sub> = S<sub>t</sub> + v (X<sub>t</sub>/&mu;)(e<sup>&mu;&Delta;t</sup> - 1)    (A22)'),
          P('(A22) is the exact solution of (A21) with fluxes frozen over the step (the static-'
            'optimization approach of Mahadevan et al., 2002); it avoids the negative-concentration '
            'artifacts of forward Euler. The Boolean catabolite-repression overlay sets '
            'V<sub>max,lac</sub> = 0 while S<sub>glc</sub> &gt; &theta;, producing diauxie.', BODY),
          P('A.18 Threshold sensitivity and bimodality', H2),
          eq('E(&tau;) = 1(r<sub>g</sub> &lt; &tau;), &nbsp; &part;|E(&tau;)|/&part;&tau; = &Sigma;<sub>g</sub> &delta;(r<sub>g</sub> - &tau;)    (A23)'),
          P('Because FBA growth fractions are almost all exactly 0 or 1, the density of r<sub>g</sub> '
            'in (0.001, 0.5) is nearly empty and |E(&tau;)| is flat across two orders of magnitude of '
            '&tau; - the robustness result of Section 3.', BODY),
          P('A.19 Lift and expected calibration', H2),
          eq('lift = AP / &pi;, &nbsp; ECE = &Sigma;<sub>b</sub> (|B<sub>b</sub>|/n) |acc(B<sub>b</sub>) - conf(B<sub>b</sub>)|    (A24)'),
          P('Lift normalizes AUPRC for prevalence, which differs between Gerdes (17.2%) and the PEC '
            'subset, making cross-yardstick comparisons meaningful.', BODY)]
    return s
