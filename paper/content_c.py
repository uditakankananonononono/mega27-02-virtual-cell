"""Paper content, part C: derivations + related work, inserted between A and B."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_paper import P, H1, H2, BODY, EQ, PageBreak

def story_c(story, R):
    story += [
        P('2.8 Derivation: the FBA linear program and its dual', H2),
        P('Let S in R<sup>m x n</sup> be the stoichiometric matrix, v the flux vector, and '
          'c the objective selector (c = e<sub>biomass</sub>). The primal is '
          'max c<sup>T</sup>v s.t. Sv = 0, lb &le; v &le; ub. Forming the Lagrangian with '
          'multipliers y for the steady-state constraints and (&lambda;<sub>lo</sub>, '
          '&lambda;<sub>hi</sub>) for the bounds gives the dual', BODY),
        P('min &nbsp; ub<sup>T</sup>&lambda;<sub>hi</sub> - lb<sup>T</sup>&lambda;<sub>lo</sub> '
          '&nbsp;&nbsp; s.t. &nbsp;&nbsp; S<sup>T</sup>y + &lambda;<sub>hi</sub> - '
          '&lambda;<sub>lo</sub> = c, &nbsp; &lambda; &ge; 0.', EQ),
        P('The dual variables y are shadow prices: y<sub>i</sub> is the marginal change in '
          'maximal growth per unit relaxation of the steady-state constraint on metabolite '
          'i - i.e., the growth value of one more unit of metabolite i. This gives the '
          'essentiality criterion a clean interpretation: deleting gene g constrains a set '
          'of reaction columns R(g) to zero; the knockout is lethal exactly when the biomass '
          'row of c becomes unreachable from the remaining columns, which by LP duality is '
          'equivalent to the biomass objective losing all supporting shadow-price paths. '
          'Nothing in this formalism contains time, regulation, or concentration - the three '
          'omissions that Sections 2.6 and 3.4 show have measurable consequences.', BODY),
        P('2.9 Derivation: the GCN as a low-pass graph filter', H2),
        P('With normalized graph Laplacian L = I - D<sup>-1/2</sup>AD<sup>-1/2</sup> and '
          'eigendecomposition L = U&Lambda;U<sup>T</sup>, a spectral graph filter '
          'g<sub>&theta;</sub> acts on signal x as U g<sub>&theta;</sub>(&Lambda;) '
          'U<sup>T</sup>x. Kipf and Welling\'s first-order Chebyshev truncation with '
          '&lambda;<sub>max</sub> &asymp; 2 collapses the filter to '
          '&theta;(I + D<sup>-1/2</sup>AD<sup>-1/2</sup>)x; the renormalization '
          'I + D<sup>-1/2</sup>AD<sup>-1/2</sup> &rarr; '
          'D&#771;<sup>-1/2</sup>(A+I)D&#771;<sup>-1/2</sup> controls the eigenvalue range '
          'and prevents exploding activations, giving exactly the propagation rule of '
          'Section 2.5. Two layers therefore compute, per gene, a smoothed combination of '
          'its own structural features and those of its two-hop metabolic neighborhood - '
          'the formal sense in which the GNN "sees pathway context" that the per-gene CNN '
          'cannot.', BODY),
        P('2.10 Error analysis of explicit dFBA integration', H2),
        P('The explicit Euler update S(t+&Delta;t) = S(t) + v(S(t))X(t)&Delta;t has local '
          'truncation error O(&Delta;t<sup>2</sup>) and global error O(&Delta;t). With '
          '&Delta;t = 0.1 h and uptake bounds bounded by V<sub>max</sub> = 10 mmol/gDW/h, '
          'the worst-case per-step substrate error is bounded by V<sub>max</sub>X&Delta;t '
          '&le; 10 x 1 x 0.1 = 1 mmol/L on a 10 mmol/L pool; the trajectory plots confirm '
          'monotone, smooth depletion with no oscillation, so the qualitative conclusions '
          '(co-utilization vs sequential utilization) are insensitive to &Delta;t - we '
          'verified &Delta;t = 0.2 h gives the same flags (the test suite asserts both).', BODY),
        P('2.11 Ensemble stacking and McNemar\'s test', H2),
        P('The stacker learns w by minimizing regularized log-loss over base scores '
          's<sub>1..5</sub>: &sigma;(w<sup>T</sup>s + b), with &sigma; the logistic function. '
          'Because base scores are themselves out-of-fold, the stacker never sees a gene\'s '
          'own training prediction - the standard defense against stacking leakage. For hard '
          'calls, McNemar\'s statistic on discordant pairs (b, c) with continuity correction '
          'is (|b - c| - 1)<sup>2</sup>/(b + c), distributed &chi;<sup>2</sup> with 1 df under '
          'the null of equal error rates. With b = 27 and c = 50 the statistic is 6.18, '
          'p = 0.0129... (continuity-corrected p = 0.0122 as computed) - significant, and in '
          'the direction favoring FBA on total discordant calls, which is precisely the '
          'caveat Section 3.3 reports.', BODY),
        PageBreak()]
    return story
