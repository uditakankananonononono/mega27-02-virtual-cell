"""Paper content: Appendix H - head-to-head on the Bernstein et al. 2023 RB-TnSeq benchmark."""
import os, sys, json, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_paper import P, H1, H2, BODY, EQ, tbl, fig, ROOT, PageBreak
from reportlab.lib.units import inch

DESC = {'iML1515_base': 'iML1515 as released (their pipeline)',
        'iML1515_bernstein_vitamins': 'iML1515 + 5 hand-picked vitamins (Bernstein Part 6)',
        'iML1515_auto_cofactors': 'iML1515 + 19 audit-derived cofactors (ours, label-free)',
        'iML1515_bernstein_allcorr': 'Bernstein all-corrections model (published SOTA)',
        'iML1515_allcorr_plus_auto': 'all-corrections + audit-derived cofactors (ours)',
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
                'they release an all-corrections model. It is the published state of the art for this task.', BODY),
              P('H.2 Protocol', H2),
              P('To rule out any protocol drift, we do not re-implement their benchmark. scripts/run_bernstein_benchmark.py '
                'executes the function cells of their released notebook (MIT licence, vendored under '
                'external/E_coli_GEM_validation) for model loading, gene and carbon-source matching, the BW25113 strain '
                'correction, removal of non-conditionally essential genes, carbon-source mapping, FBA simulation '
                '(carbon uptake -10, other media -1000 mmol/gDW/h, growth threshold 0.001) and the metric. We changed '
                'only one thing: their model_adjustments function reads a notebook-global variable, which we bind '
                'explicitly. The metric is their precision-recall AUC, where labels are the model growth calls, the '
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
    return story
