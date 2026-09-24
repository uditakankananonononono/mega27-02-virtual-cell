"""Paper content: Section 4.4 essentiality-provenance catalog (from results/corrected_benchmark.json
and results/essentiality_provenance_audit.csv)."""
import os, sys, json, ast
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_paper import P, H1, H2, BODY, tbl, ROOT, PageBreak


def story_audit(story, R):
    cb = R('corrected_benchmark.json')
    aud = pd.read_csv(os.path.join(ROOT, 'results', 'essentiality_provenance_audit.csv'))
    lab = pd.read_csv(os.path.join(ROOT, 'results', 'aligned_predictions_all_models.csv'))[['bnumber', 'gene_name', 'essential']]
    a = aud.merge(lab, left_on='gene', right_on='bnumber', how='left')
    vc = aud.verdict.value_counts().to_dict()
    t1, t2 = cb['fisher_2x2_[bf_noness,bf_ess],[nf_noness,nf_ess]'], cb['fisher_cofactor_only_vs_rest']
    b0, b1 = cb['fba_baseline_on_aligned'], cb['fba_corrected']
    story += [PageBreak(), P('Appendix E. Essentiality-provenance catalog for iJO1366', H1),
              P('E.1 Question and design', H2),
              P('Section 4.3 showed one case - the molybdopterin cluster - where an FBA essentiality call is produced '
                'by the biomass objective rather than by the network. This appendix asks the general question: for '
                f'every one of the {len(aud)} genes that iJO1366 calls essential on glucose minimal medium, is the call '
                'network_forced or biomass_forced? The rescue test (Appendix C.4, eqs. C3-C5) answers it per gene with '
                'one re-runnable LP experiment. The tool is shipped as python -m vcell audit-biomass MODEL and runs on '
                'any COBRA model file.', BODY),
              P('Prior work. The phenomenon of cofactor-driven false essentiality is known. Bernstein et al. (2023, '
                'Mol Syst Biol) found 21 vitamin/cofactor biosynthesis genes in iML1515 (biotin, R-pantothenate, '
                'thiamin, tetrahydrofolate, NAD+) that the model calls essential and RB-TnSeq finds fit, and fixed them '
                'by manually supplying those cofactors. MEMOTE tests whether biomass precursors are producible in the '
                'wild-type model but does not attribute any knockout\'s lethality to specific compounds. Our '
                'contribution is narrower and procedural: an automated, per-gene, falsifiable attribution. We found '
                'no existing tool that does this; that is a statement about our search, not a proof of absence.', BODY),
              P('E.2 Verdict census', H2)]
    story += tbl([['verdict', 'genes']] + [[k, str(v)] for k, v in vc.items()], 'Table 10. Rescue-test verdicts.')
    story += [P('E.3 A methods catch: raw biomass_forced over-calls artifacts', H2),
              P('A raw biomass_forced label means only that stripping the unproducible compounds rescues growth. For '
                'a gene such as fabD that stops lipid synthesis this is trivially true and says nothing about artifacts: '
                'lipids are really required. The artifact hypothesis only makes sense for compounds a mutant might '
                'plausibly obtain without making them (trace cofactors and vitamins by carry-over or cross-feeding, '
                'the class named by Bernstein et al.). We therefore split biomass constituents into cofactor/vitamin, '
                'inorganic and bulk classes (vcell/compound_classes.py) and label a gene cofactor_only when every '
                'unproducible compound is a cofactor or vitamin.', BODY)]
    st = cb['strata']
    story += tbl([['stratum', 'genes', 'Gerdes essential', 'fraction']] +
                 [[s['stratum'], str(s['count']), str(s['sum']), f"{s['sum']/max(s['count'],1):.2f}"] for s in st],
                 'Table 11. Gerdes 2003 essentiality by provenance stratum.')
    story += [P('E.4 The falsifiable test', H2),
              P('Prediction (Appendix C.4): cofactor_only biomass_forced genes are enriched for experimentally '
                'non-essential genes relative to all other FBA-essential genes. One-sided Fisher exact test (eqs. C18-C19) '
                f'on [[non-ess, ess], [non-ess, ess]] = {t2}: odds ratio {cb["fisher_cofactor_odds"]:.2f}, '
                f'p = {cb["fisher_cofactor_p_one_sided"]:.3g}. For the raw (unstratified) biomass_forced label the table is '
                f'{t1}: odds ratio {cb["fisher_odds_ratio"]:.2f}, p = {cb["fisher_p_one_sided"]:.3g}.', BODY),
              P('E.5 Corrected benchmark', H2)]
    story += tbl([['metric', 'FBA iJO1366 (original)', 'after cofactor_only correction']] +
                 [[k, f'{b0[k]:.3f}', f'{b1[k]:.3f}'] for k in ['auroc', 'auprc', 'precision', 'recall', 'f1', 'mcc']],
                 'Table 12. iJO1366 minimal-medium FBA vs Gerdes 2003 before and after flipping cofactor_only biomass_forced '
                 'genes to non-essential.')
    story += [P(f'Corrected AUROC 95% CI: [{b1["auroc_ci"][0]:.3f}, {b1["auroc_ci"][1]:.3f}]. This correction uses the '
                'audit only, never the Gerdes labels, so it is not fitted to the benchmark; it is still a single '
                'benchmark and should be re-tested on an independent essentiality screen.', BODY)]
    cc = cb['biomass_forced_compound_census']
    story += tbl([['biomass compound', 'genes whose knockout makes it unproducible']] + [[k, str(v)] for k, v in cc.items()],
                 'Table 13. Compound census over biomass_forced genes.')
    a = a.sort_values(['verdict', 'gene'])
    rows = [['b-number', 'gene', 'verdict', 'unproducible compounds', 'Gerdes']]
    for _, r in a.iterrows():
        u = ', '.join(ast.literal_eval(r.unproducible)) if isinstance(r.unproducible, str) else ''
        rows.append([r.gene, str(r.gene_name) if pd.notna(r.gene_name) else '', r.verdict.replace('_no_biomass_link', '*'),
                     u[:60], ('E' if r.essential == 1 else 'N') if pd.notna(r.essential) else 'n/a'])
    story += [P('E.6 Full catalog', H2),
                P('Every audited gene with its verdict, the unproducible biomass compounds (truncated to 60 characters) '
                  'and its Gerdes 2003 call (E essential, N non-essential, n/a not in the aligned set). * marks '
                  'network_forced_no_biomass_link: no individual biomass constituent becomes unproducible, yet growth '
                  'stops, which points to a coupling or energy constraint.', BODY)]
    story += tbl(rows, 'Table 14. Essentiality-provenance catalog.')
    return story
