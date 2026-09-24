"""Paper content: Appendix F (external evidence + ensemble v2), Appendix G (tools table + dataset manifest)."""
import os, sys, json
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_paper import P, H1, H2, BODY, EQ, tbl, ROOT, PageBreak
from reportlab.lib.units import inch
from reportlab.lib.styles import ParagraphStyle
CELL = ParagraphStyle('cell', fontName='Times-Roman', fontSize=8.5, leading=10)
C = lambda x: P(str(x), CELL)


def story_v2(story, R):
    e = R('ensemble_v2.json')
    fa = pd.read_csv(os.path.join(ROOT, 'results', 'feature_association.csv'))
    story += [PageBreak(), P('Appendix F. External evidence and ensemble v2', H1),
              P('F.1 Motivation', H2),
              P('Figure 6 showed that the v1 ensemble cannot recover essential genes that none of its inputs sees. '
                'The obvious remedy is new, independent evidence about each gene that uses no essentiality labels: '
                'its position in the protein-association network, its physical properties and its pathway context. '
                'We added three external resources, each snapshotted under data/external/ and processed by '
                'vcell/external_features.py: STRING v12 (organism 511145), the UniProt reference proteome '
                'UP000000625, and KEGG REST gene-pathway links.', BODY),
              P('F.2 Features', H2),
              P('From STRING we compute, per gene, total degree and weighted degree over all links, and on the '
                'high-confidence subgraph (score &ge; 700) degree, k-core number, local clustering coefficient and '
                'PageRank (networkx). PageRank is the stationary distribution of', BODY),
              P('PR(i) = (1 - d)/N + d &Sigma;<sub>j &isin; N(i)</sub> PR(j) / deg(j), &nbsp; d = 0.85    (F1)', EQ),
              P('Because STRING\'s experimental, database and text-mining channels may encode prior knowledge about '
                'essentiality (well-studied essential genes accumulate literature links), we also rebuild the network '
                'from the genomic-context and coexpression channels only, combining channel scores s<sub>c</sub> as', BODY),
              P('s = 1 - &Pi;<sub>c &isin; {neighborhood, fusion, cooccurrence, coexpression}</sub> (1 - s<sub>c</sub>/1000)    (F2)', EQ),
              P('which is STRING\'s noisy-OR combination without its prior correction. Features from this network carry '
                'the prefix strctx. From UniProt we take log protein length, transmembrane segment count, presence of '
                'a cofactor annotation, number of EC numbers, annotation score and cytoplasm/membrane localization. '
                'From KEGG we take the number of pathways and membership in a fixed list of 24 core pathways '
                '(central carbon, nucleotide, amino-acid, cofactor, cell-wall, ribosome, replication, RNA polymerase, '
                'aminoacyl-tRNA, protein export).', BODY),
              P('F.3 Univariate associations', H2),
              P('Each feature was tested against the Gerdes labels with a two-sided Mann-Whitney U test and the '
                'p-values corrected with Benjamini-Hochberg (statsmodels):', BODY),
              P('q<sub>(k)</sub> = min<sub>j &ge; k</sub> min(1, m p<sub>(j)</sub> / j)    (F3)', EQ)]
    rows = [['feature', 'AUROC', 'median ess.', 'median non-ess.', 'p', 'q (BH)']]
    for _, r in fa.iterrows():
        rows.append([r.feature, f'{r.auroc:.3f}', f'{r.median_ess:.3g}', f'{r.median_non:.3g}', f'{r.p:.2g}', f'{r.q_bh:.2g}'])
    story += tbl(rows, 'Table 15. Univariate association of external features with Gerdes 2003 essentiality.')
    top = fa.iloc[0]
    story += [P(f'Finding. The single strongest external predictor is the degree in the leakage-controlled '
                f'genomic-context/coexpression network ({top.feature}, AUROC {top.auroc:.3f}, q = {top.q_bh:.1g}). It is '
                'stronger than the degree in STRING\'s full combined network, which includes literature and database '
                'evidence. The signal that essential genes are network hubs therefore does not come from literature '
                'bias. Essential proteins are also shorter than non-essential ones (AUROC 0.383 for length, i.e. 0.617 '
                'for shortness), and less often membrane proteins.', BODY),
              P('F.4 Ensemble v2 and an honest decomposition of the gain', H2)]
    names = {'v1_published_oof': 'v1 as published (unscaled stack)', 'v1_stack_reproduced': 'v1 re-fit with feature scaling',
             'ext_only_lr': 'external features only (LR)', 'v2_lr': 'v2 LR (pre-specified primary)', 'v2_gbm': 'v2 gradient boosting',
             'v2_avg': 'v2 rank-average LR+GBM', 'v2clean_lr': 'v2 leakage-controlled LR', 'v2clean_gbm': 'v2 leakage-controlled GBM',
             'v2clean_avg': 'v2 leakage-controlled rank-average'}
    rows = [['model', 'OOF AUROC', 'OOF AUPRC']] + [[names[k], f"{e[k]['auroc']:.3f}", f"{e[k]['auprc']:.3f}"] for k in names if k in e]
    story += tbl(rows, 'Table 16. Out-of-fold performance vs Gerdes 2003 (1,249 genes, 3-fold stratified CV, seed 7).')
    p1, p2, p3 = e['paired_v2lr_vs_v1published'], e['paired_v2lr_vs_v1reproduced'], e['paired_v2cleanlr_vs_v1reproduced']
    story += [P('The improvement over the published v1 number has two sources, and we report them separately. First, a '
                'defect in v1: the k-mer scores span only 0.46-0.54, so an unscaled logistic stack gives them almost '
                f'no weight. Standardizing features before stacking lifts v1 from {e["v1_published_oof"]["auroc"]:.3f} to '
                f'{e["v1_stack_reproduced"]["auroc"]:.3f} with no new information. Second, the external evidence: '
                f'v2 LR improves on the re-fit v1 by {p2["mean"]:+.3f} AUROC (paired bootstrap 95% CI '
                f'[{p2["ci95"][0]:.3f}, {p2["ci95"][1]:.3f}], one-sided p = {p2["p_le_0"]:.3f}). With the leakage-'
                f'controlled feature set the gain is {p3["mean"]:+.3f} [{p3["ci95"][0]:.3f}, {p3["ci95"][1]:.3f}]. That is '
                f'smaller but still above zero. The primary model was fixed as LR before seeing results. The rank-average '
                'models score slightly higher but were chosen afterwards and are reported as secondary. Averaged over '
                f'five CV seeds, v2 LR scores {e.get("best_5seed_auroc", float("nan")):.3f}.', BODY),
              P(f'A further feature block (v3) was tested and did not help: PaxDb integrated protein abundance and the '
                f'codon adaptation index (CAI, Sharp and Li 1987; w<sub>c</sub> = f<sub>c</sub>/max<sub>c\' syn c</sub> f<sub>c\'</sub>, '
                f'CAI = exp(mean log w) against the ribosomal-protein reference set) with GC3. Abundance alone is '
                f'associated with essentiality (Table 15), but adding it changes v2 LR by {e["paired_v3lr_vs_v2lr"]["mean"]:+.3f} '
                f'AUROC [{e["paired_v3lr_vs_v2lr"]["ci95"][0]:.3f}, {e["paired_v3lr_vs_v2lr"]["ci95"][1]:.3f}]. Its signal is '
                'already captured by network degree and localization. v2 stays the primary model.', BODY),
              P('What this is not: it is not a state-of-the-art claim. Published machine-learning essentiality '
                'predictors for E. coli use different gene sets, labels and splits, so we have no like-for-like number '
                'to beat. The claim is only the internal, paired comparison above.', BODY)]
    m = R('tools_manifest.json')
    story += [PageBreak(), P('Appendix G. External tools and dataset manifest', H1),
              P(f'This appendix lists every external tool ({len(m["tools"])}) and every accession-level dataset '
                f'({len(m["datasets"])}) actually used in this work. A condition matrix from one study counts as one '
                'dataset (the Keio_ML9 RB-TnSeq matrix has 92 conditions but is listed once). Nothing is listed that '
                'was only considered.', BODY)]
    story += tbl([['#', 'tool', 'kind', 'used for']] + [[str(i + 1)] + [C(x) for x in t] for i, t in enumerate(m['tools'])],
                 'Table 17. External tools.', widths=[0.3 * inch, 2.3 * inch, 0.9 * inch, 3.2 * inch])
    story += tbl([['#', 'dataset', 'type']] + [[str(i + 1)] + [C(x) for x in d] for i, d in enumerate(m['datasets'])],
                 'Table 18. Dataset manifest.', widths=[0.3 * inch, 4.3 * inch, 2.0 * inch])
    return story
