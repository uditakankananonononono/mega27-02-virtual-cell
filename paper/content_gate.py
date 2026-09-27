"""Current strict gate accounting and follow-up science (no inflated manifest count)."""
import json,os
from build_paper import P,H1,H2,BODY,tbl,PageBreak,ROOT

def story_gate(story,R):
    ga=R('gate_audit.json');bg=R('bigg_model_survey.json');delta=R('bigg_accession_delta.json')
    lr=R('learner_sensitivity.json');en=R('pathway_enrichment.json');ag=R('score_agreement.json')
    story += [PageBreak(),P('Appendix I. New model survey and gate accounting',H1),
              P('I.1 Cross-model biomass-objective survey',H2),
              P(f"We fetched {bg['fetched_and_analysed']} distinct model JSON snapshots from the public BiGG v2 index "
                f"({bg['source_url']}); each model accession, source URL and SHA-256 is in results/bigg_model_survey.csv. "
                f"In {len(bg['moco_biomass_models'])} of 108 snapshots, the positive-weight biomass objective consumes "
                "a metabolite with a MoCo or molybdenum-cofactor name or ID. Four models lacked a positive-weight "
                "objective annotated in their JSON and remain in the denominator. This is structural prevalence only: "
                "we did not run knockout rescues in every model and do not call those 68 objectives wrong. Models "
                "cover heterogeneous organisms and media. The raw snapshots are privately archived in three Drive "
                "parts; results/bigg_snapshot_delivery.json gives reassembly checksums. They are model datasets, "
                "not independent wet-lab screens.",BODY),
              P('I.2 Frozen learner sensitivity',H2)]
    rows=[['Learner','OOF AUROC','OOF AUPRC','AUROC difference vs v2 LR (95% paired CI)']]
    for key in ('v2_lr','lightgbm','catboost','imblearn_oversampled_lr'):
        m=lr['metrics'][key]
        dif='reference' if key=='v2_lr' else f"{m['delta_auroc_vs_v2_lr']:+.3f} [{m['paired_ci95'][0]:+.3f}, {m['paired_ci95'][1]:+.3f}]"
        rows.append([key,f"{m['auroc']:.3f}",f"{m['auprc']:.3f}",dif])
    story+=tbl(rows,'Table I1. Frozen 3-fold out-of-fold comparison; LightGBM, CatBoost and imbalanced-learn were executed, and SHAP summarized the full-data CatBoost model.')
    top=en['top10'][0]
    story += [P('I.3 Orthogonal checks (exploratory)',H2),
              P(f"Offline gseapy overrepresentation tested {en['n_tested_pathways']} saved KEGG gene sets against "
                f"the top 100 non-essential genes ranked highest by v2. Four pathways pass BH q < 0.05; "
                f"{top['pathway_name']} ranks first (q={top['fdr_bh']:.2g}). This reflects enrichment among "
                "model errors, not an independent essentiality experiment. Pingouin Spearman correlation and "
                "top-50 Jaccard quantify how much alternative learners agree on priorities:",BODY)]
    rr=[['Scores','Spearman','Top-50 Jaccard']]
    for x in ag['comparisons']:
        if x['a']=='v2_lr':rr.append([x['a']+' vs '+x['b'],f"{x['spearman']:.3f}",f"{x['top50_jaccard']:.3f}"])
    story+=tbl(rr,'Table I2. Out-of-fold score agreement; scores share genes, labels and feature construction.')
    st=R('error_structure.json')
    significant=[x for x in st['clusters'] if x['cluster']>=0 and x['fdr_bh']<0.05]
    story += [P('I.4 Feature-space error structure (exploratory)',H2),
              P(f"UMAP embedding of the {st['n_features']} frozen v2 features and HDBSCAN yielded "
                f"{st['n_clusters_excluding_noise']} clusters plus noise. Two clusters enriched the top-100 "
                f"false essentiality priorities: cluster {significant[0]['cluster']} contains "
                f"{significant[0]['top100_false_priorities']}/{significant[0]['size']} (BH q={significant[0]['fdr_bh']:.2g}); "
                f"cluster {significant[1]['cluster']} contains {significant[1]['top100_false_priorities']}/"
                f"{significant[1]['size']} (q={significant[1]['fdr_bh']:.2g}). PyOD kNN outlier scores "
                f"were higher among these priorities (two-sided MWU p={st['pyod_outlier_score']['two_sided_mwu_p']:.2g}). "
                "No cell state or independent replication follows from this geometry. The complete cluster and "
                "gene-level results include every null cluster.",BODY),
              P('I.5 Network localization (exploratory). On the saved STRING v12 network, igraph and Leidenalg found 58 modules.',H2),
              P('The network retained 8,714 edges at score 700 or above among 1,249 aligned genes. We tested the 15 '
                'modules with at least 20 genes. Module 3 contained 17 of the top-100 false essentiality '
                'priorities among 96 genes (one-sided Fisher q=0.014 after BH), while three modules '
                'enriched experimental essentiality. The complete null-module table and assignments '
                'remain in results/network_modules.json and results/network_modules_gene_assignments.csv. '
                'This is a descriptive graph check, not a causal pathway or independent validation.',BODY),
              P('In 20,000 degree-stratified Numba permutations, module 3 had 17 false priorities '
                'versus 7.88 expected (one-sided p=0.00110; max-statistic FWER p=0.0252 over 15 modules). '
                'The JIT result matched a NumPy check. This is post-selection sensitivity on reused data, '
                'not independent validation; all modules are in results/network_degree_null.json.',BODY),
              P('A scikit-network Louvain partition matched Leiden moderately (ARI 0.582). Its '
                'closest module had 12 false priorities in 81 genes (degree-matched p=0.0298), '
                'but BH q=0.174 across 15 tested modules: a cross-algorithm null. '
                'See results/graph_algorithm_sensitivity.json.',BODY),
              P('I.6 Independent LP verification and gate accounting',H2),
              P('HiGHS rebuilt the iJO1366 LP independently of GLPK: wild type 0.9823718/h, '
                'moaD knockout zero; after removing two MoCo biomass coefficients, both grow at '
                '0.9824852/h. Maximum mass-balance residual is below 1.2e-11. This verifies the '
                'numeric model result, not an in-vivo phenotype (results/highs_moco_rescue.json).',BODY),
              P('The SymPy GPR audit checked 4,181 reaction-gene knockouts across 2,123 reactions: '
                'zero Boolean-rule disagreements with COBRApy; moaD disabled the same two reactions. '
                'This validates model rules, not cells (results/gpr_symbolic.json).',BODY),
              P('The tool audit links 48 code/result candidates: eight thin checks or overlapping '
                'wrappers and nine original exclusions do not count. The conservative distinct floor '
                'is 40. KEGG source bytes and TCDB source rows match current public tables, though '
                'historical acquisition times are unknown (tool_execution_summary.json; '
                'tool_source_snapshot_verification.json).',BODY),
              P('Dataset count: the owner accepted BiGG model accessions toward the 120-dataset gate '
                '(25 September WhatsApp, notes/gate_count_provenance_review_20260927.md). The source ledgers '
                'report 108 unique BiGG IDs and a conservative 15-source floor from 19 PaxDb files, yielding '
                '123 under that rule. A separate local check matched all 108 model and 19 PaxDb file hashes '
                'and the 19 saved PaxDb scores (results/model_inclusive_inventory_local_verification.json). '
                'Owner permission alone does not certify the count: provider freshness was not rechecked, '
                'and eight other provisional manifest rows remain unresolved. Nor are these 123 independent '
                'wet-lab studies. The old 124 count double-counted 44 files under one figshare accession. '
                'The original G2 and 50+ substantive body-page gates remain failed.',BODY)]
    story += tbl([['Gate','Observed','Decision'],
                  ['Science/data tools','40 substantial distinct uses / 48 linked','Pass: thin rows excluded; sources rechecked'],
                  ['Accessioned datasets','108 BiGG models + 15 study-collapsed PaxDb','Count basis approved; 123 inventory requires own verification'],
                  ['Independent wet-lab screens','Far fewer than 120','No claim of 120 experimental studies']],
                 'Table I3. 25 September model-inclusive inventory; owner-approved counting rule is separate from source verification.')

    return story
