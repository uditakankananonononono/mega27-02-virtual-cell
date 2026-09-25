# VC-2: A Modular Virtual Cell for E. coli K-12

MEGA-PROGRAM-27, item 2. Four-layer in-silico cell benchmarked against the
Gerdes 2003 experimental essentiality set, with a stacked ensemble and a
verified model-repair discovery (MoCo biomass artifact).

## Results (headline, honest)
| Layer | AUROC vs Gerdes 2003 |
|---|---|
| core FBA (125-gene overlap) | 0.630 |
| iJO1366 FBA minimal | 0.666 |
| iJO1366 FBA rich | 0.633 |
| GCN (metabolic graph, pure torch) | 0.660 |
| CNN (CDS sequence) | 0.643 |
| k-mer logistic baseline | 0.628 |
| **Stacked ensemble (OOF)** | **0.7225** |

McNemar caveat: ensemble improves ranking and positive-class F1 (0.456 vs
0.412) but loses total discordant hard calls 27:50 (p=0.012) - see paper 3.3.

## Discovery
Class-A residual analysis localized the molybdenum-cofactor biosynthesis
cluster (moaC/D/E, mobA, moeB) as falsely essential: the iJO1366 core biomass
objective requires MoCo compounds (bmocogdp_c, mobd_c); removing them rescues
a moaD knockout to wild-type growth (0.9825 vs 0.9824 h-1), matching the
experimental viability call. Initial conditional-essentiality hypothesis was
tested and rejected in silico - preserved as a documented negative.

## Tool: `vcell` CLI
```
python -m vcell audit-biomass MODEL.json [--genes b0001 ...] [--out audit.csv]
python -m vcell bernstein-score iML1515 [--supplement btn thf] [--out score.json]
python -m vcell rescue-audit MODEL.xml --supplement btn_c=00780 thf_c=00790,00670 \
    --gene-pathways kegg_eco_pathway.tsv [--no-known-uptake amet_c nad_c thf_c pydx5p_c] [--out rescues.csv]
```
`rescue-audit` implements the supplement-bypass finding (paper H.7,
results/pathway_concordance.json): in-silico rescues whose knocked-out gene is
not on the supplement's own KEGG pathway had median measured fitness -4.24 vs
-0.97 for on-pathway rescues on the Bernstein 2023 benchmark, so off-pathway
rescues are flagged `likely_artifact`. Caveats: confounded with supplement identity (mostly SAM), and not replicated in 8 other bacteria (paper H.8, results/cross_species_*.json: same direction, not significant; a 44-organism test with a model-internal label was uninformative, paper H.17), so treat the flag as an E. coli-validated heuristic.
`--no-known-uptake` adds a `no_known_uptake` column for supplements the organism has no known transporter for;
`vcell.uptake.uptake_systems` computes this from the TCDB substrate table and UniProt TCDB cross-references
(paper H.16, results/uptake_plausibility.json: in E. coli K-12, SAM, NAD, THF and PLP have none; 6 of 7 off-pathway rescues use them).

## Dynamics
dFBA (Mahadevan 2002 SOA, MM uptake): unregulated FBA co-utilizes acetate
with glucose (no diauxie); one Boolean catabolite-repression rule restores
sequential utilization. Both runs preserved.

## Layout
- `vcell/` - package: metabolism, data, seqcnn, graphgnn, dynamics, benchmark
- `scripts/` - pipelines: run_ijo1366_scan.py, run_ijo1366_rich_scan.py,
  run_seqcnn.py, run_gnn.py, make_figures.py
- `tests/` - 23 hermetic pytest cases, ~5 s, no network
- `results/` - every intermediate score (CSV/JSON/NPZ)
- `figures/` - 5 paper figures (regenerate: python scripts/make_figures.py)
- `paper/` - VC2_virtual_cell_paper.pdf (page count varies with source updates;
  actual embedded Times New Roman; regenerate with
  `VCELL_TNR_DIR=/path/to/extracted/fonts python paper/build_paper.py`).
  Extract the private Drive `times32.exe` locally; do not commit or redistribute its TTFs.
- `data/` - snapshots: BiGG e_coli_core + iJO1366, U00096.3 FASTA/GenBank,
  Gerdes 2003 Table S1

## Verify
```
pip install numpy pandas scipy matplotlib scikit-learn biopython pytest cobra reportlab pypdf pillow torch --index-url https://download.pytorch.org/whl/cpu
python -m pytest tests/ -q
```

Data sources: BiGG Models (e_coli_core, iJO1366), NCBI GenBank U00096.3,
Gerdes et al. 2003 J Bacteriol 185:5673 (genome.wisc.edu supplement).

## Strict evidence gates (25 September 2026 audit)
The old `results/tools_manifest.json` describes 40 tool/source entries and 153
data entries. They do **not** establish 40 science/data tools or 120 distinct
accession-level fetched-and-used datasets. Four listed entries are plainly
outside the tool count (three download/list pages and PubMed Central). Forty-four
Fitness Browser organism files share figshare article accession 25236931,
so even the old 124-study grouping falls to an **81-accession upper bound**
before remaining unaccessioned sources and use verification. `results/gate_audit.json`
and the two audit CSVs preserve the exclusions. Both gates are **unmet**; do not
report the inventory lengths as gate completion.

## Learner sensitivity (exploratory)
Frozen 3-fold v2 features were tested with LightGBM, CatBoost and a train-fold-only
imbalanced-learn oversampled logistic regression. SHAP explained CatBoost features
on a separate full-data fit. See `notes/prereg_learner_sensitivity.md`,
`scripts/run_learner_sensitivity.py`, and both result files. OOF AUROC is 0.7773,
0.7847, and 0.7975 respectively vs v2 LR 0.7967. LightGBM's paired AUROC
difference is negative (95% bootstrap CI barely below zero); CatBoost and
oversampling differences include zero. No new biological validation claim. These
four executed libraries raise the strict science-tool candidate ceiling from 36
to 40, but are not a blanket gate pass: earlier candidates still need provenance
checks and some may be excluded. The accession dataset gate remains unmet.

## BiGG cross-model MoCo objective census
`notes/prereg_bigg_model_survey.md` fixes the question and amendment. The
public BiGG v2 catalog listed 108 unique model IDs; all 108 model JSON
snapshots were fetched, parsed, SHA-256 hashed and analysed for biomass
objective constituents by `scripts/run_bigg_model_survey.py`.
`results/bigg_model_survey.csv` gives the accession ID, exact downloaded
source URL, digest, size, parsed counts and result for each model.
**68/108** models include a MoCo-named constituent in a biomass objective;
104/108 had a positive-weight objective annotated in JSON. This is a
descriptive structural finding across 85 stated organism names, not evidence
that all 68 models make a false essentiality call. Model growth was not scored.
These 108 distinct BiGG IDs overlap some models already in the manifest,
and they are metabolic-model accessions, not 108 independent wet-lab screens.
The large raw snapshot collection is retained in three private Drive parts
rather than redistributed in git. `results/bigg_snapshot_delivery.json` records
the exact links, part/archive hashes and reassembly command; per-model
hashes and fetch script allow source revalidation.

## Exploratory pathway check
An offline gseapy overrepresentation test of the top 100 v2 false essentiality
priorities against the 1,249-gene reference universe tested 93 KEGG gene sets.
Four are BH FDR < 0.05, led by folate biosynthesis (q=1.41e-6), cofactor
biosynthesis (q=1.87e-4) and sulfur relay (q=7.35e-4).
`notes/prereg_pathway_enrichment.md`, `scripts/run_pathway_enrichment.py`,
`results/pathway_enrichment.json` and CSV contain the full exploratory result,
source input hashes and caveat. This executes gseapy but does not supply a new
external dataset accession or independent validation.

## Cross-learner ranking and KEGG check
An exploratory offline gseapy enrichment of the top 100 non-essential genes
scored highest by v2 tested 93 KEGG pathways; four pass BH q<0.05, led by
folate biosynthesis (q=1.41e-6). Pingouin Spearman tests of all six learner
pairs show v2 vs oversampled LR rho=0.992 (top-50 Jaccard=0.887), but v2
vs LightGBM rho=0.816 (top-50 Jaccard=0.389): high global agreement need
not preserve the top candidates. Full paths and caveats in `results/` and
`notes/prereg_pathway_enrichment.md`, `notes/prereg_score_agreement.md`.
These are genuine science-package executions, not additional datasets.

## Exploratory error structure
UMAP embedding of standardized frozen v2 features followed by HDBSCAN
finds ten clusters and noise. Two clusters enrich the top-100 false
essentiality priorities (cluster 0: 22/68; cluster 1: 24/100; both pass BH),
while PyOD 20-nearest-neighbor outlier scores are higher for those priorities
(two-sided MWU p=1.90e-20). This is descriptive and reuses the same genes
and labels, not an independent biological result. The complete data, settings,
cluster p/q and negative clusters are in `results/error_structure.json` and
`results/error_structure_gene_scores.csv`; the question was written to
`notes/prereg_error_structure.md` before the computation.

## STRING graph modules
Using only saved high-confidence STRING edges among 1,249 aligned genes,
igraph constructed an 8,714-edge graph and Leidenalg found 58 communities.
After excluding modules with fewer than 20 genes, one of 15 modules was
enriched for top-100 false essentiality priorities after BH correction; three
modules enriched experimental essentiality. The community assignments and
all tested/non-significant rows are in `results/network_modules*.{json,csv}`
(with filenames `network_modules.json` and `network_modules_gene_assignments.csv`).
This is exploratory topology on an existing STRING accession, not another
independent dataset or evidence of cell states.
