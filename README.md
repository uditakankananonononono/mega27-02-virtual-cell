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

## Dynamics
dFBA (Mahadevan 2002 SOA, MM uptake): unregulated FBA co-utilizes acetate
with glucose (no diauxie); one Boolean catabolite-repression rule restores
sequential utilization. Both runs preserved.

## Layout
- `vcell/` - package: metabolism, data, seqcnn, graphgnn, dynamics, benchmark
- `scripts/` - pipelines: run_ijo1366_scan.py, run_ijo1366_rich_scan.py,
  run_seqcnn.py, run_gnn.py, make_figures.py
- `tests/` - 18 hermetic pytest cases, ~5 s, no network
- `results/` - every intermediate score (CSV/JSON/NPZ)
- `figures/` - 5 paper figures (regenerate: python scripts/make_figures.py)
- `paper/` - VC2_virtual_cell_paper.pdf (20 pp, Times-Roman; regenerate:
  python paper/build_paper.py)
- `data/` - snapshots: BiGG e_coli_core + iJO1366, U00096.3 FASTA/GenBank,
  Gerdes 2003 Table S1

## Verify
```
pip install numpy pandas scipy matplotlib scikit-learn biopython pytest cobra reportlab pypdf pillow torch --index-url https://download.pytorch.org/whl/cpu
python -m pytest tests/ -q
```

Data sources: BiGG Models (e_coli_core, iJO1366), NCBI GenBank U00096.3,
Gerdes et al. 2003 J Bacteriol 185:5673 (genome.wisc.edu supplement).
