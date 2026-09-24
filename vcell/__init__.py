"""vcell - a small-but-real virtual cell for Escherichia coli K-12 MG1655.

Layers
------
metabolism : genome-scale flux balance analysis on real BiGG models
data       : loaders for experimental ground truth (Gerdes 2003, Keio 2006)
             and the U00096.3 reference genome
seqcnn     : 1D-CNN on coding sequences for essentiality prediction
graphgnn   : GNN on the bipartite metabolite-reaction graph
dynamics   : dynamic FBA (dFBA) time-course simulation
benchmark  : metrics, bootstrap CIs, McNemar tests vs experimental truth
"""
__version__ = "0.1.0"
