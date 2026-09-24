"""Aggregate Fitness Browser strain fitness to gene medians (mean over used strains per experiment, median over experiments)."""
import sys, numpy as np, pandas as pd
for org in sys.argv[1:]:
    acc = n = None
    for ch in pd.read_csv(f'data/cross/{org}_strainfit.gz', sep='\t', chunksize=50000, dtype={'locusId': str}):
        ch = ch[((ch.used == True) | (ch.used == 'TRUE')) & ch.locusId.notna()]
        e = [c for c in ch.columns if c.startswith('set')]
        s = ch.groupby('locusId')[e].sum(); c = ch.groupby('locusId')[e].count()
        acc = s if acc is None else acc.add(s, fill_value=0); n = c if n is None else n.add(c, fill_value=0)
    g = acc / n.replace(0, np.nan)
    pd.DataFrame({'locusId': g.index, 'n_exp': g.notna().sum(axis=1).values, 'median_fitness': g.median(axis=1).values}).to_csv(
        f'data/cross/{org}_gene_median_fitness.tsv', sep='\t', index=False)
    print(org, g.shape, flush=True)
