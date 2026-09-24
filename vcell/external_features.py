"""Label-free gene features from external resources (STRING v12, UniProt UP000000625, KEGG REST).
All inputs are snapshots under data/external/; no essentiality labels are used here."""
from __future__ import annotations
import gzip, re
import numpy as np, pandas as pd
import networkx as nx

EXT = 'data/external/'


CONTEXT_CHANNELS = ['neighborhood', 'fusion', 'cooccurence', 'coexpression']


def string_features(min_score: int = 700, channels: list[str] | None = None, prefix: str = 'string') -> pd.DataFrame:
    """channels=None: STRING combined_score. Otherwise re-combine only the given evidence channels
    (s = 1 - prod(1 - s_i/1000), no prior correction) - used to exclude literature/database channels
    that could encode prior essentiality knowledge."""
    if channels is None:
        df = pd.read_csv(EXT + 'string_links.txt.gz', sep=' ')
    else:
        df = pd.read_csv(EXT + 'string_links_detailed.txt.gz', sep=' ')
        prod = np.ones(len(df))
        for c in channels:
            prod *= 1 - df[c].values / 1000.0
        df['combined_score'] = np.round((1 - prod) * 1000)
        df = df[df.combined_score > 0]
    df['a'] = df.protein1.str.split('.').str[1]; df['b'] = df.protein2.str.split('.').str[1]
    all_deg = df.groupby('a').combined_score.agg(['count', 'sum'])
    hi = df[df.combined_score >= min_score]
    G = nx.Graph(); G.add_edges_from(zip(hi.a, hi.b))
    core = nx.core_number(G); clus = nx.clustering(G)
    pr = nx.pagerank(G, max_iter=200)
    out = pd.DataFrame({'bnumber': all_deg.index, f'{prefix}_deg_all': all_deg['count'].values,
                        f'{prefix}_wdeg_all': all_deg['sum'].values / 1000.0})
    out[f'{prefix}_deg_hi'] = out.bnumber.map(dict(G.degree())).fillna(0)
    out[f'{prefix}_core'] = out.bnumber.map(core).fillna(0)
    out[f'{prefix}_clust'] = out.bnumber.map(clus).fillna(0)
    out[f'{prefix}_pagerank'] = out.bnumber.map(pr).fillna(0) * 1e4
    return out


def uniprot_features() -> pd.DataFrame:
    u = pd.read_csv(EXT + 'uniprot_ecoli.tsv', sep='\t')
    u['bnumber'] = u['Gene Names (ordered locus)'].astype(str).str.extract(r'(b\d{4})')[0]
    u = u.dropna(subset=['bnumber']).drop_duplicates('bnumber')
    loc = u['Subcellular location [CC]'].fillna('')
    return pd.DataFrame({
        'bnumber': u.bnumber,
        'up_len': np.log1p(u.Length.astype(float)),
        'up_tm': u['Transmembrane'].fillna('').str.count('TRANSMEM'),
        'up_cofactor': u['Cofactor'].notna().astype(int),
        'up_n_ec': u['EC number'].fillna('').map(lambda s: len([x for x in s.split(';') if x.strip()])),
        'up_annot': u['Annotation'].astype(float),
        'up_cytoplasm': loc.str.contains('Cytoplasm').astype(int),
        'up_membrane': loc.str.contains('membrane', case=False).astype(int)})


CORE_PATHWAYS = {'eco00010', 'eco00020', 'eco00030', 'eco00190', 'eco00230', 'eco00240', 'eco00061',
                 'eco00550', 'eco00540', 'eco00790', 'eco00770', 'eco00780', 'eco00730', 'eco00740',
                 'eco00760', 'eco00670', 'eco00860', 'eco00130', 'eco03010', 'eco03020', 'eco03030',
                 'eco00970', 'eco03060', 'eco03070'}


def kegg_features() -> pd.DataFrame:
    k = pd.read_csv(EXT + 'kegg_eco_pathway.tsv', sep='\t', header=None, names=['g', 'p'])
    k['bnumber'] = k.g.str.replace('eco:', ''); k['p'] = k.p.str.replace('path:', '')
    grp = k.groupby('bnumber').p
    return pd.DataFrame({'bnumber': grp.size().index, 'kegg_n_path': grp.size().values,
                         'kegg_core': grp.apply(lambda s: int(bool(set(s) & CORE_PATHWAYS))).values})


def all_features() -> pd.DataFrame:
    f = string_features().merge(string_features(channels=CONTEXT_CHANNELS, prefix='strctx'), on='bnumber', how='outer')
    f = f.merge(uniprot_features(), on='bnumber', how='outer').merge(kegg_features(), on='bnumber', how='outer')
    return f.fillna(0)
