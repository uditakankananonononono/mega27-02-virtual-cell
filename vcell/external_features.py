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
    f = f.merge(paxdb_features(), on='bnumber', how='outer').merge(codon_features(), on='bnumber', how='outer')
    f = f.merge(uniprot_features(), on='bnumber', how='outer').merge(kegg_features(), on='bnumber', how='outer')
    return f.fillna(0)


def paxdb_features() -> pd.DataFrame:
    """PaxDb integrated whole-organism protein abundance (ppm), log10-transformed."""
    p = pd.read_csv(EXT + 'paxdb_511145_integrated.txt', sep='\t', comment='#', header=None, names=['name', 'sid', 'ppm'])
    p['bnumber'] = p.sid.str.split('.').str[1]
    return pd.DataFrame({'bnumber': p.bnumber, 'pax_log_ppm': np.log10(p.ppm.astype(float) + 0.01)})


def codon_features(gb_path: str = 'data/U00096.3.gb') -> pd.DataFrame:
    """Codon Adaptation Index (Sharp & Li 1987) against ribosomal-protein reference set, plus GC3.
    w_c = f_c / max_{c' syn c} f_c'; CAI = exp(mean log w)."""
    from Bio import SeqIO
    from Bio.Data.CodonTable import standard_dna_table
    fwd = standard_dna_table.forward_table
    syn = {}
    for c, aa in fwd.items():
        syn.setdefault(aa, []).append(c)
    cds = {}
    for rec in SeqIO.parse(gb_path, 'genbank'):
        for f in rec.features:
            if f.type == 'CDS' and 'locus_tag' in f.qualifiers and 'pseudo' not in f.qualifiers:
                s = str(f.extract(rec.seq)).upper()
                if len(s) % 3 == 0 and len(s) >= 90:
                    cds[f.qualifiers['locus_tag'][0]] = (s, f.qualifiers.get('product', [''])[0])
    ref = [s for s, prod in cds.values() if prod.startswith(('50S ribosomal protein', '30S ribosomal protein'))]
    cnt = {}
    for s in ref:
        for i in range(0, len(s) - 3, 3):
            cnt[s[i:i + 3]] = cnt.get(s[i:i + 3], 0) + 1
    w = {}
    for aa, cs in syn.items():
        m = max(cnt.get(c, 0) for c in cs) or 1
        for c in cs:
            w[c] = max(cnt.get(c, 0), 0.5) / m
    rows = []
    for b, (s, _) in cds.items():
        cod = [s[i:i + 3] for i in range(0, len(s) - 3, 3)]
        lw = [np.log(w[c]) for c in cod if c in w and len(syn[fwd[c]]) > 1]
        gc3 = np.mean([c[2] in 'GC' for c in cod])
        rows.append({'bnumber': b, 'cai': float(np.exp(np.mean(lw))) if lw else 0.0, 'gc3': gc3})
    return pd.DataFrame(rows)
