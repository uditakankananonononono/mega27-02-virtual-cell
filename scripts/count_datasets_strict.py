"""Strict dataset count: one per source study (parent ruling 2026-09-25 04:13). Writes results/datasets_strict_count.json."""
import json, re
m = json.load(open('results/tools_manifest.json')); ds = m['datasets']
def key(i, name):
    n = name
    if i in (2, 5): return 'model:iML1515'
    if i in (6, 7): return 'ncbi:U00096.3'
    if i in (8, 9): return 'gerdes2003'
    if i in (11, 12, 13): return 'price2018:Keio'
    if i in (14, 15, 16): return 'string:511145'
    if i in (18, 19): return 'kegg:eco'
    if i in (17, 44, 47, 50): return 'uniprot:UP000000625'
    if i == 20: return None                       # integrated PaxDb: derived from the per-dataset files
    if i == 41: return None                       # FB aaseqs: metadata of the same figshare deposit
    if i in (24, 66): return 'price2018:MR1'
    if 'PXD014877_Mueller' in n: return 'mueller2020:PXD014877'
    if 'MSV000096603' in n: return 'MSV000096603'
    if 'PXD009705' in n: return 'PXD009705'
    if 'arike_2012' in n: return 'arike2012'
    if 'Krug_2013' in n: return 'krug2013'
    if 'ecoli1_HCD' in n or 'ecoli2_resolution' in n: return 'pride:ecoli1-2'
    if 'zhang_2006' in n: return 'zhang2006:DvH'
    if 'Albrethsen' in n: return 'albrethsen2013:PXD000111'
    if 'PXD013677' in n: return 'PXD013677'
    if 'PA_2013-7' in n or 'PA_201307' in n: return 'peptideatlas:Mtb'
    if 'atlas_build_330' in n: return 'peptideatlas:Mtb'
    return f'entry:{i}'
keys = {}
for i, e in enumerate(ds):
    k = key(i, e[0])
    if k: keys.setdefault(k, []).append(i)
out = dict(rule='one per source study; re-deposits, quantification variants and fractions of one study collapse; derived/integrated tables and deposit metadata not counted',
           n_manifest_entries=len(ds), n_strict=len(keys), collapsed={k: v for k, v in keys.items() if len(v) > 1},
           dropped=[20, 41])
json.dump(out, open('results/datasets_strict_count.json', 'w'), indent=1)
print(out['n_manifest_entries'], out['n_strict']); print(out['collapsed'])
