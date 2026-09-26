"""Pre-outcome secondary biological coherence audit, judge round 3.

Do not alter the primary five-cohort transport-stability verdict. This test
asks if a stable gene subset repeatedly maps to a Hallmark program, beyond
what is expected from this deliberately biology-enriched 70-gene panel.
The entire 70-gene panel is the enrichment universe, NOT all human genes.
MSigDB H 2025.1 Hs GMT is fetched live, not redistributed in the repo.
"""
import hashlib, json, urllib.request
from collections import Counter
from pathlib import Path
import numpy as np
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

ROOT=Path(__file__).resolve().parent
URL='https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2025.1.Hs/h.all.v2025.1.Hs.symbols.gmt'
SEED=3; B=1000; ALPHA=.05; MIN_OVERLAP=3; MIN_PROGRAM_PANEL=3

def read_programs(raw, universe):
    out={}
    for line in raw.decode().splitlines():
        cols=line.split('\t'); name=cols[0]
        g=set(cols[2:]) & universe
        if MIN_PROGRAM_PANEL <= len(g) < len(universe):out[name]=g
    assert len(out)>0, 'No Hallmark set intersects panel by at least 3 genes'
    return out

def fold_q(panel, programs):
    names=list(programs)
    p=[float(hypergeom.sf(len(panel & programs[n])-1, 70, len(programs[n]), len(panel)))
       if len(panel & programs[n])>=MIN_OVERLAP else 1. for n in names]
    q=multipletests(p,method='fdr_bh')[1]
    return dict(zip(names,map(float,q)))

def score(folds, programs, genes):
    qs=[fold_q(s,programs) for s in folds]
    counts={n:sum(q[n]<ALPHA for q in qs) for n in programs}
    return max(counts.values()), counts, qs

def analyze(j, raw):
    genes=sorted(j['gene_map'])
    assert len(genes)==70, 'Panel length changed, halt instead of silently changing null'
    universe=set(genes); programs=read_programs(raw,universe)
    folds=[set(f['stable_genes']) for f in j['folds']]
    assert len(folds)==5 and all(s<=universe for s in folds)
    if not all(folds):return {'gate_pass':False,'reason':'zero stable genes in at least one fold; no coherence claim'}
    beta=np.array([abs(j['gene_map'][g]['metabric_coef']) for g in genes])
    # A *single* gene identity permutation is shared across all five folds:
    # preserves overlap structure between folds and (via quartiles) beta scale.
    order=np.argsort(beta,kind='stable'); bins=np.zeros(70,dtype=int)
    for i,block in enumerate(np.array_split(order,4)):bins[block]=i
    obs,ct,qs=score(folds,programs,genes)
    rng=np.random.default_rng(SEED); null=[]
    for _ in range(B):
        idx=np.arange(70)
        for k in range(4):
            ix=np.flatnonzero(bins==k);idx[ix]=rng.permutation(ix)
        mapping=dict(zip(genes,[genes[i] for i in idx]))
        perm_folds=[{mapping[g] for g in s} for s in folds]
        null.append(score(perm_folds,programs,genes)[0])
    # Corrects selection over *all* Hallmark programs by max-statistic.
    p=(1+sum(x>=obs for x in null))/(B+1)
    repeated={n:ct[n] for n in programs if ct[n]>=4}
    return {'stat_max_significant_folds':obs,'repeated_programs_4of5':repeated,
            'fold_stable_sizes':list(map(len,folds)),'programs_tested':len(programs),
            'panel_background_size':70,'q_values_by_fold':qs,
            'null_max_fold_counts':null,'empirical_p':p,
            'gate_pass':bool(obs>=4 and p<.01),
            'interpretation':'Exploratory observational pathway coherence only; not proof of mechanism or equal endpoints'}

def main():
    j=json.loads((ROOT/'results/transport_stability.json').read_text())
    with urllib.request.urlopen(URL,timeout=30) as r:raw=r.read()
    z=analyze(j,raw)
    z.update({'source_url':URL,'source_sha256':hashlib.sha256(raw).hexdigest(),
              'locked_design':'MSigDB Hallmark v2025.1.Hs; 70-panel universe, >=3 panel genes per set and overlap >=3; hypergeometric upper tail, BH q<.05 per fold; max Hallmark repeat across folds >=4/5 and empirical p<.01, 1000 shared within-beta-quartile identity permutations, seed=3. Secondary only.'})
    (ROOT/'results/program_coherence.json').write_text(json.dumps(z,indent=2))
    print('coherence',z['gate_pass'],z.get('empirical_p'),flush=True)
if __name__=='__main__':main()
