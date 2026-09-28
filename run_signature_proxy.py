"""Locked same-patient, non-commercial signature proxy audit.
Preregistration: docs/PREREG_SIGNATURE_PROXY_20260928.md.
"""
import csv, hashlib, json, sys
from collections import defaultdict
from pathlib import Path
import numpy as np
import pandas as pd
import pyreadr
from scipy.stats import spearmanr
from sklearn.model_selection import train_test_split
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
sys.path.insert(0, 'src')
from recurscan.data.dataset import assemble

ROOT = Path(__file__).parent
MATRIX = ROOT/'data_cache/metabric_full_zscores.txt'
assert hashlib.sha256(MATRIX.read_bytes()).hexdigest() == 'f2ec4e1badc6f3a49c5db323e90faf3cb091cd29e880f58b3163bf24d30f6b7c'
TEMPLATE_SHA = {'sig.gene70.rda': '83362ba148219237aa282c81a6312bd7de04e55ded149fab5a6fee57e0e89f96',
                'sig.oncotypedx.rda': '1820b10e2b28d340cd7bc000f5b201e461d8318b5398f05f6523e8ee10f214d7'}
for n,h in TEMPLATE_SHA.items():
    assert hashlib.sha256((ROOT/'data/signature_proxy'/n).read_bytes()).hexdigest() == h

def signature(n):
    return list(pyreadr.read_r(str(ROOT/'data/signature_proxy'/n)).values())[0]
onco = signature('sig.oncotypedx.rda'); gene70 = signature('sig.gene70.rda')
assert len(onco) == 21 and len(gene70) == 70
oc = onco[onco['group']!='reference'].copy()
assert len(oc)==16 and len(set(oc['EntrezGene.ID']))==16
assert gene70['EntrezGene.ID'].isna().sum()==14
valid = gene70.dropna(subset=['EntrezGene.ID'])
centroid = valid.groupby(valid['EntrezGene.ID'].astype(int))['average.good.prognosis.profile'].mean()
assert len(centroid)==52 and len(valid)==56
needed = set(oc['EntrezGene.ID'].astype(int)) | set(centroid.index)

ds=assemble(); samples=set(ds.sample_ids)
rows=defaultdict(list); header=None
with MATRIX.open() as f:
    reader=csv.reader(f,delimiter='\t'); header=next(reader)
    assert len(header)==1982 and len(set(header[2:]))==1980
    sample_col = {sid: i for i,sid in enumerate(header) if sid in samples}
    assert len(sample_col)==len(ds.sample_ids)==1975
    sample_idxs=[sample_col[sid] for sid in ds.sample_ids]
    for r in reader:
        try: gid=int(r[1])
        except ValueError: continue
        if gid in needed:
            rows[gid].append([float(r[j]) if r[j] not in ('', 'NA') else np.nan for j in sample_idxs])
assert set(rows)==needed, ('missing genes',needed-set(rows))
expression={g: np.nanmean(np.array(arr,float),axis=0) for g,arr in rows.items()}
all_ids=sorted(needed)
valid_rows=np.array([all(np.isfinite(expression[g][i]) for g in all_ids) for i in range(len(ds.sample_ids))])
idx=np.flatnonzero(valid_rows)
assert len(idx)>0
print('eligibility',len(idx),'events',int(ds.event[idx].sum()),'missing',len(ds.sample_ids)-len(idx),flush=True)
tr,te=train_test_split(idx,test_size=.2,random_state=42,stratify=ds.event[idx]);
assert not set(tr)&set(te)
X=ds.X_expr.astype(float); T=ds.time.astype(float); E=ds.event.astype(int)
clinidx=None
frame=lambda ii: pd.DataFrame({**{f'x{j}': X[ii,j] for j in range(X.shape[1])}, 't':T[ii], 'e':E[ii]})
cox=CoxPHFitter(penalizer=.05).fit(frame(tr),'t','e')
coxrisk=cox.predict_partial_hazard(frame(te).drop(columns=['t','e'])).values.ravel()

# Cancer-gene mapping: the R code aliases CTSL2 -> CTSV by Entrez ID.
name_to_gid={s:int(g) for s,g in zip(oc['symbol'],oc['EntrezGene.ID'])}
assert 'CTSL2' in name_to_gid and name_to_gid['CTSL2'] in needed
scaled={}
for name,g in name_to_gid.items():
    lower=np.min(expression[g][tr]); upper=np.max(expression[g][tr]);
    assert upper>lower,(name,g)
    scaled[name]=np.clip((expression[g][te]-lower)/(upper-lower)*15,0,15)
grp=np.maximum(8,.9*scaled['GRB7']+.1*scaled['ERBB2'])
er=(.8*scaled['ESR1']+1.2*scaled['PGR']+scaled['BCL2']+scaled['SCUBE2'])/4
prolif=np.maximum(6.5,np.mean([scaled[g] for g in ['BIRC5','MKI67','MYBL2','CCNB1','AURKA']],axis=0))
inv=(scaled['CTSL2']+scaled['MMP11'])/2
oncorisk=.47*grp-.34*er+1.04*prolif+.10*inv+.05*scaled['CD68']-.08*scaled['GSTM1']-.07*scaled['BAG1']
ids=sorted(centroid.index); good=centroid.loc[ids].values
mamrisk=np.array([-spearmanr([expression[g][i] for g in ids],good).statistic for i in te])
assert np.isfinite(mamrisk).all() and np.isfinite(oncorisk).all()
sc={'cox70_expr':coxrisk,'oncotype_formula_proxy':oncorisk,'mamma_centroid_proxy':mamrisk}
ci=lambda ix,s:float(concordance_index(T[te][ix],-s[ix],E[te][ix]))
whole=np.arange(len(te)); points={n:ci(whole,s) for n,s in sc.items()}
rng=np.random.default_rng(20260928); draws={n:[] for n in sc}; diffs={n:[] for n in sc if n!='cox70_expr'}; skipped=0
for b in range(1000):
    ii=rng.integers(0,len(te),len(te))
    try:
        vals={n:ci(ii,s) for n,s in sc.items()}
    except ZeroDivisionError:
        skipped+=1;continue
    for n,v in vals.items():draws[n].append(v)
    for n in diffs:diffs[n].append(vals[n]-vals['cox70_expr'])
res={'protocol':'docs/PREREG_SIGNATURE_PROXY_20260928.md', 'dataset_sha256':'f2ec4e1badc6f3a49c5db323e90faf3cb091cd29e880f58b3163bf24d30f6b7c',
    'n_eligible':len(idx),'n_events_eligible':int(E[idx].sum()),'excluded_missing_signature':int(len(ds.sample_ids)-len(idx)),
    'n_train':len(tr),'n_holdout':len(te),'n_holdout_events':int(E[te].sum()), 'holdout_patient_ids':[ds.sample_ids[i] for i in te],
    'oncotype_cancer_genes':len(oc), 'oncotype_reference_genes_omitted':5, 'mamma_mapped_probes':len(valid), 'mamma_mapped_unique_ids':len(ids), 'mamma_unmapped_probes':14,
    'cindex':points,'cindex_ci95':{n:np.percentile(v,[2.5,97.5]).tolist() for n,v in draws.items()},
    'proxy_minus_cox_delta':{n:points[n]-points['cox70_expr'] for n in diffs},'proxy_minus_cox_ci95':{n:np.percentile(v,[2.5,97.5]).tolist() for n,v in diffs.items()},
    'bootstrap_B':1000,'bootstrap_seed':20260928,'bootstrap_skipped_degenerate':skipped,
    'limits':'Same-patient public formula/partial-centroid proxies on METABRIC z-scores. Neither is a marketed MammaPrint or Oncotype DX assay. Single random split, training-derived scaling, RFS not necessarily signature validation endpoint; no clinical or external superiority claim.'}
(ROOT/'results/signature_proxy.json').write_text(json.dumps(res,indent=1)+'\n')
print(json.dumps({k:v for k,v in res.items() if k!='holdout_patient_ids'},indent=1),flush=True)
