"""Competing-risks descriptive check (cause-specific Cox, NOT Fine-Gray).

Design locked in this commit before outcomes: code three states from committed
METABRIC clinical records - recurrence (event 1), death without recurrence
(event 2), censored (0) - then fit cause-specific penalized Cox models for
each event type on the same clinical+expression features as the primary
benchmark. This answers only: does death-before-recurrence materially
compete, and do cause-specific fits shift the panel's behavior? Fine-Gray
subdistribution hazards remain NOT DONE (needs weighted partial likelihood
with time-varying censoring weights; no validated implementation here).
"""
import json,sys
import numpy as np,pandas as pd
from lifelines import CoxPHFitter
sys.path.insert(0,'src')
from recurscan.data.dataset import assemble

def code_competing(sample_ids,clinical):
    """Return arrays time, event (0 censored, 1 recurrence, 2 death w/o recurrence)."""
    t=[];ev=[];keep=[]
    for sid in sample_ids:
        rec=clinical.get(sid)
        if rec is None: keep.append(False);continue
        try:
            rfs_m=float(rec['RFS_MONTHS']); rfs=str(rec['RFS_STATUS'])
            os_m=float(rec['OS_MONTHS']); oss=str(rec['OS_STATUS'])
        except (TypeError,ValueError,KeyError):
            keep.append(False);continue
        keep.append(True)
        if rfs=='1:Recurred':
            t.append(rfs_m);ev.append(1)          # recurrence first
        elif oss=='1:DECEASED' and os_m<=rfs_m:
            t.append(os_m);ev.append(2)           # death without recurrence
        else:
            t.append(rfs_m);ev.append(0)          # censored
    return np.array(t),np.array(ev),np.array(keep,bool)

d=assemble()
clinical=json.load(open('data_cache/clinical.json'))
t,ev,keep=code_competing(list(d.sample_ids),clinical)
n_kept=int(keep.sum());n_drop=int((~keep).sum())
X=np.concatenate([d.X_clin,d.X_expr],axis=1)[keep]
names=list(d.clin_names)+list(d.gene_names)
t=t[keep if len(t)==len(keep) else slice(None)]; 
# t, ev already filtered inside coding for kept order? they were appended only when kept
counts={'kept':n_kept,'dropped_missing_OS_or_RFS':n_drop,
        'recurrence':int((ev==1).sum()),'death_without_recurrence':int((ev==2).sum()),
        'censored':int((ev==0).sum())}
f=pd.DataFrame({f'x{j}':X[:,j].astype(float) for j in range(X.shape[1])})
out={'design':__doc__,'counts':counts,'cause_specific':{}}
for label,code in (('recurrence',1),('death_without_recurrence',2)):
    f['t']=t;f['e']=(ev==code).astype(int)
    m=CoxPHFitter(penalizer=.05).fit(f,'t','e')
    conc=float(m.concordance_index_)
    top=m.params_.sort_values(key=np.abs,ascending=False).head(8)
    out['cause_specific'][label]={'concordance_apparent':conc,
        'top_abs_coef':{names[int(k[1:])]:float(v) for k,v in top.items()}}
    print(label,'C',round(conc,4),'events',int(f['e'].sum()),flush=True)
out['limits']='Cause-specific Cox censors the other event type; it is not the Fine-Gray subdistribution model. Apparent (in-sample) concordance; coding rules and counts are the deliverable.'
with open('results/competing_risks.json','w') as fh:json.dump(out,fh,indent=2)
print('counts',counts)
