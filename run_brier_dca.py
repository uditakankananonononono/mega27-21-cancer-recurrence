"""Five-year Brier score and decision-curve audit (IPCW) for the frozen panel.

Design locked in this commit before outcomes: Brier at 60 months with
inverse-probability-of-censoring weights (Graf et al. form), null reference =
cohort KM 5-year event estimate for every patient; decision curves over
thresholds 0.05-0.60 with the same IPCW weights; treat-all and treat-none
reference lines. Cohorts: the two RFS-comparable external cohorts
(GSE7390, GSE2990). METABRIC is apparent (model trained there), reported
separately, never as validation.
"""
import json,sys
import numpy as np,pandas as pd
from lifelines import CoxPHFitter,KaplanMeierFitter
from run_transport_stability import load_cohort
sys.path.insert(0,'src')
from recurscan.data.dataset import assemble

H=60.0
THRESHOLDS=[round(0.05+0.05*i,2) for i in range(12)]

def ipcw_brier(t,e,p):
    t=np.asarray(t,float);e=np.asarray(e,int);p=np.asarray(p,float)
    n=len(t)
    km=KaplanMeierFitter().fit(t,1-e)  # censoring distribution
    G_t=km.predict(np.clip(t,0,H));G_H=float(km.predict(H))
    y=((t<=H)&(e==1)).astype(float)         # event by H
    atrisk=t>H                               # known event-free at H
    w=np.zeros(n)
    w[(t<=H)&(e==1)]=1.0/np.maximum(G_t[(t<=H)&(e==1)],1e-9)
    w[atrisk]=1.0/max(G_H,1e-9)
    brier=float(np.sum(w*(y-p)**2)/n)
    # null reference: KM 5y event prob for everyone
    p0=float(1-KaplanMeierFitter().fit(t,e).predict(H))
    brier_null=float(np.sum(w*(y-p0)**2)/n)
    n_known=int(((t<=H)&(e==1)).sum()+atrisk.sum())
    n_cens=int(((t<=H)&(e==0)).sum())
    return brier,brier_null,p0,n_known,n_cens,w,y

def decision_curves(t,e,p,w,y):
    n=len(t);rows=[]
    for pt in THRESHOLDS:
        treat=p>=pt
        tp=float(np.sum(w*treat*y));fp=float(np.sum(w*treat*(1-y)))
        nb=tp/n-fp/n*(pt/(1-pt))
        # treat-all line with same weights
        tpa=float(np.sum(w*y));fpa=float(np.sum(w*(1-y)))
        nb_all=tpa/n-fpa/n*(pt/(1-pt))
        rows.append({'pt':pt,'net_benefit':nb,'treat_all':nb_all,'treat_none':0.0})
    return rows

d=assemble();g=list(d.gene_names)
f=pd.DataFrame({f'x{j}':d.X_expr[:,j].astype(float) for j in range(len(g))})
f['t']=d.time.astype(float);f['e']=d.event.astype(int)
m=CoxPHFitter(penalizer=.05).fit(f,'t','e');beta=m.params_.values
S0=float(m.baseline_survival_.loc[m.baseline_survival_.index<=H].iloc[-1,0])
mean_eta=float(d.X_expr.mean(0)@beta)

def cohort_pack(name,t,e,X):
    eta=X@beta
    p=1-S0**np.exp(np.clip(eta-mean_eta,-15,15))
    brier,brier_null,p0,n_known,n_cens,w,y=ipcw_brier(t,e,p)
    return {'n':int(len(t)),'brier_ipcw_60m':brier,'brier_null_60m':brier_null,
            'km_5yr_event':p0,'n_known_status_60m':n_known,'n_censored_before_60m':n_cens,
            'brier_skill_vs_null':1-brier/brier_null,
            'decision_curve':decision_curves(t,e,p,w,y)}

out={'design':__doc__,'horizon_months':H,'thresholds':THRESHOLDS,'cohorts':{}}
out['cohorts']['METABRIC_apparent']=dict(cohort_pack('METABRIC',d.time,d.event,d.X_expr),role='apparent (training cohort), not validation')
for name in ('GSE7390','GSE2990'):
    c=load_cohort(name,g)
    out['cohorts'][name]=dict(cohort_pack(name,c['t'],c['e'],c['X']),role='external, RFS-comparable')
    print(name,'brier',round(out['cohorts'][name]['brier_ipcw_60m'],4),'null',round(out['cohorts'][name]['brier_null_60m'],4),flush=True)
out['limits']='Baseline survival transported without cohort refit; IPCW assumes independent censoring given time; descriptive audit, not clinical validation.'
with open('results/brier_dca.json','w') as fh:json.dump(out,fh,indent=2)
print('wrote results/brier_dca.json')
