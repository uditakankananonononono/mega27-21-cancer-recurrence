"""Pre-outcome 5-year absolute calibration audit. METABRIC 60-month RFS
baseline hazard transported without adjustment to independent cohorts.
Cohort endpoint mismatches (DMFS, metastasis, DRFS) explicitly flagged as
non-comparable to RFS; report their descriptive slope only, do not claim
absolute calibration. Requires run_transport_stability.load_cohort cache.
"""
import json, sys
import numpy as np,pandas as pd
from lifelines import CoxPHFitter,KaplanMeierFitter
from run_transport_stability import load_cohort,COHORTS
sys.path.insert(0,'src')
from recurscan.data.dataset import assemble

d=assemble();g=list(d.gene_names)
f=pd.DataFrame({f'x{j}':d.X_expr[:,j].astype(float) for j in range(len(g))})
f['t']=d.time.astype(float);f['e']=d.event.astype(int)
m=CoxPHFitter(penalizer=.05).fit(f,'t','e');beta=m.params_.values
S0=float(m.baseline_survival_.loc[m.baseline_survival_.index<=60].iloc[-1,0])
out={'metabric_baseline_survival_60m':S0,'cohorts':{},
     'warning':'Cohort-wise z-scoring and endpoint differences can invalidate transferred absolute risks; descriptive research audit only.'}
for name in COHORTS:
 c=load_cohort(name,g);eta=c['X']@beta
 # lifelines baseline survival is defined at training mean linear predictor.
 mean_eta=float(d.X_expr.mean(0)@beta)
 pred=1-S0**np.exp(np.clip(eta-mean_eta,-15,15))
 q=np.unique(np.quantile(pred,[0,.25,.5,.75,1]))
 bins=np.digitize(pred,q[1:-1],right=False)
 groups=[]
 for z in sorted(set(bins)):
  ix=bins==z
  if sum(ix)<10:continue
  km=KaplanMeierFitter().fit(c['t'][ix],c['e'][ix]);observed=1-float(km.predict(60))
  groups.append({'n':int(sum(ix)),'mean_predicted_5yr_event':float(pred[ix].mean()),
                 'km_observed_5yr_event':observed,'at_risk_5yr':int(sum(c['t'][ix]>=60))})
 # Cox slope: compare observed hazards to one unit of frozen risk score.
 slope=CoxPHFitter(penalizer=.05).fit(pd.DataFrame({'eta':eta,'t':c['t'],'e':c['e']}),'t','e')
 comparable=name in ('GSE7390','GSE2990')
 out['cohorts'][name]={'endpoint_RFS_comparable':comparable,'n':len(c['t']),
                      'slope':float(slope.params_['eta']),
                      'slope_ci95':[float(x) for x in slope.confidence_intervals_.loc['eta']],
                      'calibration_groups':groups,
                      'mean_predicted':float(pred.mean()),
                      'observed_KM_5yr_event':1-float(KaplanMeierFitter().fit(c['t'],c['e']).predict(60))}
 print(name,'comparable',comparable,'slope',round(out['cohorts'][name]['slope'],3),flush=True)
with open('results/transport_calibration.json','w') as fh:json.dump(out,fh,indent=2)
