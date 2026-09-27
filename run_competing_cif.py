"""Competing-risk cumulative incidence (not Fine-Gray), preregistered in docs/PREREG_COMPETING_CIF_20260928.md."""
import json
import numpy as np
from lifelines import AalenJohansenFitter,KaplanMeierFitter
from run_competing_risks import code_competing
from recurscan.data.dataset import assemble
D=assemble();clinical=json.load(open('data_cache/clinical.json'))
t,e,keep=code_competing(D.sample_ids,clinical)
assert len(t)==len(D.time)==1975 and keep.all()
assert (int((e==1).sum()),int((e==2).sum()),int((e==0).sum()))==(800,431,744)
h=60
out={'protocol':'docs/PREREG_COMPETING_CIF_20260928.md','horizon_months':h,
     'counts':{'recurrence':800,'death_without_recurrence':431,'censored':744},
     'events_by_60m':{'recurrence':int(((e==1)&(t<=h)).sum()),'death_without_recurrence':int(((e==2)&(t<=h)).sum())},
     'cumulative_incidence_60m':{}}
for code,name in ((1,'recurrence'),(2,'death_without_recurrence')):
    aj=AalenJohansenFitter().fit(t,e,event_of_interest=code)
    idx=aj.cumulative_density_.index[aj.cumulative_density_.index<=h][-1]
    value=float(aj.cumulative_density_.loc[idx].iloc[0]);var=float(aj.variance_.loc[idx])
    out['cumulative_incidence_60m'][name]={'estimate':value,'variance':var,'normal_ci95_clipped':[max(0,value-1.96*np.sqrt(var)),min(1,value+1.96*np.sqrt(var))]}
km=KaplanMeierFitter().fit(t,(e==1).astype(int))
out['naive_recurrence_km_60m_treat_competing_death_as_censored']=float(1-km.predict(h))
out['difference_naive_minus_aj']=out['naive_recurrence_km_60m_treat_competing_death_as_censored']-out['cumulative_incidence_60m']['recurrence']['estimate']
out['limits']='nonparametric apparent cohort CIF, no regression or Fine-Gray; normal CI approximations clipped to [0,1]; competing event coding uses RFS/OS fields'
open('results/competing_cif.json','w').write(json.dumps(out,indent=1)+'\n')
print(json.dumps(out,indent=1))
