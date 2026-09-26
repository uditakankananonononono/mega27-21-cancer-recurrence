"""Pre-outcome judge-round-1 analysis: does the 70-gene panel add information
beyond grade, ER, HER2, nodal burden and a seven-gene proliferation subset?
The nested, unpenalized Cox partial-likelihood ratio is descriptive (many
correlated genes). A fixed 80/20 held-out split provides a stronger
out-of-sample log partial-likelihood and C-index check. No clinical use claim.
"""
import json, sys
import numpy as np, pandas as pd
from scipy.stats import chi2
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from sklearn.model_selection import train_test_split
sys.path.insert(0, 'src')
from recurscan.data.dataset import assemble

PROLIF = ['MKI67','BIRC5','CCNB1','CDC20','CENPF','MYBL2','UBE2C']
CLIN = ['GRADE','ER_STATUS=Positive','HER2_STATUS=Positive',
        'LYMPH_NODES_EXAMINED_POSITIVE']
d = assemble()
base_idx = [d.clin_names.index(x) for x in CLIN]
pro_idx = [d.gene_names.index(x) for x in PROLIF]
X = np.column_stack([d.X_clin[:, base_idx], d.X_expr[:, pro_idx]])
Z = np.column_stack([d.X_clin[:, base_idx], d.X_expr])
assert X.shape[1] == 11 and Z.shape[1] == 74 and len(d.time) > 1900
tr, te = train_test_split(np.arange(len(d.time)), test_size=.2, random_state=42,
                           stratify=d.event)

def frame(arr, idx):
    f = pd.DataFrame(arr[idx], columns=[f'x{i}' for i in range(arr.shape[1])])
    f['t'] = d.time[idx].astype(float)
    f['e'] = d.event[idx].astype(int)
    return f

out = {'predeclared_proliferation_genes': PROLIF, 'clinical_covariates': CLIN,
       'n': int(len(d.time)), 'events': int(sum(d.event)),
       'n_train': len(tr), 'n_test': len(te), 'split_seed':42}
fits = {}
for name, arr in [('base',X), ('full',Z)]:
    model = CoxPHFitter(penalizer=.05)
    model.fit(frame(arr,tr),'t','e')
    fits[name] = model
    risk = model.predict_partial_hazard(frame(arr,te).drop(columns=['t','e'])).values.reshape(-1)
    out[name] = {'heldout_cindex':float(concordance_index(d.time[te],-risk,d.event[te])),
                 'train_log_likelihood':float(model.log_likelihood_)}
# A chi-square LR reference is only valid for unpenalized nested models.
# Try unpenalized models separately; nonconvergence yields a clearly marked NA.
try:
    u0 = CoxPHFitter().fit(frame(X,tr),'t','e')
    u1 = CoxPHFitter().fit(frame(Z,tr),'t','e')
    stat = max(0,2*(u1.log_likelihood_ - u0.log_likelihood_))
    out['unpenalized_LR'] = {'stat':float(stat),'df':int(Z.shape[1]-X.shape[1]),
                             'p':float(chi2.sf(stat,Z.shape[1]-X.shape[1])),
                             'caveat':'training-set nested-model test; not external validation'}
except Exception as ex:
    out['unpenalized_LR'] = {'status':'not estimable','error':str(ex)[:300]}
with open('results/beyond_proliferation.json','w') as fh: json.dump(out,fh,indent=2)
print(json.dumps(out,indent=2))
