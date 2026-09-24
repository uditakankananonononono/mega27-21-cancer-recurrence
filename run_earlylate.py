"""Round 5b: early (<=60mo) vs late (>60mo) recurrence - do predictors differ?
Early: Cox on events<=60mo (censor at 60). Late: landmark at 60mo on survivors.
Compare standardized Cox coefficients (expr panel + clinical) between windows."""
import json, sys
import numpy as np, pandas as pd
sys.path.insert(0, "src")
from lifelines import CoxPHFitter
from recurscan.data.dataset import assemble

ds = assemble()
t, e = ds.time.astype(float), ds.event.astype(int)
feat_names = ds.clin_names + ds.gene_names
X = np.hstack([ds.X_clin, ds.X_expr]).astype(float)

def fit_window(t_, e_, idx):
    df = {"t": t_[idx], "e": e_[idx]}
    for j in range(X.shape[1]):
        df[f"x{j}"] = X[idx, j]
    f = CoxPHFitter(penalizer=0.1).fit(pd.DataFrame(df), "t", "e")
    return f.summary["coef"], f.summary["p"]

CUT = 60.0
early_t = np.minimum(t, CUT); early_e = np.where(t <= CUT, e, 0)
alive = t > CUT  # landmark: survived recurrence-free past 60mo
late_t = t - CUT; late_e = e
ce, pe = fit_window(early_t, early_e, np.arange(len(t)))
cl, pl = fit_window(late_t, late_e, np.where(alive)[0])
def top(coef, p, k=10):
    idx = [i for i in np.argsort(p.values)[:k] if p.values[i] < 0.05]
    return [(feat_names[i], round(float(coef.values[i]), 3), float(p.values[i])) for i in idx]
both = np.corrcoef(ce.values, cl.values)[0, 1]
out = {"cut_months": CUT, "early_events": int(early_e.sum()),
       "late_events": int(late_e[np.where(alive)[0]].sum()),
       "late_n": int(alive.sum()),
       "coef_corr_early_vs_late": round(float(both), 4),
       "top_early": top(ce, pe), "top_late": top(cl, pl)}
json.dump(out, open("results/earlylate_round5.json", "w"), indent=2)
print(json.dumps(out, indent=1))
