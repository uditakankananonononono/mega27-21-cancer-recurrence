"""Round 7: engineered features vs clin+expr leader.
+ proliferation index (mean z of mitotic genes), NPI^2, NPI x ER interaction,
grade x proliferation. Same 5-seed protocol, leakage-free (all transforms fixed)."""
import json, sys
import numpy as np, pandas as pd
sys.path.insert(0, "src")
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from scipy import stats
from recurscan.data.dataset import assemble

ds = assemble()
t, e = ds.time.astype(float), ds.event.astype(int)
Xb = np.hstack([ds.X_clin, ds.X_expr]).astype(float)
MIT = ["CEP55", "CDC20", "BIRC5", "KIF2C", "ANLN", "MELK", "UBE2T", "PTTG1", "NDC80", "NUF2"]
gidx = [ds.gene_names.index(g) for g in MIT if g in ds.gene_names]
prolif = ds.X_expr[:, gidx].mean(1, keepdims=True)
npi = ds.X_clin[:, ds.clin_names.index("NPI")].reshape(-1, 1)
er = ds.X_clin[:, ds.clin_names.index("ER_STATUS=Positive")].reshape(-1, 1)
grade = ds.X_clin[:, ds.clin_names.index("GRADE")].reshape(-1, 1)
extra = np.hstack([prolif, npi ** 2, npi * er, grade * prolif])
extra = (extra - extra.mean(0)) / (extra.std(0) + 1e-9)
Xe = np.hstack([Xb, extra])

def run(Xtr, Xte, tr, te):
    d = {f"x{j}": Xtr[:, j] for j in range(Xtr.shape[1])}
    d["t"], d["e"] = t[tr], e[tr]
    fit = CoxPHFitter(penalizer=0.05).fit(pd.DataFrame(d), "t", "e")
    risk = fit.predict_partial_hazard(pd.DataFrame(
        {f"x{j}": Xte[:, j] for j in range(Xte.shape[1])})).values.ravel()
    return concordance_index(t[te], -risk, e[te])

res = {"base": [], "eng": []}
for seed in range(5):
    rng = np.random.default_rng(20260924 + seed)
    perm = rng.permutation(len(t))
    te = perm[: len(t) // 5]; tr = perm[len(t) // 5:]
    res["base"].append(run(Xb[tr], Xb[te], tr, te))
    res["eng"].append(run(Xe[tr], Xe[te], tr, te))
    print(f"seed {seed}: base {res['base'][-1]:.4f} eng {res['eng'][-1]:.4f}", flush=True)
out = {"base_mean": round(float(np.mean(res["base"])), 4),
       "eng_mean": round(float(np.mean(res["eng"])), 4),
       "delta": round(float(np.mean(res["eng"]) - np.mean(res["base"])), 4),
       "paired_t_p": float(stats.ttest_rel(res["eng"], res["base"]).pvalue),
       "features_added": ["prolif_index", "NPI^2", "NPIxER", "GRADExprolif"]}
json.dump(out, open("results/features_round7.json", "w"), indent=2)
print(json.dumps(out, indent=1))
