"""Round 5c: treatment-arm-specific recurrence dynamics (METABRIC).
Per-arm penalized Cox (clin+expr) in chemo / hormone / neither arms:
per-arm C-index + top coefficients. Falsifiable question: does the risk
gradient differ by systemic treatment arm?"""
import json, sys
import numpy as np, pandas as pd
sys.path.insert(0, "src")
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from recurscan.data.dataset import assemble

ds = assemble()
names = ds.clin_names + ds.gene_names
X = np.hstack([ds.X_clin, ds.X_expr]).astype(float)
t, e = ds.time.astype(float), ds.event.astype(int)
chemo = ds.X_clin[:, ds.clin_names.index("CHEMOTHERAPY=YES")].astype(bool)
horm = ds.X_clin[:, ds.clin_names.index("HORMONE_THERAPY=YES")].astype(bool)
arms = {"chemo": chemo, "hormone": horm & ~chemo, "neither": ~chemo & ~horm}
out = {}
for arm, mask in arms.items():
    idx = np.where(mask)[0]
    if mask.sum() < 150 or e[idx].sum() < 60:
        out[arm] = {"n": int(mask.sum()), "events": int(e[idx].sum()), "skipped": "too small"}
        continue
    keep = X[idx].std(0) > 1e-8
    Xk = X[:, keep]
    names_k = [n for n, k in zip(names, keep) if k]
    cis, tops = [], []
    for seed in range(5):
        rng = np.random.default_rng(7 + seed)
        perm = rng.permutation(len(idx))
        te = idx[perm[: len(idx) // 5]]; tr = idx[perm[len(idx) // 5:]]
        df = {f"x{j}": Xk[tr, j] for j in range(Xk.shape[1])}
        df["t"], df["e"] = t[tr], e[tr]
        fit = CoxPHFitter(penalizer=0.1).fit(pd.DataFrame(df), "t", "e")
        risk = fit.predict_partial_hazard(pd.DataFrame(
            {f"x{j}": Xk[te, j] for j in range(Xk.shape[1])})).values.ravel()
        cis.append(concordance_index(t[te], -risk, e[te]))
        if seed == 0:
            s_ = fit.summary
            tops = [(names_k[i], round(float(s_.iloc[i]["coef"]), 3), round(float(s_.iloc[i]["p"]), 4))
                    for i in np.argsort(s_["p"].values)[:6] if s_.iloc[i]["p"] < 0.05]
    out[arm] = {"n": int(mask.sum()), "events": int(e[idx].sum()),
                "heldout_cindex_mean": round(float(np.mean(cis)), 4),
                "heldout_cindex_std": round(float(np.std(cis)), 4),
                "heldout_cindex_runs": [round(float(c), 4) for c in cis],
                "top_predictors_seed0": tops}
json.dump(out, open("results/treatment_arms_round5c.json", "w"), indent=2)
print(json.dumps(out, indent=1))
