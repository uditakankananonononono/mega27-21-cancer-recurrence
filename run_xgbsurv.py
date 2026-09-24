"""Round 8: XGBoost AFT survival vs linear Cox leader (0.676 same-split base).
Nonlinear model class on clin+expr. Same 5 seeds/splits."""
import json, sys
import numpy as np
sys.path.insert(0, "src")
import xgboost as xgb
from lifelines.utils import concordance_index
from scipy import stats
from recurscan.data.dataset import assemble

ds = assemble()
t, e = ds.time.astype(float), ds.event.astype(int)
X = np.hstack([ds.X_clin, ds.X_expr]).astype(float)

res = {"xgb": []}
for seed in range(5):
    rng = np.random.default_rng(20260924 + seed)
    perm = rng.permutation(len(t))
    te = perm[: len(t) // 5]; tr = perm[len(t) // 5:]
    dtr = xgb.DMatrix(X[tr])
    # AFT: lower/upper bounds; events exact, censored right-open
    yl = t[tr].copy(); yu = np.where(e[tr] == 1, t[tr], np.inf)
    dtr.set_float_info("label_lower_bound", yl)
    dtr.set_float_info("label_upper_bound", yu)
    m = xgb.train({"objective": "survival:aft", "eval_metric": "aft-nloglik",
                   "aft_loss_distribution": "logistic", "max_depth": 3,
                   "eta": 0.05, "lambda": 1.0, "seed": seed},
                  dtr, num_boost_round=200)
    pred_t = m.predict(xgb.DMatrix(X[te]))
    risk = -pred_t  # longer predicted survival = lower risk
    ci = concordance_index(t[te], -risk, e[te])
    res["xgb"].append(ci)
    print(f"seed {seed}: xgb-aft {ci:.4f}", flush=True)
base = json.load(open("results/multimodal_round6.json"))["base_runs"]
out = {"xgb_aft_runs": [round(float(x), 4) for x in res["xgb"]],
       "xgb_aft_mean": round(float(np.mean(res["xgb"])), 4),
       "base_mean_same_splits": 0.676,
       "delta_vs_base": round(float(np.mean(res["xgb"])) - 0.676, 4),
       "paired_t_p_vs_base": float(stats.ttest_rel(res["xgb"], base).pvalue)}
json.dump(out, open("results/xgbsurv_round8.json", "w"), indent=2)
print(json.dumps(out, indent=1))
