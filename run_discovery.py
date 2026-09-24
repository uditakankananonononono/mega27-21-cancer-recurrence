import json, sys
sys.path.insert(0, "src")
import numpy as np
from recurscan.data.dataset import assemble
from recurscan.benchmark import standardize, stratified_event_split
from recurscan.discovery import (bootstrap_delta_ci, cox_risk, fit_cox,
                                 graph_smoothed_risk, module_scores)
from recurscan.survival import concordance_index_fast

ds = assemble()
mods = module_scores(ds)
M = np.column_stack(list(mods.values()))
X_base = np.concatenate([ds.X_clin, ds.X_expr], axis=1)

results = {}
for seed in (0, 1, 2):
    tr, te = stratified_event_split(ds.event, 0.25, seed)
    Xtr, Xte = standardize(X_base[tr], X_base[te])
    Mtr, Mte = standardize(M[tr], M[te])
    base = fit_cox(Xtr, ds.time[tr], ds.event[tr])
    risk_b_te = cox_risk(base, Xte)
    # candidate A: baseline features + module scores
    candA = fit_cox(np.hstack([Xtr, Mtr]), ds.time[tr], ds.event[tr])
    risk_a_te = cox_risk(candA, np.hstack([Xte, Mte]))
    # candidate B: graph-smoothed baseline risk (transductive)
    risk_b_tr = cox_risk(base, Xtr)
    Xall = np.vstack([Xtr, Xte])
    smooth = graph_smoothed_risk(Xall, risk_b_tr, np.arange(len(Xtr)))
    risk_s_te = 0.7 * risk_b_te + 0.3 * smooth[len(Xtr):]
    cb = concordance_index_fast(ds.time[te], ds.event[te], risk_b_te)
    ca = concordance_index_fast(ds.time[te], ds.event[te], risk_a_te)
    cs = concordance_index_fast(ds.time[te], ds.event[te], risk_s_te)
    d_a = bootstrap_delta_ci(ds.time[te], ds.event[te], risk_a_te, risk_b_te, seed=seed)
    d_s = bootstrap_delta_ci(ds.time[te], ds.event[te], risk_s_te, risk_b_te, seed=seed)
    results[seed] = {"base": cb, "candA_modules": ca, "candB_smoothed": cs,
                     "delta_A": d_a, "delta_B": d_s}
    print(f"seed {seed}: base {cb:.4f} | +modules {ca:.4f} (d {d_a[0]:+.4f} CI{d_a[1]}) | +graph {cs:.4f} (d {d_s[0]:+.4f} CI{d_s[1]})", flush=True)
json.dump(results, open("results/discovery_round1.json", "w"), indent=2, default=float)
print("SAVED", flush=True)
