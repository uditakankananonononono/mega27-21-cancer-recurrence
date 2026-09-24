"""External test of round-5c finding: recurrence predictability in an
ALL-CHEMOTHERAPY cohort (GSE25066, neoadjuvant taxane/anthracycline, n=508).
Within-cohort 5-seed held-out Cox C-index: clinical-only vs clinical+expr."""
import json, sys
import numpy as np, pandas as pd
sys.path.insert(0, "src")
import GEOparse
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index

g = GEOparse.get_GEO(geo="GSE25066", destdir="data_cache/external/geo", annotate_gpl=False, silent=True)
def chars(s):
    return dict(c.split(": ", 1) if ": " in c else (c.split(":")[0], "") for c in s.metadata.get("characteristics_ch1", []))
rows, keep_s = [], []
for s in g.gsms.values():
    c = chars(s)
    try:
        rows.append({"t": float(c["drfs_even_time_years"]) * 12,
                     "e": int(c["drfs_1_event_0_censored"]),
                     "age": float(c["age_years"]),
                     "grade": float(c["grade"]),
                     "node": 0.0 if c["clinical_nodal_status"].upper() == "N0" else 1.0,
                     "tstage": float(c["clinical_t_stage"].upper().replace("T", "")),
                     "er": 1.0 if c["er_status_ihc"].upper().startswith("P") else 0.0})
        keep_s.append(s)
    except (KeyError, ValueError):
        pass
df = pd.DataFrame(rows)
Xc = df[["age", "grade", "node", "tstage", "er"]].values
# expression of the 70-gene panel via GPL96
gpl = GEOparse.get_GEO(geo="GPL96", destdir="data_cache/external/geo", silent=True)
p2g = dict(zip(gpl.table["ID"], gpl.table["Gene Symbol"]))
from recurscan.data.dataset import assemble
panel = assemble().gene_names
gids = np.array([p2g.get(i, "") for i in keep_s[0].table["ID_REF"].values])
X = np.column_stack([s.table["VALUE"].values for s in keep_s])
Xg = []
for gene in panel:
    m = gids == gene
    if m.sum():
        row = X[m].mean(0)
        Xg.append((row - row.mean()) / (row.std() + 1e-9))
Xg = np.array(Xg).T
print("usable:", len(df), "events:", int(df.e.sum()), "panel genes mapped:", Xg.shape[1])
t, e = df.t.values, df.e.values
res = {"clin": [], "clin_expr": []}
for seed in range(5):
    rng = np.random.default_rng(100 + seed)
    perm = rng.permutation(len(df))
    te = perm[: len(df) // 5]; tr = perm[len(df) // 5:]
    for tag, Xf in [("clin", Xc), ("clin_expr", np.hstack([Xc, Xg]))]:
        keep = Xf[tr].std(0) > 1e-8
        d = {f"x{j}": Xf[tr][:, keep][:, j] for j in range(keep.sum())}
        d["t"], d["e"] = t[tr], e[tr]
        fit = CoxPHFitter(penalizer=0.1).fit(pd.DataFrame(d), "t", "e")
        risk = fit.predict_partial_hazard(pd.DataFrame(
            {f"x{j}": Xf[te][:, keep][:, j] for j in range(keep.sum())})).values.ravel()
        res[tag].append(concordance_index(t[te], -risk, e[te]))
out = {"cohort": "GSE25066", "arm": "all neoadjuvant chemotherapy", "n": int(len(df)),
       "events": int(df.e.sum()), "genes_mapped": int(Xg.shape[1]),
       "clin_ci_mean": round(float(np.mean(res["clin"])), 4), "clin_ci_std": round(float(np.std(res["clin"])), 4),
       "clin_expr_ci_mean": round(float(np.mean(res["clin_expr"])), 4),
       "clin_expr_ci_std": round(float(np.std(res["clin_expr"])), 4),
       "metabric_chemo_arm_ci": 0.5434}
json.dump(out, open("results/gse25066_chemo_round5c.json", "w"), indent=2)
print(json.dumps(out, indent=1))
