"""Round 9b: transport benchmark replication on GSE25066 (all-chemo cohort).
Comparators from cohort metadata: ggi_class, dlda30, chemosensitivity_prediction."""
import json, sys
import numpy as np, pandas as pd
sys.path.insert(0, "src")
import GEOparse
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from recurscan.data.dataset import assemble

ds = assemble()
t, e = ds.time.astype(float), ds.event.astype(int)
df = {f"x{j}": ds.X_expr[:, j].astype(float) for j in range(ds.X_expr.shape[1])}
df["t"], df["e"] = t, e
coef = CoxPHFitter(penalizer=0.05).fit(pd.DataFrame(df), "t", "e").summary["coef"].values

g = GEOparse.get_GEO(geo="GSE25066", destdir="data_cache/external/geo", annotate_gpl=False, silent=True)
gpl = GEOparse.get_GEO(geo="GPL96", destdir="data_cache/external/geo", silent=True)
p2g = dict(zip(gpl.table["ID"], gpl.table["Gene Symbol"]))
def chars(s):
    return dict(c.split(": ", 1) if ": " in c else (c.split(":")[0], "") for c in s.metadata.get("characteristics_ch1", []))
rows, keep = [], []
for s in g.gsms.values():
    c = chars(s)
    try:
        rows.append({"t": float(c["drfs_even_time_years"]) * 12,
                     "e": int(c["drfs_1_event_0_censored"]),
                     "ggi": c["ggi_class"], "dlda": c["dlda30_prediction"],
                     "sens": c["chemosensitivity_prediction"]})
        keep.append(s)
    except (KeyError, ValueError):
        pass
d = pd.DataFrame(rows)
gids = np.array([p2g.get(i, "") for i in keep[0].table["ID_REF"].values])
X = np.column_stack([s.table["VALUE"].values for s in keep])
Xg = np.zeros((len(d), len(ds.gene_names)))
mapped = 0
for j, gene in enumerate(ds.gene_names):
    m = gids == gene
    if m.sum():
        row = X[m].mean(0)
        Xg[:, j] = (row - row.mean()) / (row.std() + 1e-9)
        mapped += 1
eta = Xg @ coef
tt, ee = d.t.values, d.e.values
def cat_risk(s):
    u = {v: i for i, v in enumerate(sorted(set(s)))}
    return np.array([u[v] for v in s], dtype=float)
rng = np.random.default_rng(0)
def boot(score, n=200):
    cis = []
    for _ in range(n):
        i = rng.integers(0, len(tt), len(tt))
        if ee[i].sum() < 5: continue
        cis.append(concordance_index(tt[i], -score[i], ee[i]))
    return [round(float(np.percentile(cis, 2.5)), 4), round(float(np.percentile(cis, 97.5)), 4)]
def ev(score, cat=False):
    sc = cat_risk(score) if cat else score
    return {"ci": round(float(concordance_index(tt, -sc, ee)), 4), "boot95": boot(sc)}
out = {"cohort": "GSE25066", "arm": "all neoadjuvant chemo", "n": int(len(d)),
       "events": int(ee.sum()), "genes_mapped": mapped,
       "ours_expr_only": ev(eta), "ggi_class": ev(d.ggi, True),
       "dlda30": ev(d.dlda, True), "chemosensitivity": ev(d.sens, True),
       "note": "ours = Cox(70-gene expr) trained on METABRIC; comparators are cohort's published calls"}
json.dump(out, open("results/transport25066_round9b.json", "w"), indent=2)
print(json.dumps(out, indent=1))
