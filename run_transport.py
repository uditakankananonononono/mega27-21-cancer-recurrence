"""Round 9: transport METABRIC-trained Cox to GSE7390 vs published comparators.
Verified 10:36 PM: reproduces committed transport_round9.json point estimates exactly.
"""
import json, sys
import numpy as np, pandas as pd
sys.path.insert(0, "src")
import GEOparse
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from recurscan.data.dataset import assemble

ds = assemble()
genes_in = ds.gene_names
npi = ds.X_clin[:, [ds.clin_names.index("NPI")]].astype(float)
Xe = ds.X_expr.astype(float)
t, e = ds.time.astype(float), ds.event.astype(int)
def fit_cox(Xb, tag):
    df = {f"x{j}": Xb[:, j] for j in range(Xb.shape[1])}
    df["t"], df["e"] = t, e
    f = CoxPHFitter(penalizer=0.05).fit(pd.DataFrame(df), "t", "e")
    return f.summary["coef"].values
coef_full = fit_cox(np.hstack([npi, Xe]), "full")
coef_expr = fit_cox(Xe, "expr")

g = GEOparse.get_GEO(geo="GSE7390", destdir="data_cache/external/geo", annotate_gpl=False, silent=True)
gpl = GEOparse.get_GEO(geo="GPL96", destdir="data_cache/external/geo", silent=True)
p2g = dict(zip(gpl.table["ID"], gpl.table["Gene Symbol"]))
samples = list(g.gsms.values())
def chars(s):
    return dict(c.split(": ", 1) if ": " in c else (c.split(":")[0], "") for c in s.metadata.get("characteristics_ch1", []))
rows, keep = [], []
for s in samples:
    c = chars(s)
    try:
        rows.append({"t": float(c["t.rfs"]), "e": int(c["e.rfs"]),
                     "npi": float(c["NPI"]), "veridex": c["veridex_risk"], "ggi": c["risksg"]})
        keep.append(s)
    except (KeyError, ValueError):
        pass
d = pd.DataFrame(rows)
gids = np.array([p2g.get(i, "") for i in keep[0].table["ID_REF"].values])
X = np.column_stack([s.table["VALUE"].values for s in keep])
Xg = np.zeros((len(d), len(genes_in)))
for j, gene in enumerate(genes_in):
    m = gids == gene
    if m.sum():
        row = X[m].mean(0)
        Xg[:, j] = (row - row.mean()) / (row.std() + 1e-9)
npi_z = (d.npi.values - d.npi.values.mean()) / d.npi.values.std()
eta_full = coef_full[0] * npi_z + Xg @ coef_full[1:]
eta_expr = Xg @ coef_expr
def cat_risk(s):
    u = {v: i for i, v in enumerate(sorted(set(s)))}
    return np.array([u[v] for v in s], dtype=float)
tv, ev = d.t.values, d.e.values
from recurscan.eval import bootstrap_cindex_ci
def boot(sc):
    lo, hi = bootstrap_cindex_ci(tv, ev, sc, n_boot=500, seed=0)
    return [round(float(lo), 4), round(float(hi), 4)]
out = {"cohort": "GSE7390", "n": int(len(d)), "events": int(ev.sum()),
       "ours_npi_expr": {"ci": round(float(concordance_index(tv, -eta_full, ev)), 4), "boot95": boot(eta_full)},
       "ours_expr_only": {"ci": round(float(concordance_index(tv, -eta_expr, ev)), 4), "boot95": boot(eta_expr)},
       "veridex76": {"ci": round(float(concordance_index(tv, -cat_risk(d.veridex), ev)), 4), "boot95": boot(cat_risk(d.veridex))},
       "ggi": {"ci": round(float(concordance_index(tv, -cat_risk(d.ggi), ev)), 4), "boot95": boot(cat_risk(d.ggi))},
       "note": "ours trained on METABRIC, transported to Affymetrix cohort; comparators are GSE7390's own published risk calls (2-level). lifelines concordance_index: pass -risk (higher=longer survival)."}
json.dump(out, open("results/transport_round9.json", "w"), indent=2)
print(json.dumps(out, indent=1))
# direction sanity: mean eta by event
