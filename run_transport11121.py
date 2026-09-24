"""Round 9d: transport METABRIC expr-only Cox to GSE11121 (Schmidt 2008, Mainz,
node-negative), DMFS endpoint (t.dmfs in months). Comparator: tumor grade."""
import json, sys
import numpy as np, pandas as pd
sys.path.insert(0, "src")
import GEOparse
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from recurscan.data.dataset import assemble
from recurscan.eval import bootstrap_cindex_ci

ds = assemble()
genes_in = ds.gene_names
Xe = ds.X_expr.astype(float)
t, e = ds.time.astype(float), ds.event.astype(int)
df = {f"x{j}": Xe[:, j] for j in range(Xe.shape[1])}
df["t"], df["e"] = t, e
fit = CoxPHFitter(penalizer=0.05).fit(pd.DataFrame(df), "t", "e")
coef = fit.summary["coef"].values

g = GEOparse.get_GEO(geo="GSE11121", destdir="data_cache/external/geo", annotate_gpl=False, silent=True)
gpl = GEOparse.get_GEO(geo="GPL96", destdir="data_cache/external/geo", silent=True)
p2g = dict(zip(gpl.table["ID"], gpl.table["Gene Symbol"]))
def chars(s):
    return dict(c.split(": ", 1) if ": " in c else (c.split(":")[0], "") for c in s.metadata.get("characteristics_ch1", []))
rows, keep = [], []
for s in g.gsms.values():
    c = chars(s)
    try:
        rows.append({"t": float(c["t.dmfs"]), "e": int(c["e.dmfs"]), "grade": float(c["grade"])})
        keep.append(s)
    except (KeyError, ValueError):
        pass
d = pd.DataFrame(rows)
gids = np.array([p2g.get(i, "") for i in keep[0].table["ID_REF"].values])
X = np.column_stack([s.table["VALUE"].values for s in keep])
Xg = np.zeros((len(d), len(genes_in)))
mapped = 0
for j, gene in enumerate(genes_in):
    m = gids == gene
    if m.sum():
        row = X[m].mean(0)
        Xg[:, j] = (row - row.mean()) / (row.std() + 1e-9)
        mapped += 1
eta = Xg @ coef
tv, ev = d.t.values, d.e.values
def boot(sc):
    lo, hi = bootstrap_cindex_ci(tv, ev, sc, n_boot=500, seed=0)
    return [round(float(lo), 4), round(float(hi), 4)]
out = {"cohort": "GSE11121", "endpoint": "DMFS (months)", "n": int(len(d)), "events": int(ev.sum()),
       "genes_mapped": mapped,
       "ours_expr": {"ci": round(float(concordance_index(tv, -eta, ev)), 4), "boot95": boot(eta)},
       "grade": {"ci": round(float(concordance_index(tv, -d.grade.values, ev)), 4), "boot95": boot(d.grade.values)},
       "note": "ours = expr-only Cox (70-gene panel) trained on METABRIC, transported raw; node-negative cohort; comparator = tumor grade"}
json.dump(out, open("results/transport11121_round9d.json", "w"), indent=2)
print(json.dumps(out, indent=1))
