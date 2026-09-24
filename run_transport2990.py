"""Round 9c: transport METABRIC expr-only Cox to GSE2990 (Sotiriou 2006),
compare vs the per-sample continuous GGI published with the series.
Clinical outcomes from GSE2990_suppl_info.txt (event.rfs/time.rfs in years)."""
import json, sys
import numpy as np, pandas as pd
sys.path.insert(0, "src")
import GEOparse
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from recurscan.data.dataset import assemble

ds = assemble()
genes_in = ds.gene_names
Xe = ds.X_expr.astype(float)
t, e = ds.time.astype(float), ds.event.astype(int)
df = {f"x{j}": Xe[:, j] for j in range(Xe.shape[1])}
df["t"], df["e"] = t, e
fit = CoxPHFitter(penalizer=0.05).fit(pd.DataFrame(df), "t", "e")
coef = fit.summary["coef"].values  # 70 genes, expr-only

g = GEOparse.get_GEO(geo="GSE2990", destdir="data_cache/external/geo", annotate_gpl=False, silent=True)
gpl = GEOparse.get_GEO(geo="GPL96", destdir="data_cache/external/geo", silent=True)
p2g = dict(zip(gpl.table["ID"], gpl.table["Gene Symbol"]))
clin = pd.read_csv("data_cache/external/geo/GSE2990_suppl_info.txt", sep="\t")
clin = clin.dropna(subset=["event.rfs", "time.rfs"])
def chars(s):
    return dict(c.split(": ", 1) if ": " in c else (c.split(":")[0], "") for c in s.metadata.get("characteristics_ch1", []))
gsms = g.gsms
rows, keep = [], []
for _, r in clin.iterrows():
    s = gsms.get(r["geo_accn"])
    if s is None: continue
    c = chars(s)
    try:
        rows.append({"t": float(r["time.rfs"]), "e": int(r["event.rfs"]),
                     "ggi": float(c["ggi"])})
        keep.append(s)
    except (KeyError, ValueError):
        pass
d = pd.DataFrame(rows)
print("usable:", len(d), "events:", int(d.e.sum()), flush=True)
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
out = {"cohort": "GSE2990", "n": int(len(d)), "events": int(d.e.sum()),
       "genes_mapped": mapped,
       "ours_expr_cindex": round(float(concordance_index(d.t.values, -eta, d.e.values)), 4),
       "ggi_continuous_cindex": round(float(concordance_index(d.t.values, -d.ggi.values, d.e.values)), 4),
       "note": "ours = expr-only Cox (70-gene panel) trained on METABRIC, transported raw; ggi = per-sample continuous Genomic Grade Index shipped with GSE2990"}
json.dump(out, open("results/transport2990_round9c.json", "w"), indent=2)
print(json.dumps(out, indent=1), flush=True)
