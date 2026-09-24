"""Round 9: cross-cohort transportability benchmark vs PUBLISHED signatures.
Train Cox (NPI + 70-gene expr) on METABRIC -> score GSE7390 (Affymetrix,
GPL96-mapped) -> compare held C-index against GSE7390's own published risk
scores: veridex_risk (76-gene MammaPrint-like) and risksg (GGI)."""
import json, sys
import numpy as np, pandas as pd
sys.path.insert(0, "src")
import GEOparse
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from recurscan.data.dataset import assemble

ds = assemble()
genes_in = ds.gene_names
Xb = np.hstack([ds.X_clin[:, [ds.clin_names.index("NPI")]], ds.X_expr]).astype(float)
t, e = ds.time.astype(float), ds.event.astype(int)
df = {f"x{j}": Xb[:, j] for j in range(Xb.shape[1])}
df["t"], df["e"] = t, e
fit = CoxPHFitter(penalizer=0.05).fit(pd.DataFrame(df), "t", "e")
coef = fit.summary["coef"].values  # [NPI, 70 genes]

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
                     "npi": float(c["NPI"]),
                     "veridex": c["veridex_risk"], "ggi": c["risksg"]})
        keep.append(s)
    except (KeyError, ValueError):
        pass
d = pd.DataFrame(rows)
print("usable:", len(d), "events:", int(d.e.sum()))
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
npi_z = (d.npi.values - d.npi.values.mean()) / d.npi.values.std()
eta = coef[0] * npi_z + Xg @ coef[1:]
# comparators: veridex_risk / risksg are categorical (Good/Poor or High/Low)
def cat_risk(s):
    u = {v: i for i, v in enumerate(sorted(set(s)))}
    return np.array([u[v] for v in s], dtype=float)
out = {"cohort": "GSE7390", "n": int(len(d)), "events": int(d.e.sum()),
       "genes_mapped": mapped,
       "ours_cindex": round(float(concordance_index(d.t.values, eta, d.e.values)), 4),
       "veridex76_cindex": round(float(concordance_index(d.t.values, cat_risk(d.veridex), d.e.values)), 4),
       "ggi_cindex": round(float(concordance_index(d.t.values, cat_risk(d.ggi), d.e.values)), 4),
       "note": "ours = Cox(NPI + 70-gene expr) trained on METABRIC, transported raw to Affymetrix cohort"}
json.dump(out, open("results/transport_round9.json", "w"), indent=2)
print(json.dumps(out, indent=1))
