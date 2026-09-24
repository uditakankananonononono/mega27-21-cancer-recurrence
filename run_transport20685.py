"""Round 9e: transport METABRIC expr-only Cox to GSE20685 (Kao 2011, GPL570).
Streaming soft-file parse (GEOparse OOMs on 327x54k). Endpoint: metastasis."""
import json, sys, gzip
import numpy as np, pandas as pd
sys.path.insert(0, "src")
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from recurscan.data.dataset import assemble
from recurscan.eval import bootstrap_cindex_ci

ds = assemble()
genes_in = list(ds.gene_names)
want = set(genes_in)
Xe = ds.X_expr.astype(float)
t, e = ds.time.astype(float), ds.event.astype(int)
df = {f"x{j}": Xe[:, j] for j in range(Xe.shape[1])}
df["t"], df["e"] = t, e
fit = CoxPHFitter(penalizer=0.05).fit(pd.DataFrame(df), "t", "e")
coef = fit.summary["coef"].values

path = "data_cache/external/geo/GSE20685_family.soft.gz"
probe2gene = {}
samples = []  # (gsm, chars-dict, {gene: [sum, count]})
cur = None; mode = None; tbl_header = None
with gzip.open(path, "rt", errors="replace") as fh:
    for line in fh:
        line = line.rstrip("\n")
        if line.startswith("^SAMPLE"):
            cur = {"gsm": line.split("=")[1].strip(), "chars": {}, "genes": {}}
            samples.append(cur); mode = None
        elif line.startswith("!platform_table_begin"):
            mode = "platform"
        elif line.startswith("!platform_table_end"):
            mode = None
        elif mode == "platform":
            if tbl_header is None:
                tbl_header = line.lstrip("#").split("\t")
                continue
            parts = line.split("\t")
            try:
                pid = parts[0]
                sym = parts[tbl_header.index("Gene Symbol")].split(" ///")[0].strip()
            except (ValueError, IndexError):
                continue
            if sym in want:
                probe2gene[pid] = sym
        elif line.startswith("!Sample_characteristics_ch1"):
            v = line.split("=", 1)[1].strip()
            k, _, val = v.partition(": ")
            cur["chars"][k.strip()] = val.strip()
        elif line.startswith("!sample_table_begin"):
            mode = "sample"
        elif line.startswith("!sample_table_end"):
            mode = None
        elif mode == "sample" and not line.startswith("ID_REF"):
            parts = line.split("\t")
            if len(parts) >= 2:
                g = probe2gene.get(parts[0])
                if g:
                    try:
                        val = float(parts[1])
                    except ValueError:
                        continue
                    s = cur["genes"].setdefault(g, [0.0, 0])
                    s[0] += val; s[1] += 1
rows, Xrows = [], []
for s in samples:
    c = s["chars"]
    try:
        rows.append({"t": float(c["follow_up_duration (years)"]),
                     "e": int(c["event_metastasis"]),
                     "n_stage": float(c["n_stage"])})
        Xrows.append([s["genes"].get(g, [np.nan, 0])[0] / s["genes"][g][1] if g in s["genes"] and s["genes"][g][1] else np.nan for g in genes_in])
    except (KeyError, ValueError):
        pass
d = pd.DataFrame(rows)
Xg = np.array(Xrows)
print("usable:", len(d), "events:", int(d.e.sum()), "genes mapped:", int((~np.isnan(Xg).any(0)).sum()), flush=True)
Xg = np.nan_to_num(Xg, nan=0.0)
for j in range(Xg.shape[1]):
    col = Xg[:, j]
    sd = col.std()
    Xg[:, j] = (col - col.mean()) / (sd + 1e-9)
eta = Xg @ coef
tv, ev = d.t.values, d.e.values
def boot(sc):
    lo, hi = bootstrap_cindex_ci(tv, ev, sc, n_boot=500, seed=0)
    return [round(float(lo), 4), round(float(hi), 4)]
out = {"cohort": "GSE20685", "endpoint": "metastasis (years)", "n": int(len(d)), "events": int(ev.sum()),
       "ours_expr": {"ci": round(float(concordance_index(tv, -eta, ev)), 4), "boot95": boot(eta)},
       "n_stage": {"ci": round(float(concordance_index(tv, -d.n_stage.values, ev)), 4), "boot95": boot(d.n_stage.values)},
       "note": "ours = expr-only Cox (70-gene panel) trained on METABRIC, transported raw to GPL570 cohort; comparator = nodal stage"}
json.dump(out, open("results/transport20685_round9e.json", "w"), indent=2)
print(json.dumps(out, indent=1), flush=True)
