"""Discovery round 4: promoter methylation (RRBS) vs RFS.
(a) univariate Cox screen per gene (BH FDR);
(b) multivariate Cox: clin+expr vs clin+expr+meth, same 5 seeds/splits -> dC-index."""
import json, sys
import numpy as np
sys.path.insert(0, "src")
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from scipy import stats
from recurscan.data import cbioportal as cb
from recurscan.data.dataset import assemble

ds = assemble()
genes = cb.resolve_genes()
sym2ent = dict(genes)  # resolve_genes: symbol -> entrez
entrez_list = sorted(sym2ent.values())
idx_of = {e: i for i, e in enumerate(entrez_list)}
panel_entrez = [sym2ent[g] for g in ds.gene_names if g in sym2ent]

rows = cb.fetch_methylation(entrez_list, ds.sample_ids)
sid_idx = {s: i for i, s in enumerate(ds.sample_ids)}
M = np.full((len(entrez_list), len(ds.sample_ids)), np.nan)
for r in rows:
    e, s = r.get("entrezGeneId"), r.get("sampleId")
    if e in idx_of and s in sid_idx:
        M[idx_of[e], sid_idx[s]] = r["value"]
# keep panel genes with >=50% coverage
keep = [i for i, e in enumerate(entrez_list) if np.isfinite(M[i]).mean() >= 0.5 and e in set(panel_entrez)]
ent2sym = {v: k for k, v in sym2ent.items()}
names = [ent2sym[entrez_list[i]] for i in keep]
M = M[keep]
col_mean = np.nanmean(M, axis=1, keepdims=True)
M = np.where(np.isfinite(M), M, col_mean)  # mean-impute
Mz = (M - M.mean(1, keepdims=True)) / (M.std(1, keepdims=True) + 1e-9)
X_meth = Mz.T  # (n, k)
print("methylation matrix:", X_meth.shape, "genes kept:", len(names))

t, e = ds.time, ds.event
# (a) univariate screen
pvals = []
for j in range(X_meth.shape[1]):
    df = {"t": t, "e": e, "x": X_meth[:, j]}
    try:
        fit = CoxPHFitter(penalizer=1e-4).fit(df, "t", "e")
        pvals.append(fit.summary.loc["x", "p"])
    except Exception:
        pvals.append(1.0)
pvals = np.array(pvals)
order = np.argsort(pvals)
m = len(pvals)
bh = pvals[order] * m / (np.arange(m) + 1)
sig = [(names[i], float(pvals[i]), float(bh[r])) for r, i in enumerate(order) if bh[r] < 0.05]
print("FDR<0.05 hits:", sig[:10])

# (b) combined model, 5 seeds
rng_base = 20260924
res = {"base": [], "comb": []}
for seed in range(5):
    rng = np.random.default_rng(rng_base + seed)
    perm = rng.permutation(len(t))
    te = perm[: len(t) // 5]
    tr = perm[len(t) // 5:]
    def fit_predict(Xtr, Xte, tag):
        dtr = {"t": t[tr], "e": e[tr]}
        for j in range(Xtr.shape[1]): dtr[f"x{j}"] = Xtr[:, j]
        fit = CoxPHFitter(penalizer=0.05).fit(pd.DataFrame(dtr) if False else __import__("pandas").DataFrame(dtr), "t", "e")
        dte = {f"x{j}": Xte[:, j] for j in range(Xte.shape[1])}
        risk = fit.predict_partial_hazard(__import__("pandas").DataFrame(dte)).values.ravel()
        return concordance_index(t[te], -risk, e[te])
    Xb = np.hstack([ds.X_clin, ds.X_expr])
    Xc = np.hstack([ds.X_clin, ds.X_expr, X_meth])
    res["base"].append(fit_predict(Xb[tr], Xb[te], "b"))
    res["comb"].append(fit_predict(Xc[tr], Xc[te], "c"))
    print(f"seed {seed}: base {res['base'][-1]:.4f} comb {res['comb'][-1]:.4f}", flush=True)
out = {"n": int(len(t)), "events": int(e.sum()), "meth_genes": len(names),
       "fdr_hits": sig, "base_ci": res["base"], "comb_ci": res["comb"],
       "base_mean": float(np.mean(res["base"])), "comb_mean": float(np.mean(res["comb"])),
       "delta": float(np.mean(res["comb"]) - np.mean(res["base"])),
       "paired_t_p": float(stats.ttest_rel(res["comb"], res["base"]).pvalue)}
json.dump(out, open("results/methylation_round4.json", "w"), indent=2)
print(json.dumps({k: out[k] for k in ("base_mean", "comb_mean", "delta", "paired_t_p")}, indent=1))
