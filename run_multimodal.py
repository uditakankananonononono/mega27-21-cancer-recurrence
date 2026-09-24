"""Benchmark-break attempt: multimodal Cox (clin+expr+CNA+mutation) vs the
clin+expr leader (0.6801). Leakage-safe: CNA/mut feature screening inside
each training fold only. 5 seeds, held-out C-index."""
import json, sys
import numpy as np, pandas as pd
sys.path.insert(0, "src")
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from recurscan.data import cbioportal as cb
from recurscan.data.dataset import assemble

ds = assemble()
genes = cb.resolve_genes()
entrez_list = sorted(genes.values())
sym2ent = dict(genes)
ent2sym = {v: k for k, v in genes.items()}
idx_of = {e: i for i, e in enumerate(entrez_list)}
sid_idx = {s: i for i, s in enumerate(ds.sample_ids)}
t, e = ds.time.astype(float), ds.event.astype(int)
Xb = np.hstack([ds.X_clin, ds.X_expr]).astype(float)

def matrix(rows, discrete=True):
    M = np.zeros((len(entrez_list), len(ds.sample_ids)), dtype=np.float32)
    for r in rows:
        g, s = r.get("entrezGeneId"), r.get("sampleId")
        if g in idx_of and s in sid_idx:
            v = r.get("value", 1.0)
            if discrete:  # mutation: any non-synonymous -> 1
                v = 1.0 if r.get("mutationType", "") not in ("Silent", "") else 0.0
            M[idx_of[g], sid_idx[s]] = max(M[idx_of[g], sid_idx[s]], v) if discrete else v
    return M

Mcna = matrix(cb.fetch_cna(entrez_list, ds.sample_ids), discrete=False)
Mmut = matrix(cb.fetch_mutations(entrez_list, ds.sample_ids), discrete=True)
panel_idx = [idx_of[sym2ent[g]] for g in ds.gene_names if g in sym2ent]
# CNA: use |value|>=2 (amp/deep-del) as binary events; mutations: binary
Acna = (np.abs(Mcna[panel_idx]) >= 2).astype(float).T   # (n, 70)
Amut = Mmut[panel_idx].T                                # (n, 70)
print("CNA event rate per gene (mean):", round(float(Acna.mean()), 4),
      "| mut event rate:", round(float(Amut.mean()), 4))

def screen(Xs, tr, topk=20):
    """Univariate Cox p-screen on TRAIN only; return kept column indices."""
    p = []
    for j in range(Xs.shape[1]):
        if Xs[tr, j].std() < 1e-8 or Xs[tr, j].mean() < 0.01:
            p.append(1.0); continue
        try:
            f = CoxPHFitter(penalizer=1e-4).fit(
                pd.DataFrame({"t": t[tr], "e": e[tr], "x": Xs[tr, j]}), "t", "e")
            p.append(float(f.summary.loc["x", "p"]))
        except Exception:
            p.append(1.0)
    p = np.asarray(p)
    valid = np.where(p < 1.0)[0]
    return valid[np.argsort(p[valid])][:topk]

res = {"base": [], "multi": []}
for seed in range(5):
    rng = np.random.default_rng(20260924 + seed)
    perm = rng.permutation(len(t))
    te = perm[: len(t) // 5]; tr = perm[len(t) // 5:]
    k_cna = screen(Acna, tr); k_mut = screen(Amut, tr)
    def run(Xtr, Xte):
        ok = Xtr.std(0) > 1e-8
        Xtr, Xte = Xtr[:, ok], Xte[:, ok]
        d = {f"x{j}": Xtr[:, j] for j in range(Xtr.shape[1])}
        d["t"], d["e"] = t[tr], e[tr]
        dte = pd.DataFrame({f"x{j}": Xte[:, j] for j in range(Xte.shape[1])})
        fit = None
        for pen in (0.05, 0.2, 0.5, 1.0):
            try:
                fit = CoxPHFitter(penalizer=pen).fit(pd.DataFrame(d), "t", "e")
                break
            except Exception:
                continue
        if fit is None:
            return float("nan")
        risk = fit.predict_partial_hazard(dte).values.ravel()
        return concordance_index(t[te], -risk, e[te])
    res["base"].append(run(Xb[tr], Xb[te]))
    Xm_tr = np.hstack([Xb[tr], Acna[tr][:, k_cna], Amut[tr][:, k_mut]])
    Xm_te = np.hstack([Xb[te], Acna[te][:, k_cna], Amut[te][:, k_mut]])
    res["multi"].append(run(Xm_tr, Xm_te))
    print(f"seed {seed}: base {res['base'][-1]:.4f} multi {res['multi'][-1]:.4f}", flush=True)
from scipy import stats
out = {"base_runs": [round(x, 4) for x in res["base"]],
       "multi_runs": [round(x, 4) for x in res["multi"]],
       "base_mean": round(float(np.mean(res["base"])), 4),
       "multi_mean": round(float(np.mean(res["multi"])), 4),
       "delta": round(float(np.mean(res["multi"]) - np.mean(res["base"])), 4),
       "paired_t_p": float(stats.ttest_rel(res["multi"], res["base"]).pvalue)}
json.dump(out, open("results/multimodal_round6.json", "w"), indent=2)
print(json.dumps(out, indent=1))
