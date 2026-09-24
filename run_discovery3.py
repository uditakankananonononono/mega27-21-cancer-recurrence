"""NPI-discordant early recurrers: does a molecular signature separate
low-NPI patients who recur within 5y from those who stay disease-free?

Locked gate BEFORE looking: candidate wins only if nested-CV logistic AUC
>= 0.65 with bootstrap CI excluding 0.5, AND >= 1 univariate marker passes
BH-FDR 0.05. Otherwise honest negative.
"""
import json, sys
sys.path.insert(0, "src")
import numpy as np
from scipy import stats
from sklearn.linear_model import LogisticRegressionCV
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

from recurscan.data.dataset import assemble

ds = assemble()
clin = json.load(open("data_cache/clinical.json"))
npi = np.array([float(clin[s].get("NPI", "nan")) for s in ds.sample_ids])
low = npi <= 3.4
early_recur = low & (ds.event == 1) & (ds.time <= 60)
clean_free = low & (ds.event == 0) & (ds.time >= 60)
sel = early_recur | clean_free
y = early_recur[sel].astype(int)
X = ds.X_expr[sel]
print(f"low-NPI cohort: {sel.sum()} | early recurrers {y.sum()} | disease-free {(y==0).sum()}", flush=True)
if y.sum() < 20 or (y == 0).sum() < 20:
    print("COHORT TOO SMALL", flush=True); sys.exit(0)

# univariate screen with BH-FDR
pvals = []
for j, g in enumerate(ds.gene_names):
    a, b = X[y == 1, j], X[y == 0, j]
    pvals.append((g, stats.mannwhitneyu(a, b, alternative="two-sided").pvalue))
pvals.sort(key=lambda t: t[1])
m = len(pvals)
bh = [(g, p, p * m / (i + 1)) for i, (g, p) in enumerate(pvals)]
hits = [(g, p, q) for (g, p, q) in bh if q < 0.05]
print(f"univariate BH-FDR<0.05 hits: {len(hits)}", flush=True)
for g, p, q in hits[:15]:
    a = X[y == 1, ds.gene_names.index(g)]
    b = X[y == 0, ds.gene_names.index(g)]
    print(f"  {g}: p={p:.2e} q={q:.2e} mean recur {a.mean():+.2f} vs free {b.mean():+.2f}", flush=True)

# honest nested estimate: outer 5-fold AUC of L1-logistic
outer = StratifiedKFold(5, shuffle=True, random_state=0)
aucs = []
for tr, te in outer.split(X, y):
    sc = StandardScaler().fit(X[tr])
    clf = LogisticRegressionCV(Cs=10, penalty="l1", solver="liblinear",
                               max_iter=3000, cv=3)
    clf.fit(sc.transform(X[tr]), y[tr])
    aucs.append(roc_auc_score(y[te], clf.predict_proba(sc.transform(X[te]))[:, 1]))
aucs = np.array(aucs)
print(f"nested-CV AUC: {aucs.mean():.3f} +- {aucs.std():.3f} (folds {np.round(aucs,3)})", flush=True)
# bootstrap CI on out-of-fold predictions
sc = StandardScaler().fit(X)
oof = np.zeros(len(y))
for tr, te in outer.split(X, y):
    s2 = StandardScaler().fit(X[tr])
    clf = LogisticRegressionCV(Cs=10, penalty="l1", solver="liblinear", max_iter=3000, cv=3)
    clf.fit(s2.transform(X[tr]), y[tr])
    oof[te] = clf.predict_proba(s2.transform(X[te]))[:, 1]
rng = np.random.default_rng(0)
boots = [roc_auc_score(y[i], oof[i]) for i in (rng.integers(0, len(y), len(y)) for _ in range(500))]
print(f"OOF AUC {roc_auc_score(y, oof):.3f} CI {np.percentile(boots,[2.5,97.5]).round(3)}", flush=True)
json.dump({"cohort_n": int(sel.sum()), "recurrers": int(y.sum()),
           "fdr_hits": [(g, float(p), float(q)) for g, p, q in hits],
           "fold_aucs": aucs.tolist(), "oof_auc": float(roc_auc_score(y, oof)),
           "oof_auc_ci": np.percentile(boots, [2.5, 97.5]).tolist()},
          open("results/discovery_round3_npi_discordant.json", "w"), indent=2)
print("SAVED", flush=True)
