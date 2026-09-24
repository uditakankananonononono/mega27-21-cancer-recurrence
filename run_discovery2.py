import json, sys
sys.path.insert(0, "src")
import numpy as np
from recurscan.data import cbioportal as cb
from recurscan.data.dataset import assemble
from recurscan.benchmark import standardize, stratified_event_split
from recurscan.discovery import bootstrap_delta_ci, cox_risk, fit_cox
from recurscan.survival import concordance_index_fast

ds = assemble()
samples = json.load(open("data_cache/samples.json"))
clin = json.load(open("data_cache/clinical.json"))
genes = json.load(open("data_cache/genes.json"))
entrez = sorted({v: k for k, v in genes.items()})
gcol = {g: j for j, g in enumerate(entrez)}
sid_pos = {s: i for i, s in enumerate(ds.sample_ids)}

# CNA matrix aligned to ds samples
cna_rows = json.load(open("data_cache/cna.json"))
X_cna = np.zeros((len(ds.sample_ids), len(entrez)), dtype=np.float32)
for r in cna_rows:
    i = sid_pos.get(r["sampleId"])
    if i is not None and r["entrezGeneId"] in gcol:
        X_cna[i, gcol[r["entrezGeneId"]]] = r["value"]

# mutation indicators (nonsynonymous) aligned
mut_rows = json.load(open("data_cache/mutations.json"))
X_mut = np.zeros((len(ds.sample_ids), len(entrez)), dtype=np.float32)
for r in mut_rows:
    i = sid_pos.get(r["sampleId"])
    if i is not None and r["entrezGeneId"] in gcol:
        X_mut[i, gcol[r["entrezGeneId"]]] = 1.0
keep_mut = X_mut.sum(0) >= 20  # at least 20 mutated samples
X_mut = X_mut[:, keep_mut]
print("CNA features:", X_cna.shape[1], "| mutation features kept:", X_mut.shape[1], flush=True)

X_base = np.concatenate([ds.X_clin, ds.X_expr], axis=1)
results = {}
for seed in (0, 1, 2):
    tr, te = stratified_event_split(ds.event, 0.25, seed)
    Xtr, Xte = standardize(X_base[tr], X_base[te])
    base = fit_cox(Xtr, ds.time[tr], ds.event[tr])
    rb = cox_risk(base, Xte)
    out = {"base": concordance_index_fast(ds.time[te], ds.event[te], rb)}
    for name, Xa in (("cna", X_cna), ("mut", X_mut),
                     ("cna_mut", np.hstack([X_cna, X_mut]))):
        A_tr, A_te = standardize(Xa[tr], Xa[te])
        cand = fit_cox(np.hstack([Xtr, A_tr]), ds.time[tr], ds.event[tr])
        ra = cox_risk(cand, np.hstack([Xte, A_te]))
        d = bootstrap_delta_ci(ds.time[te], ds.event[te], ra, rb, seed=seed)
        out[name] = {"c": concordance_index_fast(ds.time[te], ds.event[te], ra),
                     "delta": d}
        print(f"seed {seed} +{name}: {out[name]['c']:.4f} (d {d[0]:+.4f} CI({d[1][0]:+.4f},{d[1][1]:+.4f}))", flush=True)
    results[seed] = out
json.dump(results, open("results/discovery_round2.json", "w"), indent=2, default=float)
print("SAVED", flush=True)
