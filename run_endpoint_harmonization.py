"""Item 3 execution: patient-level endpoint harmonization
(docs/PREREG_ENDPOINT_HARMONIZATION_20260928.md). Frozen protocol:
- model: penalized Cox (0.05), METABRIC 70-gene expr panel, fit once here
  (identical fitting code to run_transport.py, before any cohort scoring);
- endpoints derived ONLY from the committed raw pulls in
  data_cache/external/geo/clinical_*.csv (+ METABRIC committed cache);
- strata HARM-RFS / HARM-DMFS per prereg; no cross-stratum pooling;
- statistic: per-cohort C-index (unit-free) + bootstrap CI (500, seed 0,
  recurscan.eval.bootstrap_cindex_ci, identical to committed transports);
- pooled: n-weighted mean of per-cohort C-indices within a stratum, CI by
  cohort-stratified patient bootstrap (seed 0, 500 reps).
"""
import json, sys
import numpy as np, pandas as pd
sys.path.insert(0, "src")
import GEOparse
from lifelines import CoxPHFitter
from recurscan.data.dataset import assemble
from recurscan.eval import bootstrap_cindex_ci

GEO = "data_cache/external/geo"
ds = assemble()
genes_in = list(ds.gene_names)
Xe = ds.X_expr.astype(float)
df = {f"x{j}": Xe[:, j] for j in range(Xe.shape[1])}
df["t"], df["e"] = ds.time.astype(float), ds.event.astype(int)
fit = CoxPHFitter(penalizer=0.05).fit(pd.DataFrame(df), "t", "e")
coef = fit.summary["coef"].values
print("METABRIC expr-only Cox fit: n", Xe.shape[0], "genes", len(genes_in), flush=True)

def cohort_expr(gse, gpl_name):
    g = GEOparse.get_GEO(geo=gse, destdir=GEO, annotate_gpl=False, silent=True)
    gpl = GEOparse.get_GEO(geo=gpl_name, destdir=GEO, silent=True)
    p2g = dict(zip(gpl.table["ID"], gpl.table["Gene Symbol"]))
    gsms = list(g.gsms.values())
    gids = np.array([p2g.get(i, "") for i in gsms[0].table["ID_REF"].values])
    X = np.column_stack([s.table["VALUE"].values for s in gsms])
    names = [s.name for s in gsms]
    Xg = np.zeros((len(gsms), len(genes_in)))
    mapped = 0
    for j, gene in enumerate(genes_in):
        m = gids == gene
        if m.sum():
            Xg[:, j] = X[m].mean(0)
            mapped += 1
    return names, Xg, mapped

# endpoint derivations from committed raw CSVs; each entry records derivation
def load_clin(gse):
    return pd.read_csv(f"{GEO}/clinical_{gse}.csv", dtype=str).set_index("geo_accn")

def num(s):
    return pd.to_numeric(s, errors="coerce")

DERIV = {}
def endpoints(gse, clin):
    """Return dict stratum -> DataFrame(t, e) indexed by geo_accn."""
    out = {}
    if gse == "GSE2990":
        out["HARM-RFS"] = pd.DataFrame({"t": num(clin["time.rfs"]), "e": num(clin["event.rfs"])}, index=clin.index)
        out["HARM-DMFS"] = pd.DataFrame({"t": num(clin["time.dmfs"]), "e": num(clin["event.dmfs"])}, index=clin.index)
        DERIV[gse] = {"HARM-RFS": "event.rfs/time.rfs verbatim (years; unit-free C-index)",
                      "HARM-DMFS": "event.dmfs/time.dmfs verbatim (years)"}
    elif gse == "GSE7390":
        out["HARM-RFS"] = pd.DataFrame({"t": num(clin["t.rfs"]), "e": num(clin["e.rfs"])}, index=clin.index)
        out["HARM-DMFS"] = pd.DataFrame({"t": num(clin["t.dmfs"]), "e": num(clin["e.dmfs"])}, index=clin.index)
        DERIV[gse] = {"HARM-RFS": "e.rfs/t.rfs verbatim", "HARM-DMFS": "e.dmfs/t.dmfs verbatim"}
    elif gse == "GSE11121":
        out["HARM-DMFS"] = pd.DataFrame({"t": num(clin["t.dmfs"]), "e": num(clin["e.dmfs"])}, index=clin.index)
        DERIV[gse] = {"HARM-DMFS": "t.dmfs/e.dmfs verbatim (months)",
                      "HARM-RFS": "EXCLUDED: no any-recurrence field in committed raw pull"}
    elif gse == "GSE20685":
        met = num(clin["event_metastasis"])
        reg = num(clin["regional_relapse"])
        t = num(clin["follow_up_duration (years)"])
        out["HARM-DMFS"] = pd.DataFrame({"t": t, "e": met}, index=clin.index)
        DERIV[gse] = {"HARM-DMFS": "event_metastasis + follow_up_duration (years) verbatim"}
        if set(reg.dropna().unique()) <= {0, 1}:
            rfs = ((reg == 1) | (met == 1)).astype(float).where(reg.notna() & met.notna())
            out["HARM-RFS"] = pd.DataFrame({"t": t, "e": rfs}, index=clin.index)
            DERIV[gse]["HARM-RFS"] = "any recurrence = regional_relapse OR event_metastasis (both 0/1 in committed pull); same follow-up field"
        else:
            DERIV[gse]["HARM-RFS"] = "EXCLUDED: regional_relapse not a clean 0/1 indicator in committed raw pull"
    elif gse == "GSE25066":
        out["HARM-DMFS"] = pd.DataFrame({"t": num(clin["drfs_even_time_years"]), "e": num(clin["drfs_1_event_0_censored"])}, index=clin.index)
        DERIV[gse] = {"HARM-DMFS": "drfs_1_event_0_censored + drfs_even_time_years verbatim (distant relapse)",
                      "HARM-RFS": "EXCLUDED: no any-recurrence field in committed raw pull (neoadjuvant pCR cohort)"}
    return out

GPL = {"GSE2990": "GPL96", "GSE7390": "GPL96", "GSE11121": "GPL96",
       "GSE20685": "GPL570", "GSE25066": "GPL96"}

per = []
for gse in ["GSE2990", "GSE7390", "GSE11121", "GSE20685", "GSE25066"]:
    names, Xg, mapped = cohort_expr(gse, GPL[gse])
    clin = load_clin(gse)
    for stratum, ep in endpoints(gse, clin).items():
        d = ep.dropna()
        common = [s for s in names if s in d.index]
        idx = [names.index(s) for s in common]
        if len(common) < 30:
            DERIV[gse][stratum] += f" | DROPPED: only {len(common)} matched patients"
            continue
        Xs = Xg[idx]
        Xs = (Xs - Xs.mean(0)) / (Xs.std(0) + 1e-9)  # per-gene z within analysis set, as committed
        eta = Xs @ coef
        tv = d.loc[common, "t"].values.astype(float)
        ev = d.loc[common, "e"].values.astype(int)
        lo, hi = bootstrap_cindex_ci(tv, ev, eta, n_boot=500, seed=0)
        from lifelines.utils import concordance_index as cidx
        c = cidx(tv, -eta, ev)
        per.append({"cohort": gse, "stratum": stratum, "n": int(len(common)),
                    "events": int(ev.sum()), "c_index": round(float(c), 4),
                    "ci95": [round(float(lo), 4), round(float(hi), 4)],
                    "genes_mapped": int(mapped)})
        print(gse, stratum, "n", len(common), "events", int(ev.sum()),
              "C", round(float(c), 4), flush=True)

# pooled within stratum: n-weighted mean of per-cohort C; cohort-stratified bootstrap
rng = np.random.default_rng(0)
pooled = {}
for stratum in ["HARM-RFS", "HARM-DMFS"]:
    rows = [r for r in per if r["stratum"] == stratum]
    if not rows:
        continue
    ws = np.array([r["n"] for r in rows], float)
    cs = np.array([r["c_index"] for r in rows])
    point = float((ws * cs).sum() / ws.sum())
    # bootstrap: resample patients within each cohort, recompute per-cohort C
    boots = []
    caches = {}
    for r in rows:
        gse = r["cohort"]
        if gse not in caches:
            names, Xg, _ = cohort_expr(gse, GPL[gse])
            clin = load_clin(gse)
            ep = endpoints(gse, clin)[stratum].dropna()
            common = [s for s in names if s in ep.index]
            idx = [names.index(s) for s in common]
            Xs = Xg[idx]; Xs = (Xs - Xs.mean(0)) / (Xs.std(0) + 1e-9)
            caches[gse] = (ep.loc[common, "t"].values.astype(float),
                           ep.loc[common, "e"].values.astype(int), Xs @ coef)
    for _ in range(500):
        bcs, bws = [], []
        for r in rows:
            tv, ev, eta = caches[r["cohort"]]
            bi = rng.integers(0, len(tv), len(tv))
            from lifelines.utils import concordance_index as cidx
            bc = cidx(tv[bi], -eta[bi], ev[bi])
            if not np.isnan(bc):
                bcs.append(bc); bws.append(r["n"])
        if bcs:
            boots.append(float((np.array(bws) * np.array(bcs)).sum() / np.array(bws).sum()))
    lo, hi = np.percentile(boots, [2.5, 97.5])
    pooled[stratum] = {"cohorts": [r["cohort"] for r in rows],
                       "pooled_c_index": round(point, 4),
                       "ci95": [round(float(lo), 4), round(float(hi), 4)],
                       "definition": "n-weighted mean of per-cohort C-index; cohort-stratified patient bootstrap (500, seed 0)"}
    print("pooled", stratum, round(point, 4), [round(float(lo),4), round(float(hi),4)], flush=True)

out = {"prereg": "docs/PREREG_ENDPOINT_HARMONIZATION_20260928.md",
       "model": "penalized Cox 0.05, METABRIC 70-gene expr-only, fit once (run_transport.py protocol)",
       "per_cohort": per, "pooled": pooled, "derivation_audit": DERIV,
       "metabric_discovery": {"stratum": "HARM-RFS", "n": int(Xe.shape[0]),
                              "events": int(ds.event.sum()),
                              "note": "training cohort; in-sample fit, not a transport result"}}
json.dump(out, open("results/endpoint_harmonization.json", "w"), indent=1)
print("wrote results/endpoint_harmonization.json")
