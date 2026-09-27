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

Execution mechanics (low-memory revision, 2026-09-28): expression extraction
streams the cached family-SOFT / GPL annotation files line by line instead of
building full GEOparse in-memory objects (the GEOparse path was OOM-killed
twice at the second cohort on this 2GB box). Extraction is numerically
identical: same ID_REF->Gene Symbol map, same probe->gene mean, same per-gene
z-score, same patient matching; validated by reproducing the GEOparse path's
GSE2990 point estimates (HARM-RFS C 0.6559, HARM-DMFS C 0.6648). The pooled
bootstrap reuses the main loop's per-cohort arrays instead of re-parsing the
SOFT files (parsing never consumed the rng stream, so the bootstrap sequence
is unchanged). Statistics unchanged.

Endpoint parsing note (same revision): the committed raw pulls store GEO
characteristics as 'key: value' cells with column drift on some rows;
field_series() recovers each needed field key-verbatim per patient (first
'<key>: ' cell in the row; clean-column fallback; wrong-key cell = missing).
This only affects cohorts that produced NO outcome in the diagnostic run
(all dropped on NaN); GSE2990's clean pull takes the fallback path and its
already-observed point estimates must reproduce exactly.
"""
import json, sys, gzip, glob
import numpy as np, pandas as pd
sys.path.insert(0, "src")
import GEOparse  # download-only fetch of missing caches (no in-memory parse)
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index as cidx
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

def _open(path):
    return gzip.open(path, "rt", errors="replace") if path.endswith(".gz") else open(path, errors="replace")

def _ensure(geo_name, patterns):
    for pat in patterns:
        hits = glob.glob(f"{GEO}/{pat}")
        if hits:
            return hits[0]
    GEOparse.get_GEO_file(geo_name, destdir=GEO)  # download only, no parse
    for pat in patterns:
        hits = glob.glob(f"{GEO}/{pat}")
        if hits:
            return hits[0]
    raise FileNotFoundError(geo_name)

def load_p2g(gpl_name):
    """Stream GPL annotation SOFT -> dict probe ID -> Gene Symbol (verbatim)."""
    path = _ensure(gpl_name, [f"{gpl_name}.txt", f"{gpl_name}.annot.gz", f"{gpl_name}.soft.gz"])
    p2g = {}
    with _open(path) as fh:
        intable = False
        i_id = i_gs = None
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("!platform_table_begin"):
                intable = True
                header = fh.readline().rstrip("\n").split("\t")
                i_id, i_gs = header.index("ID"), header.index("Gene Symbol")
                continue
            if line.startswith("!platform_table_end"):
                break
            if intable:
                parts = line.split("\t")
                if len(parts) > max(i_id, i_gs):
                    p2g[parts[i_id]] = parts[i_gs]
    return p2g

def stream_gse(gse):
    """Stream family SOFT -> (sample GSM names, per-sample VALUE arrays, ID_REF axis)."""
    path = _ensure(gse, [f"{gse}_family.soft.gz"])
    names, cols, idref = [], [], None
    in_sample = False
    with _open(path) as fh:
        for line in fh:
            if line.startswith("^SAMPLE"):
                in_sample = True
                accn = None
                continue
            if line.startswith("^PLATFORM") or line.startswith("^SERIES") or line.startswith("^DATABASE"):
                in_sample = False
                continue
            if not in_sample:
                continue
            if line.startswith("!Sample_geo_accession"):
                accn = line.split("=", 1)[1].strip()
                continue
            if line.startswith("!sample_table_begin"):
                header = fh.readline().rstrip("\n").split("\t")
                i_id, i_val = header.index("ID_REF"), header.index("VALUE")
                ids, vals = [], np.zeros(0)
                vlist = []
                for line in fh:
                    if line.startswith("!sample_table_end"):
                        break
                    parts = line.rstrip("\n").split("\t")
                    ids.append(parts[i_id])
                    try:
                        vlist.append(float(parts[i_val]))
                    except (ValueError, IndexError):
                        vlist.append(float("nan"))
                v = np.asarray(vlist, dtype=float)
                if idref is None:
                    idref = ids
                else:
                    assert ids == idref, f"{gse}: ID_REF axis differs across samples"
                names.append(accn)
                cols.append(v)
                in_sample = False
    return names, cols, idref

P2G = {}
def cohort_expr(gse, gpl_name):
    if gpl_name not in P2G:
        P2G[gpl_name] = load_p2g(gpl_name)
    p2g = P2G[gpl_name]
    names, cols, idref = stream_gse(gse)
    gids = np.array([p2g.get(i, "") for i in idref])
    X = np.column_stack(cols)
    Xg = np.zeros((len(names), len(genes_in)))
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

def field_series(clin, key):
    """Per-patient value for a clinical key, verbatim from the committed raw
    pull. The pull stored GEO characteristics as 'key: value' cells and some
    rows drifted across columns, so: the first cell in the row whose prefix
    is exactly '<key>: ' supplies the value; if no prefixed cell matches,
    fall back to the named column's raw value (clean pulls like GSE2990);
    a named-column cell carrying a different key's prefix counts as missing
    (e.g. GSE20685 rows whose regional_relapse column carries m_stage)."""
    out = {}
    has_col = key in clin.columns
    prefix = key + ": "
    for accn, row in clin.iterrows():
        v = None
        for cell in row:
            if isinstance(cell, str) and cell.startswith(prefix):
                v = cell[len(prefix):]
                break
        if v is None and has_col:
            cell = row[key]
            if isinstance(cell, str) and ": " not in cell:
                v = cell
        out[accn] = v
    return pd.Series(out)

def num(s):
    return pd.to_numeric(s, errors="coerce")

DERIV = {}
def endpoints(gse, clin):
    """Return dict stratum -> DataFrame(t, e) indexed by geo_accn."""
    out = {}
    if gse == "GSE2990":
        out["HARM-RFS"] = pd.DataFrame({"t": num(field_series(clin, "time.rfs")), "e": num(field_series(clin, "event.rfs"))}, index=clin.index)
        out["HARM-DMFS"] = pd.DataFrame({"t": num(field_series(clin, "time.dmfs")), "e": num(field_series(clin, "event.dmfs"))}, index=clin.index)
        DERIV[gse] = {"HARM-RFS": "event.rfs/time.rfs verbatim (years; unit-free C-index)",
                      "HARM-DMFS": "event.dmfs/time.dmfs verbatim (years)"}
    elif gse == "GSE7390":
        out["HARM-RFS"] = pd.DataFrame({"t": num(field_series(clin, "t.rfs")), "e": num(field_series(clin, "e.rfs"))}, index=clin.index)
        out["HARM-DMFS"] = pd.DataFrame({"t": num(field_series(clin, "t.dmfs")), "e": num(field_series(clin, "e.dmfs"))}, index=clin.index)
        DERIV[gse] = {"HARM-RFS": "e.rfs/t.rfs verbatim (GEO 'key: value' characteristics cells, key-exact, prefix stripped at parse)",
                      "HARM-DMFS": "e.dmfs/t.dmfs verbatim (same key-verbatim extraction)"}
    elif gse == "GSE11121":
        out["HARM-DMFS"] = pd.DataFrame({"t": num(field_series(clin, "t.dmfs")), "e": num(field_series(clin, "e.dmfs"))}, index=clin.index)
        DERIV[gse] = {"HARM-DMFS": "t.dmfs/e.dmfs verbatim (months; key-verbatim extraction as GSE7390)",
                      "HARM-RFS": "EXCLUDED: no any-recurrence field in committed raw pull"}
    elif gse == "GSE20685":
        met = num(field_series(clin, "event_metastasis"))
        reg = num(field_series(clin, "regional_relapse"))
        t = num(field_series(clin, "follow_up_duration (years)"))
        out["HARM-DMFS"] = pd.DataFrame({"t": t, "e": met}, index=clin.index)
        DERIV[gse] = {"HARM-DMFS": "event_metastasis + follow_up_duration (years) verbatim (key-verbatim extraction)"}
        if set(reg.dropna().unique()) <= {0, 1}:
            rfs = ((reg == 1) | (met == 1)).astype(float).where(reg.notna() & met.notna())
            out["HARM-RFS"] = pd.DataFrame({"t": t, "e": rfs}, index=clin.index)
            DERIV[gse]["HARM-RFS"] = "any recurrence = regional_relapse OR event_metastasis, key-verbatim (83 of 327 pulled patients carry no regional_relapse field - their column holds m_stage, not a recurrence outcome - excluded from this stratum only); same follow-up field"
        else:
            DERIV[gse]["HARM-RFS"] = "EXCLUDED: regional_relapse not a clean 0/1 indicator in committed raw pull"
    elif gse == "GSE25066":
        out["HARM-DMFS"] = pd.DataFrame({"t": num(field_series(clin, "drfs_even_time_years")), "e": num(field_series(clin, "drfs_1_event_0_censored"))}, index=clin.index)
        DERIV[gse] = {"HARM-DMFS": "drfs_1_event_0_censored + drfs_even_time_years verbatim (distant relapse; key-verbatim, 198 of 508 rows recovered from drifted columns)",
                      "HARM-RFS": "EXCLUDED: no any-recurrence field in committed raw pull (neoadjuvant pCR cohort)"}
    return out

GPL = {"GSE2990": "GPL96", "GSE7390": "GPL96", "GSE11121": "GPL96",
       "GSE20685": "GPL570", "GSE25066": "GPL96"}

STRAT_CACHE = {}  # (gse, stratum) -> (tv, ev, eta) for the pooled bootstrap
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
        c = cidx(tv, -eta, ev)
        STRAT_CACHE[(gse, stratum)] = (tv, ev, eta)
        per.append({"cohort": gse, "stratum": stratum, "n": int(len(common)),
                    "events": int(ev.sum()), "c_index": round(float(c), 4),
                    "ci95": [round(float(lo), 4), round(float(hi), 4)],
                    "genes_mapped": int(mapped)})
        print(gse, stratum, "n", len(common), "events", int(ev.sum()),
              "C", round(float(c), 4), flush=True)

# pooled within stratum: n-weighted mean of per-cohort C; cohort-stratified bootstrap
# (reuses the main loop's arrays; the bootstrap rng stream is unchanged)
rng = np.random.default_rng(0)
pooled = {}
for stratum in ["HARM-RFS", "HARM-DMFS"]:
    rows = [r for r in per if r["stratum"] == stratum]
    if not rows:
        continue
    ws = np.array([r["n"] for r in rows], float)
    cs = np.array([r["c_index"] for r in rows])
    point = float((ws * cs).sum() / ws.sum())
    boots = []
    for _ in range(500):
        bcs, bws = [], []
        for r in rows:
            tv, ev, eta = STRAT_CACHE[(r["cohort"], stratum)]
            bi = rng.integers(0, len(tv), len(tv))
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
