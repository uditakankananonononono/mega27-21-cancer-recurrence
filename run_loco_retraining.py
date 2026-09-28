"""Item 11 execution: leave-cohort-out (LOCO) retraining across the harmonized
external cohorts (docs/PREREG_LOCO_RETRAINING_20260928.md). Frozen protocol:
- frozen comparator: METABRIC expr-only penalized Cox (0.05), identical fitting
  code to run_endpoint_harmonization.py; integrity gate aborts the run before
  any LOCO fit if recomputed frozen per-cohort C-indices, n, events or
  genes_mapped differ from the committed results/endpoint_harmonization.json;
- extraction helpers below are COPIED VERBATIM from run_endpoint_harmonization.py
  (same committed raw pulls, same key-verbatim endpoints, same per-gene z);
- LOCO arm per cohort: CoxPHFitter(penalizer=0.05) on the stratum's other
  cohorts pooled (each z-scored within itself), zero-variance genes dropped;
- statistic: concordance_index(t, -eta, e), the committed lane convention;
- delta_C = LOCO - committed frozen C-index; paired patient bootstrap B=1000,
  one rng stream seeded 0, cohorts in the fixed prereg order.
Numbers only; the verdict mapping lives in the queue and the paper."""
import json, sys
import numpy as np, pandas as pd
sys.path.insert(0, "src")
import GEOparse  # download-only fetch of missing caches (no in-memory parse)
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index as cidx
from recurscan.data.dataset import assemble

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
            DERIV[gse]["HARM-RFS"] = "any recurrence = regional_relapse OR event_metastasis, key-verbatim (20 of 327 pulled patients carry regional_relapse: NA, no usable value - excluded from this stratum only, leaving n=307); same follow-up field"
        else:
            DERIV[gse]["HARM-RFS"] = "EXCLUDED: regional_relapse not a clean 0/1 indicator in committed raw pull"
    elif gse == "GSE25066":
        out["HARM-DMFS"] = pd.DataFrame({"t": num(field_series(clin, "drfs_even_time_years")), "e": num(field_series(clin, "drfs_1_event_0_censored"))}, index=clin.index)
        DERIV[gse] = {"HARM-DMFS": "drfs_1_event_0_censored + drfs_even_time_years verbatim (distant relapse; key-verbatim, 198 of 508 rows recovered from drifted columns)",
                      "HARM-RFS": "EXCLUDED: no any-recurrence field in committed raw pull (neoadjuvant pCR cohort)"}
    return out

GPL = {"GSE2990": "GPL96", "GSE7390": "GPL96", "GSE11121": "GPL96",
       "GSE20685": "GPL570", "GSE25066": "GPL96"}
STRATA = {"HARM-RFS": ["GSE2990", "GSE7390", "GSE20685"],
          "HARM-DMFS": ["GSE2990", "GSE7390", "GSE11121", "GSE20685", "GSE25066"]}

# ---- build per-(cohort, stratum) analysis sets exactly as the committed path ----
data = {}
for gse in ["GSE2990", "GSE7390", "GSE11121", "GSE20685", "GSE25066"]:
    names, Xg, mapped = cohort_expr(gse, GPL[gse])
    clin = load_clin(gse)
    for stratum, ep in endpoints(gse, clin).items():
        d = ep.dropna()
        common = [s for s in names if s in d.index]
        idx = [names.index(s) for s in common]
        if len(common) < 30:
            continue
        Xs = Xg[idx]
        Xs = (Xs - Xs.mean(0)) / (Xs.std(0) + 1e-9)
        tv = d.loc[common, "t"].values.astype(float)
        ev = d.loc[common, "e"].values.astype(int)
        data[(gse, stratum)] = {"X": Xs, "t": tv, "e": ev, "n": int(len(common)),
                                "events": int(ev.sum()), "mapped": int(mapped),
                                "eta_frozen": Xs @ coef}

# ---- integrity gate: reproduce the committed frozen record before any LOCO fit ----
committed = json.load(open("results/endpoint_harmonization.json"))
for r in committed["per_cohort"]:
    g = data[(r["cohort"], r["stratum"])]
    assert g["n"] == r["n"], (r["cohort"], r["stratum"], "n")
    assert g["events"] == r["events"], (r["cohort"], r["stratum"], "events")
    assert g["mapped"] == r["genes_mapped"], (r["cohort"], r["stratum"], "mapped")
    c = cidx(g["t"], -g["eta_frozen"], g["e"])
    assert abs(c - r["c_index"]) < 5e-5, (r["cohort"], r["stratum"], c, r["c_index"])
print("integrity gate passed: frozen record reproduced", flush=True)

# ---- LOCO retrains + paired bootstrap deltas (single rng stream, seed 0) ----
rng = np.random.default_rng(0)
per = []
for stratum, cohorts in STRATA.items():
    for held in cohorts:
        train = [c for c in cohorts if c != held]
        Xtr = np.vstack([data[(c, stratum)]["X"] for c in train])
        ttr = np.concatenate([data[(c, stratum)]["t"] for c in train])
        etr = np.concatenate([data[(c, stratum)]["e"] for c in train])
        keep = Xtr.var(0) > 0
        df = {f"x{j}": Xtr[:, j] for j in np.where(keep)[0]}
        df["t"], df["e"] = ttr, etr
        f = CoxPHFitter(penalizer=0.05).fit(pd.DataFrame(df), "t", "e")
        b = f.summary["coef"].values
        g = data[(held, stratum)]
        eta_loco = g["X"][:, keep] @ b
        g["eta_loco"] = eta_loco
        c_loco = float(cidx(g["t"], -eta_loco, g["e"]))
        c_froz = float(cidx(g["t"], -g["eta_frozen"], g["e"]))
        dlt, bl, bf = [], [], []
        n = len(g["t"])
        for _ in range(1000):
            bi = rng.integers(0, n, n)
            cl = cidx(g["t"][bi], -eta_loco[bi], g["e"][bi])
            cf = cidx(g["t"][bi], -g["eta_frozen"][bi], g["e"][bi])
            if not (np.isnan(cl) or np.isnan(cf)):
                dlt.append(cl - cf); bl.append(cl); bf.append(cf)
        dlo, dhi = np.percentile(dlt, [2.5, 97.5])
        llo, lhi = np.percentile(bl, [2.5, 97.5])
        row = {"cohort": held, "stratum": stratum, "n": g["n"], "events": g["events"],
               "genes_mapped_heldout": g["mapped"], "genes_fit": int(keep.sum()),
               "train_cohorts": train, "train_n": int(len(ttr)), "train_events": int(etr.sum()),
               "c_loco": round(c_loco, 4), "c_loco_ci95": [round(float(llo), 4), round(float(lhi), 4)],
               "c_frozen": round(c_froz, 4),
               "delta": round(c_loco - c_froz, 4),
               "delta_ci95": [round(float(dlo), 4), round(float(dhi), 4)]}
        per.append(row)
        print(stratum, held, "LOCO", row["c_loco"], "frozen", row["c_frozen"],
              "delta", row["delta"], row["delta_ci95"], flush=True)

# ---- secondary descriptive: per-stratum n-weighted mean LOCO C, stratified bootstrap ----
pooled = {}
for stratum, cohorts in STRATA.items():
    rows = [r for r in per if r["stratum"] == stratum]
    ws = np.array([r["n"] for r in rows], float)
    cs = np.array([r["c_loco"] for r in rows])
    point = float((ws * cs).sum() / ws.sum())
    boots = []
    for _ in range(500):
        bcs, bws = [], []
        for r in rows:
            g = data[(r["cohort"], stratum)]
            n = len(g["t"])
            bi = rng.integers(0, n, n)
            bc = cidx(g["t"][bi], -g["eta_loco"][bi], g["e"][bi])
            if not np.isnan(bc):
                bcs.append(bc); bws.append(r["n"])
        if bcs:
            boots.append(float((np.array(bws) * np.array(bcs)).sum() / np.array(bws).sum()))
    lo, hi = np.percentile(boots, [2.5, 97.5])
    pooled[stratum] = {"pooled_loco_c_index": round(point, 4),
                       "ci95": [round(float(lo), 4), round(float(hi), 4)],
                       "definition": "n-weighted mean of per-cohort LOCO C-index; cohort-stratified patient bootstrap (500, same rng stream); descriptive, no gate"}
    print("pooled LOCO", stratum, pooled[stratum]["pooled_loco_c_index"], pooled[stratum]["ci95"], flush=True)
out = {"prereg": "docs/PREREG_LOCO_RETRAINING_20260928.md",
       "integrity": "frozen per-cohort record reproduced to 4 decimals before any LOCO fit",
       "per_cohort": per, "pooled_loco_point": pooled}
json.dump(out, open("results/loco_retraining.json", "w"), indent=1)
print("wrote results/loco_retraining.json", flush=True)
