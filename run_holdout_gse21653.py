"""Lane 21 item 5: locked prospective holdout evaluation on GSE21653.
Per docs/PREREG_PROSPECTIVE_HOLDOUT_20260927.md and
docs/HOLDOUT_FETCH_LOG_20260928.md. Frozen model results/risk_model.json
(SHA-256 8ae3992542e60468...). No refitting, no recalibration, no feature
changes, no threshold tuning.

Declared imputation (locked in this script, reported in paper/queue):
TUMOR_SIZE, LYMPH_NODES_EXAMINED_POSITIVE, NPI and the three treatment flags
are unavailable as numeric/recorded fields in the GSE21653 series matrix
(only T-stage and binary nodal status exist); those features are set to the
training mean (standardized contribution 0). Per-sample missing IHC values
are likewise imputed to the training mean. Comparator uses covariates
available in this cohort, as the prereg specifies.
"""
import gzip, json
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index as cidx

MATRIX = "data_cache/geo/GSE21653_series_matrix.txt.gz"
GPL = "data_cache/external/geo/GPL570.txt"
MODEL = "results/risk_model.json"
OUT = "results/holdout_gse21653.json"

def parse_matrix(path):
    """Stream the series matrix: samples, per-sample clinical dict, expression."""
    samples = None
    clin_rows = []
    expr = []
    with gzip.open(path, "rt") as fh:
        in_table = False
        n_cols = 0
        for ln in fh:
            if ln.startswith("!Sample_geo_accession"):
                samples = [s.strip('"') for s in ln.rstrip("\n").split("\t")[1:]]
            elif ln.startswith("!Sample_characteristics_ch1"):
                clin_rows.append(ln)
            elif ln.startswith("!series_matrix_table_begin"):
                in_table = True
                header = fh.readline()
                n_cols = len(header.split("\t")) - 1
            elif ln.startswith("!series_matrix_table_end"):
                in_table = False
            elif in_table:
                parts = ln.rstrip("\n").split("\t")
                parts[0] = parts[0].strip('"')
                expr.append(parts)
    ids = [r[0] for r in expr]
    X = np.empty((len(expr), len(samples)), dtype=np.float32)
    for i, r in enumerate(expr):
        for j, v in enumerate(r[1:]):
            try:
                X[i, j] = float(v)
            except ValueError:
                X[i, j] = np.nan
    del expr
    clin = {}
    for row in clin_rows:
        cells = [c.strip('"') for c in row.rstrip("\n").split("\t")[1:]]
        for i, cell in enumerate(cells):
            if i >= len(samples) or ": " not in cell:
                continue
            k, v = cell.split(": ", 1)
            clin.setdefault(samples[i], {})[k.strip().lower()] = v.strip()
    return samples, clin, ids, X


def _f(v):
    """float or nan; GEO clinical strings use NA/empty."""
    try:
        return float(str(v).strip())
    except (TypeError, ValueError):
        return np.nan

def load_p2g(path):
    p2g = {}
    with open(path) as fh:
        intable = False
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("!platform_table_begin"):
                intable = True
                h = fh.readline().rstrip("\n").split("\t")
                i_id, i_gs = h.index("ID"), h.index("Gene Symbol")
                continue
            if line.startswith("!platform_table_end"):
                break
            if intable:
                p = line.split("\t")
                if len(p) > max(i_id, i_gs):
                    p2g[p[i_id]] = p[i_gs]
    return p2g

def main():
    model = json.load(open(MODEL))
    feats = model["features"]
    genes = feats[17:]
    coef = np.array([model["coef"][f"x{i}"] for i in range(len(feats))])
    means = np.array([model["means"][f"x{i}"] if isinstance(model["means"], dict)
                      else model["means"][i] for i in range(len(feats))])
    stds = np.array([model["stds"][f"x{i}"] if isinstance(model["stds"], dict)
                     else model["stds"][i] for i in range(len(feats))])
    samples, clin, ids, X = parse_matrix(MATRIX)
    p2g = load_p2g(GPL)
    gsym = np.array([p2g.get(i, "") for i in ids])
    Xg = np.full((len(samples), len(genes)), np.nan)
    mapped = 0
    for j, g in enumerate(genes):
        m = gsym == g
        if m.sum():
            Xg[:, j] = np.nanmean(X[m], axis=0)
            mapped += 1
    print(f"genes mapped: {mapped}/{len(genes)}", flush=True)
    assert mapped >= 60, "below locked 60-gene eligibility minimum"

    rows = []
    IMPUTE = {"TUMOR_SIZE", "LYMPH_NODES_EXAMINED_POSITIVE", "NPI",
              "CHEMOTHERAPY=YES", "HORMONE_THERAPY=YES", "RADIO_THERAPY=YES"}
    SUB = {"LuminalA": "CLAUDIN_SUBTYPE=LumA", "LuminalB": "CLAUDIN_SUBTYPE=LumB",
           "Basal": "CLAUDIN_SUBTYPE=Basal", "ERBB2": "CLAUDIN_SUBTYPE=Her2",
           "Normal": "CLAUDIN_SUBTYPE=Normal", "claudin-low": "CLAUDIN_SUBTYPE=claudin-low"}
    n_imputed_ihc = 0
    for i, s in enumerate(samples):
        c = clin.get(s, {})
        def fget(k):
            v = c.get(k, "")
            try:
                return float(v)
            except ValueError:
                return np.nan
        evt = fget("dfs evt"); tim = fget("dfs time (months)")
        row = {"sample": s, "dfs_evt": evt, "dfs_time": tim}
        vals = {}
        vals["AGE_AT_DIAGNOSIS"] = fget("age at diagnosis")
        vals["GRADE"] = fget("sbr grade")
        for feat, key in (("ER_STATUS=Positive", "er ihc"),
                          ("PR_STATUS=Positive", "pr ihc"),
                          ("HER2_STATUS=Positive", "erbb2 ihc")):
            v = fget(key)
            if np.isnan(v):
                n_imputed_ihc += 1
            vals[feat] = v
        sub = c.get("molecular subtype", "")
        for f in SUB.values():
            vals[f] = np.nan
        if sub in SUB:
            for f in SUB.values():
                vals[f] = 1.0 if f == SUB[sub] else 0.0
        row["vals"] = vals
        rows.append(row)
    df = pd.DataFrame(rows)
    n_endpoint = int((~df["dfs_evt"].isna() & ~df["dfs_time"].isna()).sum())
    df = df[~df["dfs_evt"].isna() & ~df["dfs_time"].isna()].reset_index(drop=True)
    print(f"samples with endpoint: {len(df)}/{len(samples)}", flush=True)

    F = np.zeros((len(df), len(feats)))
    for j, feat in enumerate(feats[:17]):
        F[:, j] = np.nan
        for i in range(len(df)):
            v = df.loc[i, "vals"].get(feat, np.nan)
            F[i, j] = v
    # expression matrix aligned to df sample order
    idx = [samples.index(s) for s in df["sample"]]
    Xg_df = Xg[idx]
    F[:, 17:] = Xg_df
    # impute: clinical unavailable fields and any NaN -> training mean
    for j, feat in enumerate(feats):
        if feat in IMPUTE:
            F[:, j] = means[j]
        else:
            nanmask = np.isnan(F[:, j])
            F[nanmask, j] = means[j]
    Xs = (F - means) / stds
    lp = Xs @ coef
    t = df["dfs_time"].to_numpy(float)
    e = df["dfs_evt"].to_numpy(float)
    c_point = float(cidx(t, -lp, e))
    rng = np.random.default_rng(0)
    boots = np.empty(1000)
    n = len(df)
    for b in range(1000):
        bi = rng.integers(0, n, n)
        if e[bi].sum() < 2:
            boots[b] = np.nan
            continue
        boots[b] = cidx(t[bi], -lp[bi], e[bi])
    ci = [float(v) for v in np.nanpercentile(boots, [2.5, 97.5])]
    # calibration slope: Cox of endpoint on lp
    cdf = pd.DataFrame({"t": t, "e": e, "lp": lp})
    cph = CoxPHFitter().fit(cdf, "t", "e")
    slope = float(cph.params_["lp"])
    slope_ci = [float(v) for v in cph.confidence_intervals_.loc["lp"].tolist()]
    # clinical-only comparator (covariates available in this cohort)
    comp = pd.DataFrame({"t": t, "e": e,
                         "age": F[:, feats.index("AGE_AT_DIAGNOSIS")],
                         "grade": F[:, feats.index("GRADE")],
                         "pn": np.array([df.loc[i, "vals"].get("ER_STATUS=Positive", np.nan)
                                         for i in range(len(df))])})
    # use raw available covariates, not imputed: age, grade, binary nodal pn
    pn = np.array([_f(clin.get(s, {}).get("pn")) for s in df["sample"]])
    comp = pd.DataFrame({"t": t, "e": e,
                         "age": np.array([_f(clin.get(s, {}).get("age at diagnosis")) for s in df["sample"]]),
                         "grade": np.array([_f(clin.get(s, {}).get("sbr grade")) for s in df["sample"]]),
                         "pn": pn}).dropna()
    cph_c = CoxPHFitter().fit(comp, "t", "e")
    lp_c = cph_c.predict_partial_hazard(comp).to_numpy()
    c_comp = float(cidx(comp["t"].to_numpy(), -np.log(lp_c), comp["e"].to_numpy()))
    verdict = bool(ci[0] > 0.5 and c_point > c_comp)
    out = {
        "cohort": "GSE21653", "platform": "GPL570",
        "endpoint": "DFS (dfs evt / dfs time months), declared natively in series matrix",
        "n_samples_matrix": len(samples), "n_endpoint": len(df),
        "n_events": int(e.sum()),
        "genes_mapped": mapped,
        "imputed_to_training_mean": sorted(IMPUTE),
        "n_ihc_values_imputed": n_imputed_ihc,
        "c_index_frozen": c_point,
        "bootstrap_B1000_seed0_ci95": ci,
        "calibration_slope": slope, "calibration_slope_ci": slope_ci,
        "comparator": {"covariates": ["age", "sbr grade", "pn (binary)"],
                       "n": int(len(comp)), "c_index_in_sample": c_comp},
        "brier": "NOT COMPUTED - frozen model carries no baseline hazard, so absolute risk is unavailable (declared in prereg secondary as conditional)",
        "decision_rule": "supported if CI excludes 0.5 AND point estimate > clinical-only comparator",
        "transport_supported": verdict,
    }
    json.dump(out, open(OUT, "w"), indent=1)
    print(json.dumps(out, indent=1), flush=True)

if __name__ == "__main__":
    main()
