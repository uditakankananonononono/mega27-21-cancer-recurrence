"""Lane 21 item 5 follow-up: locked prospective transport evaluation on GSE17705.
Per docs/PREREG_TRANSPORT_GSE17705_20260928.md (committed 5ece112, pre-outcome)
and docs/TRANSPORT_FRESH_FETCH_LOG_20260928.md. Frozen model
results/risk_model.json (SHA-256 recorded in the execution commit). No refitting,
no recalibration, no feature changes, no threshold tuning. Endpoint: the cohort's
NATIVE DRFS fields ("distant relapse (1=dr, 0 censored)" / "event time (years)"),
matching the prereg's RFS/DRFS/DMFS wording; no DFS substitution.

Declared field handling (locked in this script, reported in paper/queue):
- ER_STATUS=Positive = 1.0 for all samples and HORMONE_THERAPY=YES = 1.0 for all
  samples: both are series-level RECORDED fields (all 298 are ER+, all uniformly
  tamoxifen-treated 5 years) - set, not imputed.
- AGE_AT_DIAGNOSIS, GRADE, TUMOR_SIZE, NPI, PR_STATUS=Positive,
  HER2_STATUS=Positive, the six CLAUDIN_SUBTYPE dummies, CHEMOTHERAPY=YES and
  RADIO_THERAPY=YES are absent from the series matrix: set to the training mean
  (standardized contribution 0).
- LYMPH_NODES_EXAMINED_POSITIVE (a count) is NOT mappable from the cohort's
  binary nodal status: imputed to training mean. The binary nodal field is used
  ONLY for the clinical-only comparator, as locked in the prereg.
- Comparator (locked): Cox on the binary nodal field alone (grade/age absent;
  ER constant), 287 samples with recorded nodal status.
- Descriptive only, no gate: per-profiling-lab C-index (MDA vs JBI).
"""
import gzip, json
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index as cidx

MATRIX = "data_cache/geo/GSE17705_series_matrix.txt.gz"
GPL = "data_cache/external/geo/GPL96.txt"
MODEL = "results/risk_model.json"
OUT = "results/transport_gse17705.json"

EVT_KEY = "distant relapse (1=dr, 0 censored)"
TIME_KEY = "event time (years)"
NODE_KEY = "nodal status (0=negative, 1=positive, na=not applicable)"
LAB_KEY = "profiling lab"

def parse_matrix(path):
    samples = None
    clin_rows = []
    expr = []
    with gzip.open(path, "rt") as fh:
        in_table = False
        for ln in fh:
            if ln.startswith("!Sample_geo_accession"):
                samples = [s.strip('"') for s in ln.rstrip("\n").split("\t")[1:]]
            elif ln.startswith("!Sample_characteristics_ch1"):
                clin_rows.append(ln)
            elif ln.startswith("!series_matrix_table_begin"):
                in_table = True
                fh.readline()
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

def _f(v):
    try:
        return float(str(v).strip())
    except (TypeError, ValueError):
        return np.nan

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

    # endpoint + clinical alignment
    keep, t_all, e_all, node_all, lab_all = [], [], [], [], []
    for i, s in enumerate(samples):
        c = clin.get(s, {})
        ev = _f(c.get(EVT_KEY)); tm = _f(c.get(TIME_KEY))
        if np.isnan(ev) or np.isnan(tm):
            continue
        keep.append(i); t_all.append(tm); e_all.append(ev)
        node_all.append(_f(c.get(NODE_KEY)))
        lab_all.append(c.get(LAB_KEY, "?"))
    t = np.array(t_all, float); e = np.array(e_all, float)
    node = np.array(node_all, float); lab = np.array(lab_all)
    n_events = int(e.sum())
    print(f"samples with endpoint: {len(t)}/{len(samples)}, events: {n_events}", flush=True)
    assert n_events >= 40, "below locked 40-event eligibility minimum"

    F = np.full((len(t), len(feats)), np.nan)
    # recorded series-level constants (see module docstring)
    F[:, feats.index("ER_STATUS=Positive")] = 1.0
    F[:, feats.index("HORMONE_THERAPY=YES")] = 1.0
    # expression
    F[:, 17:] = Xg[keep]
    # everything else clinical -> training mean (declared above)
    for j in range(len(feats)):
        nanmask = np.isnan(F[:, j])
        F[nanmask, j] = means[j]
    Xs = (F - means) / stds
    lp = Xs @ coef

    c_point = float(cidx(t, -lp, e))
    rng = np.random.default_rng(0)
    boots = np.empty(1000)
    n = len(t)
    for b in range(1000):
        bi = rng.integers(0, n, n)
        if e[bi].sum() < 2:
            boots[b] = np.nan
            continue
        boots[b] = cidx(t[bi], -lp[bi], e[bi])
    ci = [float(v) for v in np.nanpercentile(boots, [2.5, 97.5])]

    cdf = pd.DataFrame({"t": t, "e": e, "lp": lp})
    cph = CoxPHFitter().fit(cdf, "t", "e")
    slope = float(cph.params_["lp"])
    slope_ci = [float(v) for v in cph.confidence_intervals_.loc["lp"].tolist()]

    # clinical-only comparator (locked: binary nodal only)
    m = ~np.isnan(node)
    comp = pd.DataFrame({"t": t[m], "e": e[m], "node": node[m]})
    cph_c = CoxPHFitter().fit(comp, "t", "e")
    lp_c = cph_c.predict_partial_hazard(comp).to_numpy()
    c_comp = float(cidx(comp["t"].to_numpy(), -np.log(lp_c), comp["e"].to_numpy()))

    # descriptive per-lab breakdown (no gate)
    per_lab = {}
    for L in sorted(set(lab)):
        mL = lab == L
        if mL.sum() >= 20 and e[mL].sum() >= 2:
            per_lab[L] = {"n": int(mL.sum()), "events": int(e[mL].sum()),
                          "c_index": float(cidx(t[mL], -lp[mL], e[mL]))}

    verdict = bool(ci[0] > 0.5 and c_point > c_comp)
    out = {
        "cohort": "GSE17705", "platform": "GPL96",
        "endpoint": "DRFS (native fields 'distant relapse (1=dr, 0 censored)' / 'event time (years)')",
        "n_samples_matrix": len(samples), "n_endpoint": int(len(t)),
        "n_events": n_events,
        "genes_mapped": mapped,
        "set_from_series_records": ["ER_STATUS=Positive=1.0 (all ER+)",
                                    "HORMONE_THERAPY=YES=1.0 (uniform tamoxifen 5y)"],
        "imputed_to_training_mean": ["AGE_AT_DIAGNOSIS", "GRADE", "TUMOR_SIZE", "NPI",
                                     "PR_STATUS=Positive", "HER2_STATUS=Positive",
                                     "CLAUDIN_SUBTYPE dummies (6)", "CHEMOTHERAPY=YES",
                                     "RADIO_THERAPY=YES", "LYMPH_NODES_EXAMINED_POSITIVE"],
        "c_index_frozen": c_point,
        "bootstrap_B1000_seed0_ci95": ci,
        "calibration_slope": slope, "calibration_slope_ci": slope_ci,
        "comparator": {"covariates": ["binary nodal status"], "n": int(len(comp)),
                       "c_index_in_sample": c_comp},
        "per_profiling_lab_descriptive": per_lab,
        "brier": "NOT COMPUTED - frozen model carries no baseline hazard (declared conditional in prereg)",
        "decision_rule": "supported if CI excludes 0.5 AND point estimate > clinical-only comparator",
        "transport_supported": verdict,
    }
    json.dump(out, open(OUT, "w"), indent=1)
    print(json.dumps(out, indent=1), flush=True)

if __name__ == "__main__":
    main()
