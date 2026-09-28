"""Lane 21 item 5 pivot A: locked OS-endpoint transport of the frozen METABRIC
panel in the designated cohorts GSE7390 + GSE20685.
Per docs/PREREG_OS_TRANSPORT_20260928.md (48b72e4) and
docs/OS_TRANSPORT_FETCH_LOG_20260928.md (35c8d96), both pre-outcome.
Frozen model results/risk_model.json (87 features). No refitting, no
recalibration, no feature changes, no threshold tuning. Endpoint: each
cohort's NATIVE OS fields, verbatim keys. Field mappings and imputations
exactly as declared in the fetch log; the GSE7390 size conversion (cm x10 to
mm) is the declared one. Decision rule (locked): SUPPORTED only if the
bootstrap CI excludes 0.5 AND the point estimate beats the in-cohort
clinical-only comparator in BOTH cohorts.
"""
import gzip, json, hashlib
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index as cidx

GEO = "data_cache/external/geo"
MODEL = "results/risk_model.json"
OUT = "results/os_transport.json"

COHORTS = {
    "GSE7390": {
        "soft": f"{GEO}/GSE7390_family.soft.gz", "gpl": f"{GEO}/GPL96.txt",
        "t_key": "t.os", "e_key": "e.os", "time_unit": "days (as stored)",
        "model_map": {  # feature -> (characteristics key, transform)
            "AGE_AT_DIAGNOSIS": ("age", lambda v: v),
            "TUMOR_SIZE": ("size", lambda v: v * 10.0),  # cm -> mm (declared)
            "GRADE": ("grade", lambda v: v),
            "LYMPH_NODES_EXAMINED_POSITIVE": ("node", lambda v: v),
            "ER_STATUS=Positive": ("er", lambda v: 1.0 if v == 1.0 else 0.0),
            "NPI": ("NPI", lambda v: v),
        },
        "model_const": {"HORMONE_THERAPY=YES": 0.0, "CHEMOTHERAPY=YES": 0.0},
        "comparator": ["age", "size", "grade", "er"],
    },
    "GSE20685": {
        "soft": f"{GEO}/GSE20685_family.soft.gz", "gpl": f"{GEO}/GPL570.txt",
        "t_key": "follow_up_duration (years)", "e_key": "event_death",
        "time_unit": "years (as stored)",
        "model_map": {
            "AGE_AT_DIAGNOSIS": ("age at diagnosis", lambda v: v),
            "CHEMOTHERAPY=YES": ("adjuvant_chemotherapy",
                                 lambda v: 1.0 if str(v).lower() == "yes" else 0.0),
        },
        "model_const": {},
        "comparator": ["age at diagnosis", "t_stage", "n_stage"],
    },
}

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
                    p2g[p[i_id]] = p[i_gs].split(" ///")[0].strip()
    return p2g

def stream_soft(path, p2g, want):
    """Per-sample: characteristics dict + mean expression per wanted gene."""
    samples = []
    cur = None; mode = None
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("^SAMPLE"):
                cur = {"gsm": line.split("=")[1].strip(), "chars": {}, "acc": {}}
                samples.append(cur); mode = None
            elif cur is None:
                continue
            elif line.startswith("!Sample_characteristics_ch1"):
                v = line.split("=", 1)[1].strip()
                k, sep, val = v.partition(":")
                if sep:
                    cur["chars"][k.strip()] = val.strip()
            elif line.startswith("!sample_table_begin"):
                mode = "table"
            elif line.startswith("!sample_table_end"):
                mode = None
            elif mode == "table":
                parts = line.split("\t")
                if len(parts) < 2 or parts[0] == "ID_REF":
                    continue
                g = p2g.get(parts[0])
                if not g or g not in want:
                    continue
                try:
                    v = float(parts[1])
                except ValueError:
                    continue
                s, c = cur["acc"].get(g, (0.0, 0))
                cur["acc"][g] = (s + v, c + 1)
    return samples

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
    sha = hashlib.sha256(open(MODEL, "rb").read()).hexdigest()
    out = {"prereg": "docs/PREREG_OS_TRANSPORT_20260928.md",
           "fetch_log": "docs/OS_TRANSPORT_FETCH_LOG_20260928.md",
           "model_sha256": sha, "cohorts": {}}
    for name, cfg in COHORTS.items():
        p2g = load_p2g(cfg["gpl"])
        samples = stream_soft(cfg["soft"], p2g, set(genes))
        keep, t_l, e_l, X_l, comp_l = [], [], [], [], []
        for i, s in enumerate(samples):
            t = _f(s["chars"].get(cfg["t_key"])); e = _f(s["chars"].get(cfg["e_key"]))
            if not (np.isfinite(t) and t > 0 and e in (0, 1)):
                continue
            keep.append(i); t_l.append(t); e_l.append(e)
            X_l.append([s["acc"].get(g, (np.nan, 0))[0] / s["acc"][g][1]
                        if g in s["acc"] and s["acc"][g][1] else np.nan
                        for g in genes])
            comp_l.append([_f(s["chars"].get(k)) for k in cfg["comparator"]])
        t = np.array(t_l, float); e = np.array(e_l, float)
        Xg = np.array(X_l, float); C = np.array(comp_l, float)
        mapped = int(np.isfinite(Xg).any(axis=0).sum())
        print(f"{name}: n={len(t)} events={int(e.sum())} genes_mapped={mapped}/70", flush=True)
        assert mapped >= 60, "below locked 60-gene eligibility minimum"
        assert int(e.sum()) >= 30, "below locked 30-death scouting gate"

        F = np.full((len(t), len(feats)), np.nan)
        for feat, (key, tf) in cfg["model_map"].items():
            for r, i in enumerate(keep):
                v = _f(samples[i]["chars"].get(key))
                if np.isfinite(v):
                    F[r, feats.index(feat)] = tf(v)
        for feat, val in cfg["model_const"].items():
            F[:, feats.index(feat)] = val
        F[:, 17:] = Xg
        for j in range(len(feats)):
            nm = np.isnan(F[:, j])
            F[nm, j] = means[j]
        lp = ((F - means) / stds) @ coef

        c_point = float(cidx(t, -lp, e))
        rng = np.random.default_rng(0)
        n = len(t); boots = np.empty(1000)
        for b in range(1000):
            bi = rng.integers(0, n, n)
            boots[b] = cidx(t[bi], -lp[bi], e[bi]) if e[bi].sum() >= 2 else np.nan
        ci = [float(v) for v in np.nanpercentile(boots, [2.5, 97.5])]

        mcomp = np.isfinite(C).all(axis=1)
        cph = CoxPHFitter().fit(
            pd.DataFrame({"t": t[mcomp], "e": e[mcomp],
                          **{f"c{k}": C[mcomp, k] for k in range(C.shape[1])}}), "t", "e")
        lp_c = np.log(cph.predict_partial_hazard(
            pd.DataFrame({f"c{k}": C[mcomp, k] for k in range(C.shape[1])})).to_numpy())
        c_comp = float(cidx(t[mcomp], -lp_c, e[mcomp]))
        # frozen score restricted to the comparator's complete-case set, for a
        # like-for-like comparison set (declared descriptive; gate uses c_point)
        c_frozen_cc = float(cidx(t[mcomp], -lp[mcomp], e[mcomp]))

        out["cohorts"][name] = {
            "endpoint": f"native OS ({cfg['t_key']} / {cfg['e_key']}; {cfg['time_unit']})",
            "n": int(len(t)), "deaths": int(e.sum()), "genes_mapped": mapped,
            "c_index_frozen": c_point, "c_index_frozen_comparator_cases": c_frozen_cc,
            "bootstrap_B1000_seed0_ci95": ci,
            "comparator": {"covariates": cfg["comparator"], "n": int(mcomp.sum()),
                           "c_index_in_sample": c_comp},
            "rule_met": bool(ci[0] > 0.5 and c_point > c_comp),
        }
        print(name, json.dumps(out["cohorts"][name], indent=1), flush=True)
    out["brier"] = "NOT COMPUTED - frozen model carries no baseline hazard (prereg)"
    out["decision_rule"] = "SUPPORTED only if rule_met in BOTH cohorts"
    out["os_transport_supported"] = bool(all(c["rule_met"] for c in out["cohorts"].values()))
    json.dump(out, open(OUT, "w"), indent=1)
    print("SUPPORTED:", out["os_transport_supported"], flush=True)

if __name__ == "__main__":
    main()
