"""recurscan CLI - usable tool built on the item-21 work.

  python -m recurscan fit [--out results/risk_model.json]   train Cox risk model (METABRIC)
  python -m recurscan risk --model M --features '{"NPI": 3.5, "AGE_AT_DIAGNOSIS": 60, ...}'
  python -m recurscan benchmark                             print verified benchmark numbers
  python -m recurscan validate-geo --gse GSE2034 --genes CEP55,CDC20 [--label bone_relapse]
"""
import argparse, json, os, sys
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np


def cmd_fit(a):
    import pandas as pd
    from lifelines import CoxPHFitter
    from recurscan.data.dataset import assemble
    from recurscan.benchmark import stratified_event_split, standardize
    from recurscan.eval import survival_metrics
    ds = assemble()
    names = ds.clin_names + ds.gene_names
    X = np.hstack([ds.X_clin, ds.X_expr]).astype(float)
    tr, te = stratified_event_split(ds.event, 0.2, seed=0)
    Xtr, Xte = standardize(X[tr], X[te])
    df = {f"x{j}": Xtr[:, j] for j in range(X.shape[1])}
    df["t"], df["e"] = ds.time[tr], ds.event[tr]
    fit = CoxPHFitter(penalizer=0.05).fit(pd.DataFrame(df), "t", "e")
    risk = fit.predict_partial_hazard(pd.DataFrame(
        {f"x{j}": Xte[:, j] for j in range(X.shape[1])})).values.ravel()
    m = survival_metrics(ds.time[te], ds.event[te], risk)
    # cohort risk terciles for calibration labels
    risk_tr = fit.predict_partial_hazard(pd.DataFrame(df).drop(columns=["t", "e"])).values.ravel()
    q = np.quantile(np.log(risk_tr), [1/3, 2/3])
    model = {"features": names, "coef": {f"x{j}": float(fit.summary.iloc[j]["coef"]) for j in range(len(names))},
             "means": X[tr].mean(0).tolist(), "stds": (X[tr].std(0) + 1e-9).tolist(),
             "tercile_cut_log": q.tolist(),
             "heldout_cindex": round(float(m["c_index"]), 4),
             "trained_on": "METABRIC n=%d (held-out n=%d)" % (len(tr), len(te))}
    json.dump(model, open(a.out, "w"), indent=1)
    print(json.dumps({"heldout_cindex": model["heldout_cindex"], "saved": a.out}, indent=1))


def cmd_risk(a):
    model = json.load(open(a.model))
    feats = json.loads(a.features)
    x = np.array([feats.get(n, np.nan) for n in model["features"]], dtype=float)
    means, stds = np.array(model["means"]), np.array(model["stds"])
    x = np.where(np.isnan(x), means, x)  # missing -> cohort mean
    z = (x - means) / stds
    logh = sum(model["coef"][f"x{j}"] * z[j] for j in range(len(z)))
    terc = int(np.digitize(logh, model["tercile_cut_log"])) + 1
    print(json.dumps({"log_hazard": round(float(logh), 3), "risk_tercile": terc,
                      "band": ["low", "intermediate", "high"][terc - 1],
                      "model": model["trained_on"],
                      "model_heldout_cindex": model["heldout_cindex"]}, indent=1))


def cmd_benchmark(a):
    d = json.load(open("results/metabric_rfs_benchmark.json"))
    s = d["summary"] if "summary" in d else d
    print(json.dumps(s, indent=1))


def cmd_validate_geo(a):
    import GEOparse
    from sklearn.metrics import roc_auc_score
    from scipy import stats
    genes = a.genes.split(",")
    g = GEOparse.get_GEO(geo=a.gse, destdir="data_cache/external/geo",
                         annotate_gpl=False, silent=True)
    gpl = GEOparse.get_GEO(geo="GPL96", destdir="data_cache/external/geo", silent=True)
    p2g = dict(zip(gpl.table["ID"], gpl.table["Gene Symbol"]))
    samples = list(g.gsms.values())
    key = [c for c in samples[0].metadata.get("characteristics_ch1", []) if "1=yes" in c][0]
    y = np.array([int(s.metadata["characteristics_ch1"][
        [i for i, c in enumerate(s.metadata["characteristics_ch1"]) if "1=yes" in c][0]
    ].split(":")[1]) for s in samples])
    gids = np.array([p2g.get(i, "") for i in samples[0].table["ID_REF"].values])
    X = np.column_stack([s.table["VALUE"].values for s in samples])
    scores = []
    for gene in genes:
        m = gids == gene
        if m.sum():
            row = X[m].mean(0)
            scores.append((row - row.mean()) / (row.std() + 1e-9))
    S = np.mean(scores, axis=0)
    auc = roc_auc_score(y, S)
    t, p = stats.ttest_ind(S[y == 1], S[y == 0])
    print(json.dumps({"gse": a.gse, "n": int(len(y)), "events": int(y.sum()),
                      "genes_mapped": len(scores), "endpoint": key.split(":")[0],
                      "sig_auc": round(float(auc), 4), "ttest_p": float(p)}, indent=1))


def main(argv=None):
    ap = argparse.ArgumentParser(prog="recurscan")
    sub = ap.add_subparsers(dest="cmd", required=True)
    f = sub.add_parser("fit"); f.add_argument("--out", default="results/risk_model.json")
    r = sub.add_parser("risk"); r.add_argument("--model", default="results/risk_model.json")
    r.add_argument("--features", required=True)
    sub.add_parser("benchmark")
    v = sub.add_parser("validate-geo"); v.add_argument("--gse", required=True)
    v.add_argument("--genes", required=True)
    a = ap.parse_args(argv)
    {"fit": cmd_fit, "risk": cmd_risk, "benchmark": cmd_benchmark,
     "validate-geo": cmd_validate_geo}[a.cmd](a)


if __name__ == "__main__":
    main()
