"""Discovery experiments: candidate recurrence signatures beyond NPI/PAM50.

Gate (locked before looking): a candidate WINS only if it adds >= +0.01
C-index over the coxph_full baseline on held-out data, bootstrap CI on the
delta excluding 0. Losers are reported as honest negatives.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter

from .benchmark import standardize, stratified_event_split
from .graphs import knn_similarity_graph, normalize_adjacency

PROLIFERATION = ["MKI67", "CCNB1", "CDC20", "CENPF", "BIRC5", "UBE2C",
                 "CCNE1", "MYC", "TYMS", "RRM2"]
ER_SIGNAL = ["ESR1", "PGR", "FOXA1", "GATA3"]
IMMUNE_STROMAL = ["STAT3", "VEGFA", "HIF1A"]


def module_scores(ds) -> dict:
    """Biology-driven module scores from the 70-gene panel."""
    g = {name: i for i, name in enumerate(ds.gene_names)}
    out = {}
    for mod, genes in (("proliferation", PROLIFERATION),
                       ("er_signaling", ER_SIGNAL),
                       ("hypoxia_immune", IMMUNE_STROMAL)):
        idx = [g[x] for x in genes if x in g]
        out[mod] = ds.X_expr[:, idx].mean(axis=1)
    # ratio feature: proliferation minus ER signaling (aggressive phenotype)
    out["prolif_minus_er"] = out["proliferation"] - out["er_signaling"]
    return out


def graph_smoothed_risk(X_feats, y_risk_train, idx_train, k=8, alpha=0.5):
    """Label-spreading of train risk scores through the patient graph.

    Returns a per-node smoothed risk; for test nodes this is transductive
    (uses only train risk values, graph structure includes all nodes).
    """
    A = normalize_adjacency(knn_similarity_graph(X_feats, k=k))
    n = X_feats.shape[0]
    F = np.zeros(n)
    F[idx_train] = y_risk_train
    F0 = F.copy()
    for _ in range(50):
        F = alpha * (A @ F) + (1 - alpha) * F0
    return F


def fit_cox(X, time, event, penalizer=0.1):
    df = pd.DataFrame(X, columns=[f"x{i}" for i in range(X.shape[1])])
    df["time"], df["event"] = time, event
    cph = CoxPHFitter(penalizer=penalizer)
    cph.fit(df, "time", "event")
    return cph


def cox_risk(cph, X):
    df = pd.DataFrame(X, columns=[f"x{i}" for i in range(X.shape[1])])
    return np.log(cph.predict_partial_hazard(df).values.ravel())


def bootstrap_delta_ci(time, event, risk_a, risk_b, n_boot=500, seed=0):
    """CI for C(risk_a) - C(risk_b) on the SAME resamples (paired)."""
    from .survival import concordance_index_fast
    rng = np.random.default_rng(seed)
    deltas = []
    n = len(time)
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        ca = concordance_index_fast(time[idx], event[idx], risk_a[idx])
        cb = concordance_index_fast(time[idx], event[idx], risk_b[idx])
        if not (np.isnan(ca) or np.isnan(cb)):
            deltas.append(ca - cb)
    return float(np.mean(deltas)), tuple(np.percentile(deltas, [2.5, 97.5]))
