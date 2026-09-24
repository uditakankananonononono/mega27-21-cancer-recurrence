"""Assemble the survival dataset from cached cBioPortal payloads. Pure code."""
from __future__ import annotations

import json
import os
from dataclasses import dataclass

import numpy as np

from . import cbioportal as cb


@dataclass
class SurvivalDataset:
    name: str
    X_expr: np.ndarray        # (n, g) float32 expression z-scores
    X_clin: np.ndarray        # (n, c) float32 encoded clinical covariates
    time: np.ndarray          # (n,) float32 months
    event: np.ndarray         # (n,) int64 1=recurred, 0=censored
    gene_names: list
    clin_names: list
    sample_ids: list

    def __post_init__(self):
        n = len(self.time)
        assert self.X_expr.shape[0] == n == self.X_clin.shape[0]
        assert len(self.event) == n
        assert (self.time > 0).all()


NUMERIC_CLIN = ["AGE_AT_DIAGNOSIS", "TUMOR_SIZE", "GRADE",
                "LYMPH_NODES_EXAMINED_POSITIVE", "NPI"]
CATEGORICAL_CLIN = {
    "ER_STATUS": ["Positive"],
    "PR_STATUS": ["Positive"],
    "HER2_STATUS": ["Positive"],
    "CLAUDIN_SUBTYPE": ["Basal", "Her2", "LumA", "LumB", "Normal", "claudin-low"],
    "CHEMOTHERAPY": ["YES"],
    "HORMONE_THERAPY": ["YES"],
    "RADIO_THERAPY": ["YES"],
}


def encode_clinical(records: list) -> tuple:
    """One record dict per sample -> (X, names), median-imputed numerics."""
    names = list(NUMERIC_CLIN)
    for attr, levels in CATEGORICAL_CLIN.items():
        names.extend(f"{attr}={lvl}" for lvl in levels)
    X = np.full((len(records), len(names)), np.nan, dtype=np.float32)
    for i, rec in enumerate(records):
        for j, attr in enumerate(NUMERIC_CLIN):
            try:
                X[i, j] = float(rec.get(attr, ""))
            except (TypeError, ValueError):
                pass
        k = len(NUMERIC_CLIN)
        for attr, levels in CATEGORICAL_CLIN.items():
            val = str(rec.get(attr, "")).strip()
            for lvl in levels:
                X[i, k] = 1.0 if val == lvl else 0.0
                k += 1
    # median impute numerics only
    for j in range(len(NUMERIC_CLIN)):
        col = X[:, j]
        mask = ~np.isnan(col)
        X[~mask, j] = np.median(col[mask]) if mask.any() else 0.0
    return X, names


def assemble(cache_dir: str = None) -> SurvivalDataset:
    """Join expression panel + clinical on samples with a valid RFS endpoint."""
    cache_dir = cache_dir or cb.CACHE_DIR

    def load(name):
        with open(os.path.join(cache_dir, name)) as fh:
            return json.load(fh)

    genes = load("genes.json")
    clinical = load("clinical.json")
    expr_rows = load("expression.json")
    entrez_to_symbol = {v: k for k, v in genes.items()}
    gene_list = sorted(entrez_to_symbol)
    gcol = {g: j for j, g in enumerate(gene_list)}

    expr_by_sample = {}
    for row in expr_rows:
        g = row.get("entrezGeneId")
        if g in gcol:
            expr_by_sample.setdefault(row["sampleId"], {})[gcol[g]] = row["value"]

    sample_ids, times, events = [], [], []
    for sid, rec in clinical.items():
        status = str(rec.get("RFS_STATUS", ""))
        if not (status.startswith("0") or status.startswith("1")):
            continue
        try:
            t = float(rec.get("RFS_MONTHS", ""))
        except (TypeError, ValueError):
            continue
        if t <= 0 or sid not in expr_by_sample:
            continue
        sample_ids.append(sid)
        times.append(t)
        events.append(int(status[0]))

    n, g = len(sample_ids), len(gene_list)
    X_expr = np.zeros((n, g), dtype=np.float32)
    for i, sid in enumerate(sample_ids):
        for j, v in expr_by_sample[sid].items():
            X_expr[i, j] = v
    records = [clinical[sid] for sid in sample_ids]
    X_clin, clin_names = encode_clinical(records)

    return SurvivalDataset(
        name="metabric_brca_rfs",
        X_expr=X_expr, X_clin=X_clin,
        time=np.array(times, dtype=np.float32),
        event=np.array(events, dtype=np.int64),
        gene_names=[entrez_to_symbol[g] for g in gene_list],
        clin_names=clin_names,
        sample_ids=sample_ids,
    )


def make_synthetic_survival(n: int = 400, d: int = 10, seed: int = 0,
                            effect: float = 1.5, censor_rate: float = 0.3):
    """Weibull survival with log-linear hazards; administrative censoring."""
    rng = np.random.default_rng(seed)
    X = rng.normal(size=(n, d)).astype(np.float32)
    beta = np.zeros(d, dtype=np.float32)
    beta[:3] = effect
    eta = X @ beta
    scale = np.exp(-eta)  # higher risk -> shorter survival
    t_event = rng.weibull(1.5, n) * scale * 10.0
    t_censor = rng.weibull(1.5, n) * np.quantile(t_event, 1 - censor_rate + 0.05) * 2
    time = np.minimum(t_event, t_censor).astype(np.float32) + 1e-3
    event = (t_event <= t_censor).astype(np.int64)
    return X, time, event, beta
