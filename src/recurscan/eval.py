"""Survival metrics with bootstrap CIs."""
from __future__ import annotations

import numpy as np

from .survival import concordance_index_fast


def survival_metrics(time: np.ndarray, event: np.ndarray,
                     risk: np.ndarray) -> dict:
    return {
        "n": int(len(time)),
        "events": int(np.sum(event)),
        "event_rate": float(np.mean(event)),
        "c_index": float(concordance_index_fast(time, event, risk)),
    }


def bootstrap_cindex_ci(time, event, risk, n_boot: int = 500, seed: int = 0):
    rng = np.random.default_rng(seed)
    stats = []
    n = len(time)
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        ci = concordance_index_fast(time[idx], event[idx], risk[idx])
        if not np.isnan(ci):
            stats.append(ci)
    lo, hi = np.percentile(stats, [2.5, 97.5])
    return float(lo), float(hi)


def summarize(runs: list) -> dict:
    keys = [k for k in runs[0] if isinstance(runs[0][k], float)]
    return {k: {"mean": float(np.mean([r[k] for r in runs])),
                "std": float(np.std([r[k] for r in runs]))} for k in keys}
