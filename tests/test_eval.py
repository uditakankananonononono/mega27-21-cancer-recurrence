"""Bootstrap CI helper used by transport and endpoint-harmonization runs."""
import numpy as np

from recurscan.eval import bootstrap_cindex_ci, survival_metrics, summarize
from recurscan.survival import concordance_index_fast


def test_bootstrap_ci_deterministic_and_brackets_point():
    rng = np.random.default_rng(7)
    n = 300
    risk = rng.normal(size=n)
    t = rng.exponential(10, n) + 2 * (risk < 0)
    e = rng.binomial(1, 0.6, n).astype(float)
    lo1, hi1 = bootstrap_cindex_ci(t, e, risk, n_boot=200, seed=0)
    lo2, hi2 = bootstrap_cindex_ci(t, e, risk, n_boot=200, seed=0)
    assert (lo1, hi1) == (lo2, hi2)  # fixed seed reproduces exactly
    point = concordance_index_fast(t, e, risk)
    assert lo1 <= point <= hi1
    assert 0.0 <= lo1 < hi1 <= 1.0


def test_survival_metrics_keys():
    t = np.array([1.0, 2.0, 3.0, 4.0])
    e = np.array([1.0, 1.0, 0.0, 1.0])
    r = np.array([3.0, 2.0, 1.0, 0.5])
    m = survival_metrics(t, e, r)
    assert m["n"] == 4 and m["events"] == 3
    assert 0.0 <= m["c_index"] <= 1.0


def test_summarize_mean_std():
    runs = [{"c_index": 0.6, "brier": 0.2}, {"c_index": 0.7, "brier": 0.3}]
    s = summarize(runs)
    assert abs(s["c_index"]["mean"] - 0.65) < 1e-12
    assert abs(s["brier"]["std"] - 0.05) < 1e-12
