"""Cox math verified against lifelines and hand computations."""
import numpy as np
import torch
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index as ll_ci
import pandas as pd

from recurscan.data.dataset import make_synthetic_survival
from recurscan.survival import (breslow_baseline_hazard,
                                concordance_index,
                                concordance_index_fast,
                                cox_partial_loglik)


def test_cindex_matches_lifelines():
    rng = np.random.default_rng(0)
    X, t, e, _ = make_synthetic_survival(n=200, seed=1)
    risk = X @ np.array([1.5, 1.5, 1.5] + [0] * 7)
    ours = concordance_index(t, e, risk)
    ref = ll_ci(t, -risk, e)  # lifelines: higher score = longer survival
    assert abs(ours - ref) < 1e-9


def test_cindex_fast_matches_slow():
    X, t, e, _ = make_synthetic_survival(n=150, seed=2)
    risk = X[:, 0]
    assert abs(concordance_index(t, e, risk)
               - concordance_index_fast(t, e, risk)) < 1e-12


def test_cindex_hand_computed():
    # i=0 events at t=1, comparable to j=1,2; risk says 0 dies first
    t = np.array([1.0, 2.0, 3.0])
    e = np.array([1, 1, 0])
    risk = np.array([3.0, 2.0, 1.0])
    assert concordance_index(t, e, risk) == 1.0
    assert concordance_index(t, e, risk[::-1]) == 0.0


def test_cox_loglik_gradient_descent_recovers_betas():
    X, t, e, beta = make_synthetic_survival(n=400, d=5, seed=3, effect=2.0)
    Xt = torch.tensor(X)
    tt = torch.tensor(t)
    et = torch.tensor(e, dtype=torch.float32)
    w = torch.zeros(5, requires_grad=True)
    opt = torch.optim.Adam([w], lr=0.05)
    for _ in range(500):
        loss = -cox_partial_loglik(Xt @ w, tt, et)
        opt.zero_grad(); loss.backward(); opt.step()
    est = w.detach().numpy()
    # estimated coefficients align with truth direction and beat chance C
    c = concordance_index_fast(t, e, X @ est)
    assert c > 0.8
    assert np.corrcoef(est, beta)[0, 1] > 0.9


def test_cox_loglik_matches_lifelines_fit():
    X, t, e, beta = make_synthetic_survival(n=300, d=4, seed=4, effect=1.2)
    df = pd.DataFrame(X, columns=[f"x{i}" for i in range(4)])
    df["time"], df["event"] = t, e
    cph = CoxPHFitter()
    cph.fit(df, "time", "event")
    eta_ll = df[[f"x{i}" for i in range(4)]].values @ cph.params_.values
    ours = cox_partial_loglik(torch.tensor(eta_ll), torch.tensor(t),
                              torch.tensor(e, dtype=torch.float32)).item()
    # lifelines reports the same Breslow partial log-lik (unnormalized)
    assert abs(ours * e.sum() - cph.log_likelihood_) < 1e-3


def test_breslow_baseline_monotone():
    X, t, e, _ = make_synthetic_survival(n=200, seed=5)
    times, H = breslow_baseline_hazard(t, e, X[:, 0])
    assert (np.diff(H) >= 0).all()
    assert len(times) == e.sum()
