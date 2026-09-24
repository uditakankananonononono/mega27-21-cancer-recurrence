"""Cox proportional-hazards math: partial likelihood, C-index, baseline hazard."""
from __future__ import annotations

import numpy as np
import torch


def cox_partial_loglik(risk: torch.Tensor, time: torch.Tensor,
                       event: torch.Tensor) -> torch.Tensor:
    """Breslow partial log-likelihood. risk = linear predictor eta (n,).

    L = sum_{i: event} [ eta_i - log( sum_{j: t_j >= t_i} exp(eta_j) ) ]
    Computed by sorting on descending time and cumulative logsumexp.
    """
    order = torch.argsort(time, descending=True)
    eta = risk[order]
    ev = event[order]
    # cumulative logsumexp over the sorted (descending time) axis gives, at
    # position i, log(sum_{j >= i in sorted order} exp(eta_j)) = risk set of i.
    # One GLOBAL max shift keeps every term on the same scale (a running max
    # would rescale earlier terms and corrupt the cumsum).
    m = eta.max()
    cum_exp = torch.cumsum(torch.exp(eta - m), dim=0)
    log_risk = m + torch.log(cum_exp)
    terms = (eta - log_risk) * ev
    n_events = ev.sum()
    return terms.sum() / torch.clamp(n_events, min=1.0)


def concordance_index(time: np.ndarray, event: np.ndarray,
                      risk: np.ndarray) -> float:
    """Harrell's C-index with 0.5 credit for tied risk scores."""
    time = np.asarray(time, dtype=float)
    event = np.asarray(event, dtype=int)
    risk = np.asarray(risk, dtype=float)
    n = len(time)
    conc, disc, tied = 0.0, 0.0, 0.0
    for i in range(n):
        if event[i] != 1:
            continue
        for j in range(n):
            if time[i] < time[j]:
                if risk[i] > risk[j]:
                    conc += 1
                elif risk[i] < risk[j]:
                    disc += 1
                else:
                    tied += 1
    total = conc + disc + tied
    return (conc + 0.5 * tied) / total if total > 0 else float("nan")


def concordance_index_fast(time: np.ndarray, event: np.ndarray,
                           risk: np.ndarray) -> float:
    """O(n log n) version via sort + Fenwick-free counting on unique times."""
    time = np.asarray(time, dtype=float)
    event = np.asarray(event, dtype=int)
    risk = np.asarray(risk, dtype=float)
    order = np.argsort(time)
    t_sorted, e_sorted, r_sorted = time[order], event[order], risk[order]
    n = len(t_sorted)
    conc = tied_n = total = 0.0
    for i in range(n):
        if e_sorted[i] != 1:
            continue
        later = t_sorted > t_sorted[i]
        if not later.any():
            continue
        rj = r_sorted[later]
        total += len(rj)
        conc += (r_sorted[i] > rj).sum()
        tied_n += (r_sorted[i] == rj).sum()
    return (conc + 0.5 * tied_n) / total if total > 0 else float("nan")


def breslow_baseline_hazard(time: np.ndarray, event: np.ndarray,
                            eta: np.ndarray) -> tuple:
    """Breslow estimator: cumulative baseline hazard at each event time."""
    order = np.argsort(time)
    t, e, eta_s = time[order], event[order], eta[order]
    exp_eta = np.exp(eta_s)
    risk_sum = np.cumsum(exp_eta[::-1])[::-1]
    ev_times, dH = [], []
    cum = 0.0
    for i in range(len(t)):
        if e[i] == 1:
            cum += 1.0 / max(risk_sum[i], 1e-12)
            ev_times.append(t[i])
            dH.append(cum)
    return np.array(ev_times), np.array(dH)
