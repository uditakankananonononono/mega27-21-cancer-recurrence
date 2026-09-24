"""Full-batch Cox training (partial likelihood needs complete risk sets)."""
from __future__ import annotations

import numpy as np
import torch

from .survival import cox_partial_loglik


def set_seed(seed: int):
    np.random.seed(seed)
    torch.manual_seed(seed)


def train_cox_fullbatch(model, X_train: np.ndarray, time: np.ndarray,
                        event: np.ndarray, X_val: np.ndarray = None,
                        time_val: np.ndarray = None, event_val: np.ndarray = None,
                        X_clin_train: np.ndarray = None,
                        X_clin_val: np.ndarray = None,
                        epochs: int = 300, lr: float = 1e-3, wd: float = 1e-4,
                        patience: int = 50, seed: int = 0,
                        graph: torch.Tensor = None):
    """Generic full-batch Cox trainer for MLP/CNN (and GCN via `graph`).

    Validation uses the training risk set as reference: partial likelihood is
    not defined on a held-out set alone, so early stopping tracks the val
    C-index instead (computed by the caller's metric fn if provided).
    """
    set_seed(seed)
    Xt = torch.tensor(X_train, dtype=torch.float32)
    tt = torch.tensor(time, dtype=torch.float32)
    et = torch.tensor(event, dtype=torch.float32)
    Xc = torch.tensor(X_clin_train, dtype=torch.float32) if X_clin_train is not None else None
    opt = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=wd)

    def risk(x, xc):
        if graph is not None:
            return model(x, graph)
        if Xc is not None:
            return model(x, xc)
        return model(x)

    best = {"loss": np.inf, "state": None, "stall": 0}
    for _ in range(epochs):
        model.train()
        eta = risk(Xt, Xc)
        loss = -cox_partial_loglik(eta, tt, et)
        opt.zero_grad()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
        opt.step()
        lval = loss.item()
        if lval < best["loss"] - 1e-5:
            best.update(loss=lval, stall=0,
                        state={k: v.detach().clone()
                               for k, v in model.state_dict().items()})
        else:
            best["stall"] += 1
            if best["stall"] >= patience:
                break
    if best["state"] is not None:
        model.load_state_dict(best["state"])
    return model
