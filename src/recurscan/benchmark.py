"""Benchmark: Cox baselines vs deep Cox models on METABRIC RFS."""
from __future__ import annotations

import json
import time

import numpy as np
import torch
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index as ll_ci
import pandas as pd

from .data.dataset import SurvivalDataset
from .eval import bootstrap_cindex_ci, summarize, survival_metrics
from .graphs import knn_similarity_graph, normalize_adjacency
from .models.cnn_surv import CoxCNN
from .models.deepsurv import DeepSurv
from .models.gnn_surv import CoxGCN
from .train import set_seed, train_cox_fullbatch

MODELS = ("coxph_clin", "coxph_full", "deepsurv", "coxcnn", "coxgcn")


def stratified_event_split(event: np.ndarray, test_frac: float, seed: int):
    rng = np.random.default_rng(seed)
    tr, te = [], []
    for cls in (0, 1):
        idx = np.where(event == cls)[0]
        rng.shuffle(idx)
        n_te = max(1, int(round(test_frac * len(idx))))
        te.extend(idx[:n_te].tolist())
        tr.extend(idx[n_te:].tolist())
    return np.array(sorted(tr)), np.array(sorted(te))


def standardize(Xtr, Xte):
    mu = Xtr.mean(0, keepdims=True)
    sd = np.where(Xtr.std(0, keepdims=True) < 1e-8, 1.0, Xtr.std(0, keepdims=True))
    return (Xtr - mu) / sd, (Xte - mu) / sd


def run_one_seed(ds: SurvivalDataset, model: str, seed: int,
                 test_frac: float = 0.25, k: int = 8) -> dict:
    set_seed(seed)
    tr, te = stratified_event_split(ds.event, test_frac, seed)
    t_tr, t_te = ds.time[tr], ds.time[te]
    e_tr, e_te = ds.event[tr], ds.event[te]

    if model in ("coxph_clin", "coxph_full"):
        feats = ["clin"] if model == "coxph_clin" else ["clin", "expr"]
        blocks = [ds.X_clin] + ([ds.X_expr] if "expr" in feats else [])
        X = np.concatenate(blocks, axis=1)
        Xtr, Xte = standardize(X[tr], X[te])
        names = [f"x{i}" for i in range(X.shape[1])]
        df = pd.DataFrame(Xtr, columns=names)
        df["time"], df["event"] = t_tr, e_tr
        cph = CoxPHFitter(penalizer=0.1)
        cph.fit(df, "time", "event")
        risk_te = cph.predict_partial_hazard(pd.DataFrame(Xte, columns=names)).values.ravel()
        risk = np.log(risk_te)
    else:
        Xtr_e, Xte_e = standardize(ds.X_expr[tr], ds.X_expr[te])
        Ctr, Cte = standardize(ds.X_clin[tr], ds.X_clin[te])
        d_expr, d_clin = ds.X_expr.shape[1], ds.X_clin.shape[1]
        if model == "deepsurv":
            net = DeepSurv(d_expr + d_clin)
            net = train_cox_fullbatch(net, np.hstack([Xtr_e, Ctr]), t_tr, e_tr,
                                      epochs=300, lr=1e-3, seed=seed)
            net.eval()
            with torch.no_grad():
                risk = net(torch.tensor(np.hstack([Xte_e, Cte]),
                                        dtype=torch.float32)).numpy()
        elif model == "coxcnn":
            net = CoxCNN(d_expr, d_clin=d_clin)
            net = train_cox_fullbatch(net, Xtr_e, t_tr, e_tr,
                                      X_clin_train=Ctr, epochs=300, lr=1e-3,
                                      seed=seed)
            net.eval()
            with torch.no_grad():
                risk = net(torch.tensor(Xte_e, dtype=torch.float32),
                           torch.tensor(Cte, dtype=torch.float32)).numpy()
        elif model == "coxgcn":
            # transductive over train+test nodes, Cox loss on train only
            Xall = np.vstack([Xtr_e, Xte_e])
            Call = np.vstack([Ctr, Cte])
            feats = np.hstack([Xall, Call])
            A = torch.tensor(normalize_adjacency(
                knn_similarity_graph(feats, k=k)), dtype=torch.float32)
            ntr = len(Xtr_e)
            tall = np.concatenate([t_tr, t_te])
            eall = np.concatenate([e_tr, e_te])
            net = CoxGCN(feats.shape[1])
            net = _train_gcn_cox(net, feats, tall, eall, ntr, A, seed)
            net.eval()
            with torch.no_grad():
                risk = net(torch.tensor(feats, dtype=torch.float32), A)[ntr:].numpy()
        else:
            raise KeyError(model)

    m = survival_metrics(t_te, e_te, risk)
    m["c_index_ci95"] = bootstrap_cindex_ci(t_te, e_te, risk, seed=seed)
    return m


def _train_gcn_cox(net, feats, time, event, n_train, A, seed,
                   epochs: int = 400, lr: float = 5e-3):
    from .survival import cox_partial_loglik
    set_seed(seed)
    Xt = torch.tensor(feats, dtype=torch.float32)
    tt = torch.tensor(time, dtype=torch.float32)
    et = torch.tensor(event, dtype=torch.float32)
    opt = torch.optim.Adam(net.parameters(), lr=lr, weight_decay=1e-4)
    best, stall, state = np.inf, 0, None
    for _ in range(epochs):
        net.train()
        eta = net(Xt, A)[:n_train]
        loss = -cox_partial_loglik(eta, tt[:n_train], et[:n_train])
        opt.zero_grad()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(net.parameters(), 5.0)
        opt.step()
        if loss.item() < best - 1e-5:
            best, stall = loss.item(), 0
            state = {k_: v.detach().clone() for k_, v in net.state_dict().items()}
        else:
            stall += 1
            if stall >= 60:
                break
    if state is not None:
        net.load_state_dict(state)
    return net


def run_benchmark(ds: SurvivalDataset, models=MODELS, seeds=(0, 1, 2, 3, 4)):
    t0 = time.time()
    out = {"dataset": ds.name, "n": len(ds.time),
           "events": int(ds.event.sum()), "models": {}, "seeds": list(seeds)}
    for m in models:
        runs = [run_one_seed(ds, m, s) for s in seeds]
        out["models"][m] = {"per_seed": runs, "summary": summarize(runs)}
        s = out["models"][m]["summary"]["c_index"]
        print(f"{m:12s} C-index {s['mean']:.3f}+-{s['std']:.3f}", flush=True)
    out["wall_seconds"] = round(time.time() - t0, 2)
    return out


def save(result: dict, path: str):
    with open(path, "w") as fh:
        json.dump(result, fh, indent=2)
