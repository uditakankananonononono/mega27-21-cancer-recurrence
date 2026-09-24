import numpy as np
import torch

from recurscan.data.dataset import make_synthetic_survival
from recurscan.models.cnn_surv import CoxCNN
from recurscan.models.deepsurv import DeepSurv
from recurscan.models.gnn_surv import CoxGCN
from recurscan.survival import concordance_index_fast
from recurscan.train import train_cox_fullbatch
from recurscan.graphs import knn_similarity_graph, normalize_adjacency


def test_deepsurv_learns_synthetic():
    X, t, e, _ = make_synthetic_survival(n=500, d=8, seed=11, effect=2.0,
                                         censor_rate=0.2)
    net = DeepSurv(8, hidden=(32,))
    net = train_cox_fullbatch(net, X, t, e, epochs=250, lr=5e-3, seed=0)
    net.eval()
    with torch.no_grad():
        risk = net(torch.tensor(X)).numpy()
    assert concordance_index_fast(t, e, risk) > 0.8


def test_coxcnn_forward_and_grad():
    net = CoxCNN(20, channels=(8,), d_clin=3)
    x = torch.randn(16, 20)
    xc = torch.randn(16, 3)
    out = net(x, xc)
    assert out.shape == (16,)
    out.sum().backward()
    assert all(p.grad is not None for p in net.parameters())


def test_coxgcn_learns_clustered_risk():
    rng = np.random.default_rng(12)
    n = 200
    X = np.vstack([rng.normal(0, 0.3, (n // 2, 6)), rng.normal(3, 0.3, (n // 2, 6))])
    eta = np.array([2.0] * (n // 2) + [0.0] * (n // 2))
    t = rng.weibull(1.5, n) * np.exp(-eta) * 10 + 0.1
    e = np.ones(n, dtype=np.int64)
    e[rng.random(n) < 0.25] = 0  # random censoring
    Xs = (X - X.mean(0)) / X.std(0)
    A = torch.tensor(normalize_adjacency(knn_similarity_graph(Xs, k=8)),
                     dtype=torch.float32)
    net = CoxGCN(6, hidden=16)
    net = train_cox_fullbatch(net, Xs.astype(np.float32), t, e,
                              epochs=300, lr=5e-3, seed=0, graph=A)
    net.eval()
    with torch.no_grad():
        risk = net(torch.tensor(Xs, dtype=torch.float32), A).numpy()
    # the Bayes C-index for this noise level is ~0.72: the model must reach
    # the oracle ceiling, not an arbitrary absolute bar
    oracle = concordance_index_fast(t, e, eta)
    assert concordance_index_fast(t, e, risk) >= oracle - 0.05
