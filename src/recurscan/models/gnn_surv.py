"""GCN Cox model over a patient-similarity graph (transductive)."""
from __future__ import annotations

import torch
import torch.nn as nn


class GCNLayer(nn.Module):
    def __init__(self, d_in: int, d_out: int):
        super().__init__()
        self.weight = nn.Parameter(torch.empty(d_in, d_out))
        self.bias = nn.Parameter(torch.zeros(d_out))
        nn.init.xavier_uniform_(self.weight)

    def forward(self, x, a_norm):
        return a_norm @ (x @ self.weight) + self.bias


class CoxGCN(nn.Module):
    def __init__(self, d_in: int, hidden: int = 32, dropout: float = 0.3):
        super().__init__()
        self.gc1 = GCNLayer(d_in, hidden)
        self.gc2 = GCNLayer(hidden, 1)
        self.dropout = nn.Dropout(dropout)

    def forward(self, x, a_norm):
        h = torch.relu(self.gc1(x, a_norm))
        h = self.dropout(h)
        return self.gc2(h, a_norm).squeeze(-1)
