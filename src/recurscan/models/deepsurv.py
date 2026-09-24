"""DeepSurv: MLP risk scorer trained with the Cox partial likelihood."""
from __future__ import annotations

import torch
import torch.nn as nn


class DeepSurv(nn.Module):
    def __init__(self, d_in: int, hidden: tuple = (64, 32), dropout: float = 0.2):
        super().__init__()
        layers, prev = [], d_in
        for h in hidden:
            layers += [nn.Linear(prev, h), nn.BatchNorm1d(h), nn.ELU(),
                       nn.Dropout(dropout)]
            prev = h
        layers.append(nn.Linear(prev, 1))
        self.net = nn.Sequential(*layers)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x).squeeze(-1)  # risk eta (n,)
