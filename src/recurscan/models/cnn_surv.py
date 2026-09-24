"""1D-CNN Cox model over expression ordered along the panel's gene axis."""
from __future__ import annotations

import torch
import torch.nn as nn


class CoxCNN(nn.Module):
    def __init__(self, d_in: int, channels: tuple = (16, 32), kernel: int = 5,
                 dropout: float = 0.2, d_clin: int = 0):
        super().__init__()
        blocks, prev = [], 1
        for c in channels:
            blocks += [nn.Conv1d(prev, c, kernel, padding=kernel // 2),
                       nn.BatchNorm1d(c), nn.ELU()]
            prev = c
        self.features = nn.Sequential(*blocks)
        self.head = nn.Sequential(
            nn.Linear(prev + d_clin, 32), nn.ELU(), nn.Dropout(dropout),
            nn.Linear(32, 1))
        self.d_clin = d_clin

    def forward(self, x: torch.Tensor, x_clin: torch.Tensor = None):
        h = self.features(x.unsqueeze(1)).mean(dim=-1)
        if self.d_clin and x_clin is not None:
            h = torch.cat([h, x_clin], dim=1)
        return self.head(h).squeeze(-1)
