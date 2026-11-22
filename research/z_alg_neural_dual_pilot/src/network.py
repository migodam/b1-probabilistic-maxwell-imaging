"""Physics-aligned proposals for coarse piecewise constant dual multipliers."""
from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import numpy as np
import torch
from torch import nn
from torch.nn import functional as F

from chart import LAMBDA_BIAS, LAMBDA_SLOPE, MU_FLOOR, MU_OFFSET, MU_SCALE

PAPER_PUBLIC = (Path(__file__).resolve().parents[2]
                / 'z_alg_optical_inertia/synthetic/public')
FEATURE_SIDE = 8
COARSE_SIDE = 32
WIDTH = 12
FEATURE_CHANNELS = 12
NETWORK_VERSION = 1


def block_side(side):
    if side in (16, 64):
        side = int(np.sqrt(side))
    if side not in (4, 8):
        raise ValueError('Expected a 4×4 or 8×8 block grid.')
    return side


@lru_cache(maxsize=1)
def feature_geometry():
    with np.load(PAPER_PUBLIC / 'pw360_10x20.npz', allow_pickle=False) as data:
        s, b = data['S'], data['B']
    t = np.linspace(-1 + 1/COARSE_SIDE, 1 - 1/COARSE_SIDE, COARSE_SIDE)
    row, col = np.meshgrid(t, t, indexing='ij')
    receiver_norm = np.sqrt(np.mean(np.abs(s)**2, axis=0))
    incident_norm = np.sqrt(np.mean(np.abs(b)**2, axis=1))
    known = np.stack((col, row, receiver_norm.reshape(row.shape)/receiver_norm.mean() - 1,
                      incident_norm.reshape(row.shape)/incident_norm.mean() - 1))
    return s.conj().T, b.conj(), incident_norm, known


def observation_features(Y, u, noise_radius=None, *, use_observation=True):
    """Pool source-normalized fields in row-major blocks of the public grid."""
    y = np.asarray(Y)
    if y.shape != (20, 10):
        raise ValueError(f'Expected complex 20×10 port observations, received {y.shape}.')
    adjoint, incident_conj, incident_norm, known = feature_geometry()
    shape = (COARSE_SIDE, COARSE_SIDE)
    if use_observation:
        ynorm = float(np.linalg.norm(y))
        q = adjoint @ y
        local_amplitude = np.sqrt(np.mean(np.abs(q)**2, axis=1))
        coarse_scale = max(float(np.sqrt(np.mean(local_amplitude**2))), 1e-100)
        c = np.mean(incident_conj * q, axis=1) / (coarse_scale*incident_norm)
        phase_lock = np.abs(c) / np.maximum(local_amplitude/coarse_scale, 1e-8)
        source_axis_variance = np.std(np.abs(q), axis=1)/coarse_scale
        response = np.stack((c.real, c.imag, np.log1p(local_amplitude/coarse_scale),
                             phase_lock, np.log1p(source_axis_variance)))
        scalars = (np.log(max(ynorm, 1e-100))/4,
                   float(noise_radius)/max(ynorm, 1e-100)*10 if noise_radius is not None else 0.)
    else:
        response = np.zeros((5, COARSE_SIDE**2))
        scalars = (0., 0.)
    maps = np.stack([np.full(shape, value) for value in ((float(u)-0.75)/0.25, *scalars)])
    full = np.concatenate((response.reshape(5, *shape), known, maps))
    cells = COARSE_SIDE // FEATURE_SIDE
    pooled = full.reshape(FEATURE_CHANNELS, FEATURE_SIDE, cells, FEATURE_SIDE, cells).mean(axis=(2, 4))
    return np.ascontiguousarray(pooled, dtype=np.float32)


class DualProposal(nn.Module):
    """A shared spatial encoder with individual 4×4 and 8×8 readout heads."""

    def __init__(self, *, use_observation=True):
        super().__init__()
        self.use_observation = bool(use_observation)
        self.encoder = nn.Sequential(
            nn.Conv2d(FEATURE_CHANNELS, WIDTH, 3, padding=1),
            nn.SiLU(),
            nn.Conv2d(WIDTH, WIDTH, 3, padding=1),
            nn.SiLU(),
        )
        self.heads = nn.ModuleDict({str(side): nn.Conv2d(WIDTH, 2, 1) for side in (4, 8)})
        for head in self.heads.values():
            nn.init.zeros_(head.weight)
            nn.init.zeros_(head.bias)

    def proposals(self, features, side):
        return self.from_encoded(self.encoder(features), side)

    def from_encoded(self, encoded, side):
        side = block_side(side)
        z = F.adaptive_avg_pool2d(encoded, side)
        raw = self.heads[str(side)](z).flatten(start_dim=2)
        if raw.dtype != torch.float64:
            raw = raw.to(torch.float64)
        a, b = raw.unbind(dim=1)
        mu = MU_SCALE * (MU_FLOOR + F.softplus(MU_OFFSET + a))
        lam = MU_SCALE * (LAMBDA_BIAS + LAMBDA_SLOPE * b)
        return mu, lam, raw

    def checkpoint(self):
        return dict(network_version=NETWORK_VERSION, use_observation=self.use_observation,
                    chart=dict(ref_mu=MU_SCALE, mu_floor_factor=MU_FLOOR,
                               lambda_ref_ratio=LAMBDA_BIAS, lambda_coord_scale=LAMBDA_SLOPE),
                    state_dict=self.state_dict())


def load_network(path):
    with torch.inference_mode():
        state = torch.load(path, map_location='cpu', weights_only=True)
        if state['network_version'] != NETWORK_VERSION:
            raise ValueError('Checkpoint construction differs from the current proposal network.')
        model = DualProposal(use_observation=state['use_observation'])
        if model.checkpoint()['chart'] != state['chart']:
            raise ValueError('Checkpoint multiplier chart differs from the current implementation.')
        model.load_state_dict(state['state_dict'])
    return model.eval()


def predict(model, Y, u, side, noise_radius=None):
    """Return detached coarse physical µ and λ arrays in row-major block order."""
    features = observation_features(Y, u, noise_radius, use_observation=model.use_observation)
    with torch.inference_mode():
        mu, lam, _ = model.proposals(torch.from_numpy(features[None]), side)
    return mu[0].numpy().copy(), lam[0].numpy().copy()
