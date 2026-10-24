"""Explicit observation+anchor conditional densities, separate from raw-y fits.

The anchor is a frozen observation-derived condition, not a hidden GT endpoint.
All66 stochastic coordinates remain modeled. Residual centering is a
translation of determinant1; the learned law still uses raw Gaussianized dx66.
"""
from __future__ import annotations

import numpy as np
import torch

from .density import DensityPredictor

FORMAT = 'p3r4_anchor_full66_density_v1'


def context_features(y, context_x, context_J=None, j_summary='none'):
    y, context_x = np.asarray(y), np.asarray(context_x)
    if y.ndim != 2 or y.shape[1] != 512 or context_x.shape != (len(y), 66):
        raise ValueError('raw y[B,512] and observation-derived context_x[B,66] required')
    if any(np.iscomplexobj(a) or not np.isfinite(a).all() for a in (y, context_x)):
        raise ValueError('finite real raw y and Gaussianized anchor context required')
    fields = [y, context_x]
    if j_summary == 'log_column_norm':
        J = np.asarray(context_J)
        if J.shape != (len(y), 512, 66) or np.iscomplexobj(J) or not np.isfinite(J).all():
            raise ValueError('cached real whitened native Jq[B,512,66] required')
        J = J.astype(np.float64, copy=True)
        # Cached source Jacobians differentiate native q=(z,u), not eta.
        J[:, :, 64:] *= (np.exp(-.5*context_x[:, None, 64:]**2)/np.sqrt(2*np.pi))
        fields.append(np.log1p(np.sum(J*J, axis=1)))
    elif j_summary != 'none':
        raise ValueError('j_summary must be none or log_column_norm')
    return np.concatenate(fields, axis=1)


def base_parameters(context_x, base_mean=None, base_chol=None, anchor_scale=1.):
    mean = np.asarray(context_x if base_mean is None else base_mean)
    if mean.shape != np.asarray(context_x).shape or np.iscomplexobj(mean) or not np.isfinite(mean).all():
        raise ValueError('finite observation-derived base_mean[B,66] required')
    if base_chol is None:
        chol = np.broadcast_to(np.eye(66)*anchor_scale, (len(mean), 66, 66))
    else:
        chol = np.asarray(base_chol)
    if chol.shape != (len(mean), 66, 66) or np.iscomplexobj(chol) or not np.isfinite(chol).all():
        raise ValueError('finite real fullrank base_chol[B,66,66] required')
    if np.any(np.diagonal(chol, axis1=1, axis2=2) <= 0) or not np.allclose(chol, np.tril(chol)):
        raise ValueError('positive-diagonal lower-triangular fullrank Gaussian base required')
    return mean, chol


class AnchorDensityPredictor:
    """Inference consumes only raw y and explicitly y-derived anchor arrays."""
    def __init__(self, payload, device='cpu'):
        if payload.get('format') != FORMAT:
            raise ValueError('explicit anchor-conditioned checkpoint required')
        self.payload = payload
        self.inner = DensityPredictor(payload, device)
        self.model = self.inner.model
        self.device = self.inner.device
        self.contract = payload['anchor_contract']

    def _prepare(self, y, context_x, base_mean, base_chol, context_J):
        features = context_features(y, context_x, context_J, self.contract['j_summary'])
        mean, chol = base_parameters(context_x, base_mean, base_chol, self.contract['anchor_scale'])
        center = mean if self.contract['residual_center'] else np.zeros_like(mean)
        extra = dict(base_mean=mean-center, base_chol=chol) if self.model.kind == 'cfm' else {}
        return features, center, extra

    def log_prob(self, x, y, context_x, *, base_mean=None, base_chol=None,
                 context_J=None, batch_size=16, **kwargs):
        features, center, base = self._prepare(y, context_x, base_mean, base_chol, context_J)
        x = np.asarray(x)
        if x.shape != center.shape:
            raise ValueError('aligned raw Gaussianized x[B,66] required')
        return self.inner.log_prob(x-center, features, batch_size=batch_size, **base, **kwargs)

    def sample(self, y, context_x, *, n_samples=1, seed=0, batch_size=16,
               base_mean=None, base_chol=None, context_J=None, **kwargs):
        features, center, base = self._prepare(y, context_x, base_mean, base_chol, context_J)
        return self.inner.sample(features, n_samples, seed, batch_size, **base, **kwargs)+center[:, None]


def load_anchor_density(path, device='cpu'):
    return AnchorDensityPredictor(torch.load(path, map_location='cpu', weights_only=False), device)
