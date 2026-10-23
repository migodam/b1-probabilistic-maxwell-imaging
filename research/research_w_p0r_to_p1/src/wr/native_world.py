"""ORACLE physical-existence chart around a known reference material world.

This is not an online inversion model. It changes the search representation,
not the original noise law, fixed wall support, or spatial material bounds.
No solver is constructed or executed by this module.
"""
from __future__ import annotations

import numpy as np

from .physical_witness import ModulatedWorld, wall_values


class LinearWorld(ModulatedWorld):
    """Affine wall-material edits using 54 nonnegative partition-of-unity hats.

    q_contrast and q_conductivity each have 54 entries in [-1,1]. Tensor
    ordering is y-major / x-minor. Coordinates outside the node rectangle are
    clamped to its edges, so the basis stays nonnegative and sums to one on
    the entire plane. There is no nonlinear coefficient-chart saturation.

    The inherited fixed probe and explicitly declared cavity overrides remain
    identical to ModulatedWorld. Wall contrast <=6.4 and conductivity <=.063
    follow from c0<=4, sigma0<=.036. A fixed probe may have contrast 7 (within
    the original material bound 8); it is not part of the wall modulation.
    """

    def __init__(self, seed, probe=None, cavities=()):
        super().__init__(seed, probe=probe, shape=(6, 9), cavities=cavities)
        self.xnodes = np.linspace(-1.6, 1.6, 9)
        self.ynodes = np.linspace(.4, 1.4, 6)
        self.lower = np.full(self.n, -1.)
        self.upper = np.full(self.n, 1.)

    @staticmethod
    def _hats(coordinates, nodes):
        coordinates = np.clip(coordinates, nodes[0], nodes[-1])
        indices = np.clip(np.searchsorted(nodes, coordinates, side='right') - 1,
                          0, len(nodes) - 2)
        fraction = ((coordinates - nodes[indices])
                    / (nodes[indices + 1] - nodes[indices]))
        hats = np.zeros((len(coordinates), len(nodes)))
        rows = np.arange(len(coordinates))
        hats[rows, indices] = 1. - fraction
        hats[rows, indices + 1] = fraction
        return hats

    def basis(self, pts):
        if np.iscomplexobj(pts):
            raise ValueError('points must be real')
        points = np.asarray(pts, dtype=float)
        if points.ndim != 2 or points.shape[1] != 2 or not np.all(np.isfinite(points)):
            raise ValueError('points must be a finite (N,2) real array')
        hx = self._hats(points[:, 0], self.xnodes)
        hy = self._hats(points[:, 1], self.ynodes)
        return (hy[:, :, None] * hx[:, None, :]).reshape(len(points), self.m)

    def _parameters(self, x):
        if np.iscomplexobj(x):
            raise ValueError('material parameters must be real')
        q = np.asarray(x, dtype=float)
        if q.shape != (self.n,) or not np.all(np.isfinite(q)):
            raise ValueError('expected 108 finite material parameters')
        if np.any(q < -1.) or np.any(q > 1.):
            raise ValueError('linear material chart requires q in [-1,1]')
        return q

    @staticmethod
    def _frequency(f):
        frequency = float(f)
        if not np.isfinite(frequency) or frequency <= 0:
            raise ValueError('frequency must be finite and positive')
        return frequency

    def components(self, x, pts):
        q = self._parameters(x)
        b = self.basis(pts)
        er, sigma0 = wall_values(pts, self.seed)
        # Same explicitly declared cavity overrides as the parent chart.
        for center, radius in self.cavities:
            active = np.linalg.norm(np.asarray(pts) - np.asarray(center), axis=1) < radius
            er[active] = 1.
            sigma0[active] = 0.
        c0 = er - 1.
        a = b @ q[:self.m]
        c = b @ q[self.m:]
        re = 1. + c0 * (1. + .6*a)
        loss = sigma0 * (1. + .75*c)
        return re, loss, b, c0, sigma0, a, c

    def coefficient(self, x, f):
        # Parent applies exactly the same fixed probe material override.
        return super().coefficient(self._parameters(x).copy(), self._frequency(f))

    def derivative(self, x, pts, f):
        frequency = self._frequency(f)
        _, _, b, c0, sigma0, _, _ = self.components(x, pts)
        derivative = np.c_[.6*c0[:, None]*b,
                           (1j*.75/frequency)*sigma0[:, None]*b]
        if self.probe is not None:
            center, radius, _ = self.probe
            active = np.linalg.norm(np.asarray(pts) - np.asarray(center), axis=1) < radius
            derivative[active] = 0.
        return derivative

    def material_metadata(self, x):
        self._parameters(x)
        return dict(model='linear_fixed_support_tensor_hat_partition_of_unity',
                    identity='ORACLE_PHYSICAL_EXISTENCE_KNOWN_REFERENCE_NOT_ONLINE',
                    seed=self.seed, n=self.n, m=self.m,
                    xnodes=self.xnodes.copy(), ynodes=self.ynodes.copy(),
                    tensor_order='y_major_x_minor', parameter_box=[-1., 1.],
                    probe=self.probe, cavities=self.cavities,
                    coordinate_extension='clamp each spatial coordinate to end nodes',
                    derivative='exact affine material derivative; no saturation',
                    same_support_as_parent_world=True,
                    same_support_as_unmodified_reference=not bool(self.cavities) and self.probe is None,
                    modulated_wall_bounds=dict(contrast_max=6.4, conductivity_max=.063,
                                               passive=True),
                    fixed_probe_override='identical inherited override; probe contrast may exceed wall-chart 6.4',
                    original_material_bounds='contrast<=8, conductivity<=.15',
                    original_noise_unchanged=True,
                    prior_probability_status='none assigned; no posterior or prior-probability claim',
                    new_features='none beyond explicitly inherited cavity/probe configuration')
