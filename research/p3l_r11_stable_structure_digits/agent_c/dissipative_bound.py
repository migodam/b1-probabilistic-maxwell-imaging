"""Float64 discrete dissipative residual bounds; no PDE/interval certificate.

The total-field residual is E - (I - delta C diag(a)) e_approx.
No inverse of sqrt(a) is used, including at support birth.
"""
from __future__ import annotations

import time
import numpy as np
from scipy.linalg import eigh


def spectral_norm_rows(X):
    """2-norm through the smaller row Gram; X has only 32 receiver rows."""
    G = X @ X.conj().T
    ev = float(eigh(G, eigvals_only=True,
                    subset_by_index=(G.shape[0] - 1, G.shape[0] - 1),
                    check_finite=False)[0])
    return float(np.sqrt(max(0.0, ev)))


class GlobalDissipativeBound:
    """Uniform-in-a,delta stability with either current or uniform port norm.

    mode='global': beta_star=min(0,lambda_min(Im_H C)), b_star is uniform
    over 0<=a<=1 and the declared delta rectangle.
    mode='local': use lambda_min(Im_H(P C P)) and the current delta.
    The local comparator does not alter the global definition.
    """

    def __init__(self, C, O, sigma, *, mode='global', uniform_port=False,
                 rectangle=((.8, 1.8), (.03, .18)), white_gate=.01,
                 potential_gate=.1):
        t0 = time.perf_counter()
        self.C = np.asarray(C, dtype=complex)
        self.O = np.asarray(O, dtype=complex)
        self.sigma = float(sigma)
        self.mode = mode
        self.uniform_port = bool(uniform_port)
        self.rectangle = np.asarray(rectangle, dtype=float)
        self.white_gate = float(white_gate)
        self.potential_gate = float(potential_gate)
        if mode not in ('global', 'local'):
            raise ValueError('mode must be global or local')
        if self.C.ndim != 2 or self.C.shape[0] != self.C.shape[1]:
            raise ValueError('C must be square')
        if self.O.shape[1] != self.C.shape[0] or not self.sigma > 0:
            raise ValueError('O shape or noise sigma invalid')
        if not np.isfinite(self.C).all() or not np.isfinite(self.O).all():
            raise ValueError('nonfinite fixed operators')
        (xmin, xmax), (ymin, ymax) = self.rectangle
        if not (0 < ymin <= ymax < xmin <= xmax):
            raise ValueError('rectangle requires 0<ymin<=ymax<xmin<=xmax')
        self.amplitude_loss_min = float(ymin / (xmax*xmax + ymin*ymin))
        self.lambda_min_ImC = float(eigh((self.C-self.C.conj().T)/(2j),
                    eigvals_only=True, subset_by_index=(0, 0),
                    check_finite=False)[0])
        self.beta_star = min(0.0, self.lambda_min_ImC)
        self.b_star = self.amplitude_loss_min + self.beta_star
        self.O_norm = spectral_norm_rows(self.O)
        self.setup_seconds = time.perf_counter() - t0
        self.calls = 0
        self.gate_seconds = 0.0
        self.shape_setup_seconds = 0.0
        self.shape_setups = 0
        self._a = None
        self._port_norm = None
        self._local_beta = None

    def _shape(self, a):
        if self._a is not None and np.array_equal(self._a, a):
            return
        t0 = time.perf_counter()
        p = np.sqrt(a)
        self._port_norm = self.O_norm if self.uniform_port else spectral_norm_rows(self.O*p)
        if self.mode == 'local':
            B = p[:, None] * self.C * p[None, :]
            self._local_beta = float(eigh((B-B.conj().T)/(2j),
                        eigvals_only=True, subset_by_index=(0, 0),
                        check_finite=False)[0])
        self._a = a.copy()
        self.shape_setups += 1
        self.shape_setup_seconds += time.perf_counter() - t0

    def __call__(self, a, delta, residual, approxF, y=None):
        t0 = time.perf_counter()
        self.calls += 1
        a = np.asarray(a, dtype=float).reshape(-1)
        delta = complex(delta)
        residual = np.asarray(residual)
        approxF = np.asarray(approxF)
        (xmin, xmax), (ymin, ymax) = self.rectangle
        valid = bool(a.shape == (self.C.shape[0],) and np.isfinite(a).all()
                     and np.all(a >= 0) and np.all(a <= 1)
                     and np.isfinite(delta) and xmin <= delta.real <= xmax
                     and ymin <= delta.imag <= ymax
                     and residual.shape[0] == self.C.shape[0]
                     and np.isfinite(residual).all() and np.isfinite(approxF).all())
        result = dict(accepted=False, qualified=False, gate=False,
                      mode=self.mode, valid_domain=valid,
                      beta_star=self.beta_star, b_star=self.b_star,
                      lambda_min_ImC=self.lambda_min_ImC,
                      finite_float_conditional=True,
                      interval_or_PDE_certificate=False,
                      white_bound=None, potential_bound=None,
                      numerical_coercivity=None, failure=None)
        if not valid:
            result['failure'] = 'OUTSIDE_DECLARED_DOMAIN_OR_NONFINITE'
        else:
            self._shape(a)
            b = self.b_star if self.mode == 'global' else delta.imag/abs(delta)**2+self._local_beta
            rn = float(np.linalg.norm(np.sqrt(a)[:, None]*residual, 'fro'))
            result.update(numerical_coercivity=float(b),
                          residual_weighted_norm=rn,
                          residual_norm=float(np.linalg.norm(residual, 'fro')),
                          port_norm=self._port_norm,
                          beta_local=self._local_beta)
            if not np.isfinite(b) or b <= 0:
                result['failure'] = 'NONPOSITIVE_DISSIPATIVE_BOUND'
            else:
                ew = float(np.sqrt(2)*self._port_norm*rn/(self.sigma*b))
                result['white_bound'] = ew
                if y is None:
                    result['failure'] = 'OBSERVATION_REQUIRED_FOR_POTENTIAL_GATE'
                else:
                    y = np.asarray(y)
                    if y.shape != approxF.shape or not np.isfinite(y).all():
                        result['failure'] = 'INVALID_OBSERVATION'
                    else:
                        rw = float(np.sqrt(2)/self.sigma*np.linalg.norm(approxF-y, 'fro'))
                        pot = float(rw*ew + .5*ew*ew)
                        result.update(potential_bound=pot, approx_white_data_residual=rw)
                        ok = bool(np.isfinite(ew) and np.isfinite(pot)
                                  and ew <= self.white_gate and pot <= self.potential_gate)
                        result.update(accepted=ok, qualified=ok, gate=ok,
                                      failure=None if ok else 'TASK_BOUND_GATE')
        dt = time.perf_counter()-t0
        self.gate_seconds += dt
        result['gate_seconds'] = dt
        return result

    def summary(self):
        return dict(mode=self.mode, uniform_port=self.uniform_port,
                    rectangle=self.rectangle.tolist(), sigma=self.sigma,
                    lambda_min_ImC=self.lambda_min_ImC,
                    beta_star=self.beta_star,
                    amplitude_loss_min=self.amplitude_loss_min,
                    b_star=self.b_star, global_bound_positive=self.b_star > 0,
                    O_norm=self.O_norm, setup_seconds=self.setup_seconds,
                    shape_setups=self.shape_setups,
                    shape_setup_seconds=self.shape_setup_seconds,
                    calls=self.calls, gate_seconds=self.gate_seconds,
                    white_gate=self.white_gate, potential_gate=self.potential_gate,
                    qualification='NUMERICAL_DISCRETE_FLOAT64_CONDITIONAL_BOUND',
                    interval_or_PDE_certificate=False)


class LocalDissipativeBound(GlobalDissipativeBound):
    def __init__(self, C, O, sigma, **kwargs):
        kwargs.setdefault('mode', 'local')
        super().__init__(C, O, sigma, **kwargs)
