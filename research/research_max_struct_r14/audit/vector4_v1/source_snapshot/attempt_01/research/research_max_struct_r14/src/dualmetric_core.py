"""Stateful R14 metric factorial over the inherited finite FFT/PR loop.

The current step is the R8 native Fourier/DPCS, real PR+, exact quadratic
line search.  Material steps call the original scalar readout or R7's full
Möbius/Arc readout.  This module neither reads truth nor audits Maxwell states;
the caller owns paid physical audits and candidate selection.
"""
from collections import Counter
from copy import deepcopy
from pathlib import Path
import json
import re
import time

import numpy as np
from scipy.linalg import svd

from research.research_max_struct_r2.src.cie_global import (
    _source_metric, _inner, _regularizer, _baseline_readout)
from research.research_max_struct_r4.src.cie_onehop import (
    OneHopAbsorption, onehop_material_value_gradient)
from research.research_max_struct_r6.src.fft_tsom_ops import (
    NativeFourierAPCS, DiskPulseToeplitzGreen)
from research.research_max_struct_r7.src.aplus_readout import mobius_material_readout


PAIRS = ('00', '01', '10', '11')
_RUNTIME_KEYS = {
    '_counter', '_aplus_prior_context', '_initial_physical_j',
    'budget_check', 'iteration_callback', 'stage_callback', 'cost_callback',
    'physical_audit_due', 'physical_whiten', 'regularizer', 'aplus_regularizer',
    'material_readout', 'receipt_path', 'max_wall_seconds', 'max_iterations',
    'physical_sigma',
}


def _save_state_npz(path, state):
    """Safe recursive tree plus named plain arrays; no object-array payloads."""
    arrays = {}
    def encode(value, name='state'):
        if isinstance(value, np.ndarray):
            if value.dtype.hasobject:
                raise TypeError('Checkpoint arrays must not contain Python objects')
            base = re.sub(r'[^A-Za-z0-9_]', '_', name)
            key = base
            serial = 1
            while key in arrays or key == '__state_tree_json__':
                key = base+'__'+str(serial)
                serial += 1
            arrays[key] = value
            return {'type': 'ndarray', 'key': key}
        if isinstance(value, np.generic):
            return encode(value.item(), name)
        if value is None or isinstance(value, (str, bool, int)):
            return value
        if isinstance(value, float):
            return value if np.isfinite(value) else {'type': 'float', 'value': repr(value)}
        if isinstance(value, complex):
            return {'type': 'complex', 'real': encode(value.real), 'imag': encode(value.imag)}
        if isinstance(value, (list, tuple)):
            return {'type': 'tuple' if isinstance(value, tuple) else 'list',
                'items': [encode(v, name+'__'+str(i)) for i, v in enumerate(value)]}
        if isinstance(value, dict):
            return {'type': 'dict', 'items': [[encode(k, name+'__key_'+str(i)),
                encode(v, str(k) if name == 'state' and isinstance(k, str)
                    else name+'__'+(str(k) if isinstance(k, str) else 'item_'+str(i)))]
                for i, (k, v) in enumerate(value.items())]}
        raise TypeError('Unsupported checkpoint value type: '+type(value).__name__)
    tree = encode(state)
    arrays['__state_tree_json__'] = np.asarray(json.dumps(tree, allow_nan=False, separators=(',', ':')))
    with Path(path).open('wb') as stream:
        np.savez_compressed(stream, **arrays)


def _load_state_npz(path):
    """Decode only the fixed JSON tags and numeric/string NPZ arrays."""
    with np.load(path, allow_pickle=False) as archive:
        tree = json.loads(str(archive['__state_tree_json__'].item()))
        def decode(value):
            if not isinstance(value, dict):
                return value
            tag = value['type']
            if tag == 'ndarray':
                return archive[value['key']].copy()
            if tag == 'float':
                if value['value'] not in ('inf', '-inf', 'nan'):
                    raise ValueError('Unsupported nonfinite checkpoint float')
                return float(value['value'])
            if tag == 'complex':
                return complex(decode(value['real']), decode(value['imag']))
            if tag == 'list':
                return [decode(v) for v in value['items']]
            if tag == 'tuple':
                return tuple(decode(v) for v in value['items'])
            if tag == 'dict':
                return {decode(k): decode(v) for k, v in value['items']}
            raise ValueError('Unsupported checkpoint JSON tag: '+str(tag))
        return decode(tree)


def _plain_config(config):
    """Freeze numerical flags without pickling callbacks or prior objects."""
    def plain(value):
        if value is None or isinstance(value, (str, bool, int, float, complex, np.generic)):
            return deepcopy(value)
        if isinstance(value, np.ndarray):
            return value.copy()
        if isinstance(value, (list, tuple)):
            converted = [plain(v) for v in value]
            return tuple(converted) if isinstance(value, tuple) else converted
        if isinstance(value, dict):
            return {k: plain(v) for k, v in value.items()}
        raise TypeError
    result = {}
    for key, value in config.items():
        if key in _RUNTIME_KEYS or callable(value):
            continue
        try:
            result[key] = plain(value)
        except TypeError:
            # Runtime objects are supplied again by the caller on restore.
            pass
    return result


class StatefulDualMetric:
    """One completed current/material update per :meth:`step`.

    Cold DPCS initialization is the default, as in historical A+-Arc.  An
    explicitly requested physical-current warm start reuses
    ``config['_initial_physical_j']`` or the matching paid ``model._current``;
    only when neither exists does it call the original ``model.state(m0)``.

    Stage-end bookkeeping happens in the step that completes that stage.
    The next step transports the Fourier coefficients, resets PR memory, and
    enters the next stage.  A saved stage-end state therefore has unambiguous
    next-iteration semantics.  Exhaustion raises ``StopIteration``.  Callbacks
    are optional; returned rows contain copied state arrays, while retained
    history contains only the native scalar diagnostics.
    """

    checkpoint_version = 1

    def __init__(self, model, Y, m0, config, pair='01'):
        self._bind(model, Y, config)
        self._validate_pair(pair)
        self.pair = pair
        self.cost = Counter()
        self.config['_counter'] = self.cost
        self.m = np.asarray(m0, float).copy()
        if self.m.shape != (self.n,) or not np.all(np.isfinite(self.m)):
            raise ValueError('Invalid full real initial material')
        if np.any(self.m < self.lo) or np.any(self.m > self.hi):
            raise ValueError('Initial material must belong to the common box')
        self._check()
        geometry = self.config['fft_geometry']
        self.green = DiskPulseToeplitzGreen(geometry['coords'], self.shape,
            float(geometry['radius']), float(geometry['k']), self.cost)
        self.blocks = OneHopAbsorption(self.G, self.shape, self.cost)
        if self.config.get('onehop_self_only', False):
            self.blocks.goff[:] = 0
        self.Od, self.Ld = _source_metric(self.config, 'data_metric', 'data_weights',
            self.ns, 1/np.maximum(np.sum(abs(self.Y)**2, axis=0), 1e-300), self.cost)
        self.Os, self.Ls = _source_metric(self.config, 'state_metric', 'state_weights',
            self.ns, 1/np.maximum(np.sum(abs(self.E)**2, axis=0), 1e-300), self.cost)
        self._bind_weights()
        self.initial_physical_j = None
        if self.config.get('initialize_physical_current', False):
            paid = self.config.get('_initial_physical_j')
            if paid is None:
                state = getattr(model, '_current', None)
                if state is None or not np.array_equal(state.get('m'), self.m):
                    state = model.state(self.m)
                    self.cost['physical_validations'] += 1
                paid = state['j']
            self.initial_physical_j = np.asarray(paid, complex).copy()
            if self.initial_physical_j.shape != self.E.shape:
                raise ValueError('Initial physical current shape mismatch')
        Us, sig, Vh = svd(self.S, full_matrices=False, check_finite=False)
        self.cost['receiver_SVD'] += 1
        if 'dpcs_rank' in self.config:
            self.receiver_rank = int(self.config['dpcs_rank'])
        elif 'dpcs_rtol' in self.config:
            self.receiver_rank = int(np.count_nonzero(sig > float(self.config['dpcs_rtol'])*sig[0]))
        else:
            raise ValueError('Freeze dpcs_rank or dpcs_rtol from public noise/calibration')
        if not 0 <= self.receiver_rank <= len(sig):
            raise ValueError('Invalid DPCS rank')
        if self.receiver_rank and sig[self.receiver_rank-1] <= np.finfo(float).eps*max(self.S.shape)*sig[0]:
            raise ValueError('DPCS attempts to invert a numerically null receiver mode')
        self.receiver_left = Us[:, :self.receiver_rank].copy()
        self.receiver_singular_values = sig.copy()
        self.Vd = Vh[:self.receiver_rank].conj().T
        self.JD = self.Vd @ ((self.receiver_left.conj().T @ self.Y)/sig[:self.receiver_rank, None])
        if self.config.get('current_anchor') is not None:
            raise ValueError('R14 main factorial does not change the historical affine current anchor')
        self.GJD, self.SJD = self.green.apply(self.JD), self.S@self.JD
        self.cost['G_columns'] += self.ns
        self.cost['S_columns'] += self.ns
        self.old_coefficients = {}
        self.old_grad = self.old_direction = None
        self.iterations = 0
        self.stage = -1
        self.iteration = 0
        self.history, self.stage_states = [], []
        self.stage_pending = False
        self.done = False
        self.status = 'SCHEDULE_COMPLETE'
        self.stage_stop = 'MAXITER'
        self._enter_stage(0)
        self._check()

    def _bind(self, model, Y, config):
        self.model = model
        self.config = dict(config)
        self.config.pop('_aplus_prior_context', None)
        self.G, self.E, self.S = [np.asarray(a, complex) for a in (model.G, model.E, model.S)]
        self.Y = np.asarray(Y, complex)
        self.n, self.ns = self.E.shape
        self.shape = tuple(map(int, model.shape))
        if self.G.shape != (self.n, self.n) or self.S.shape[1] != self.n or self.Y.shape != (self.S.shape[0], self.ns):
            raise ValueError('Inconsistent public operators/observations')
        if len(self.shape) != 2 or np.prod(self.shape) != self.n or any(s % 2 for s in self.shape):
            raise ValueError('R14 requires the original even 2D grid')
        if self.config.get('physics_family', 'original_disk_pulse') != 'original_disk_pulse':
            raise ValueError('R14 core reuses only the original finite disk-pulse physics')
        if self.config.get('method', 'cie').lower() != 'cie' or self.config.get('current_solver', 'pr') != 'pr':
            raise ValueError('R14 preserves the original CIE PR solver')
        self.stages = [dict(s) for s in self.config['stages']]
        if not self.stages:
            raise ValueError('At least one preregistered stage required')
        self.lo, self.hi = [np.broadcast_to(np.asarray(x, float), (self.n,)).copy()
            for x in self.config.get('material_box', (.02, 20.))]
        if np.any(self.lo > self.hi):
            raise ValueError('Invalid real material box')
        self.reg = _regularizer(self.config)
        self.has_reg = self.config.get('regularizer') is not None
        self.rho = float(self.config.get('rho', 1.))
        if self.rho <= 0:
            raise ValueError('Positive state penalty required')
        self.sourcewise = bool(self.config.get('sourcewise_cg', False))
        self.started = time.perf_counter()
        self.elapsed_before = 0.

    def _bind_weights(self):
        if self.sourcewise and (not np.allclose(self.Od, np.diag(np.diag(self.Od)), rtol=0, atol=0)
                or not np.allclose(self.Os, np.diag(np.diag(self.Os)), rtol=0, atol=0)):
            raise ValueError('Historical sourcewise PR requires diagonal source metrics')
        self.wd, self.ws = np.diag(self.Od).real, np.diag(self.Os).real

    @staticmethod
    def _validate_pair(pair):
        if pair not in PAIRS:
            raise ValueError('Metric pair must be one of '+repr(PAIRS))

    @property
    def global_iteration(self):
        return self.iterations

    @property
    def counts(self):
        return self.cost

    @property
    def coefficients(self):
        return self.A

    @property
    def old_gradient(self):
        return self.old_grad

    def _elapsed(self):
        return self.elapsed_before+time.perf_counter()-self.started

    def _check(self):
        check = self.config.get('budget_check')
        if check is not None:
            check()
        if self._elapsed() >= float(self.config.get('max_wall_seconds', np.inf)):
            self.status = 'WALL_BUDGET'
            raise TimeoutError('R14 inherited wall budget exhausted')

    def _enter_stage(self, index):
        self._check()
        self.stage = int(index)
        stage = self.stages[self.stage]
        self.beta1 = float(stage.get('beta_material', stage.get('beta', 1.)))
        self.beta0 = float(stage.get('beta_current', self.beta1))
        if self.beta0 <= 0 or self.beta1 <= 0 or np.any(1+self.beta0*self.lo <= 0) or np.any(1+self.beta1*self.lo <= 0):
            raise ValueError('CIE requires positive beta on the original pole-free branch')
        self.U = NativeFourierAPCS(self.shape, stage['mf'], self.Vd, self.cost)
        self.basis_info = dict(mf=int(stage['mf']), fourier_columns=self.U.k,
            apcs_rank=None, native_action_without_roundoff_SVD=True)
        self.frequencies = [tuple(map(int, f)) for f in self.U.frequencies]
        self.A = np.array([self.old_coefficients.get(k, np.zeros(self.ns, complex)) for k in self.frequencies])
        if index == 0 and self.initial_physical_j is not None:
            UU = self.U.apply(np.eye(self.U.k, dtype=complex))
            self.A, _, warm_rank, _ = np.linalg.lstsq(UU, self.initial_physical_j-self.JD, rcond=1e-12)
            self.cost['initial_current_projection_SVD'] += 1
            self.cost['initial_projection_columns'] += self.U.k
            self.cost['initial_projection_rank'] = int(warm_rank)
        self.cost['peak_retained_FFT_spectrum_bytes'] = max(self.cost['peak_retained_FFT_spectrum_bytes'], self.green.spectrum.nbytes)
        self.cost['peak_native_coefficient_bytes'] = max(self.cost['peak_native_coefficient_bytes'], self.U.k*self.ns*16)
        self.cost['explicit_IDFT_columns_allocated'] = 0
        self.cost['explicit_GU_columns_allocated'] = 0
        self.cost['explicit_H_columns_allocated'] = 0
        self.J = self.JD+self.U.apply(self.A)
        self.Q = self.E+self.GJD+self.green.apply(self.J-self.JD)
        self.old_grad = self.old_direction = None
        self.iteration = 0
        self.stage_pending = False
        self.stage_stop = 'MAXITER'
        self._check()

    def _P(self, m, rhs, family):
        if family == '0':
            p = self.beta0/(1+self.beta0*m)
            return p[:, None]*rhs
        return self.blocks.apply(m, self.beta1, rhs)

    def _PH(self, m, rhs, family):
        if family == '0':
            p = self.beta0/(1+self.beta0*m)
            return p.conj()[:, None]*rhs
        _, d, a = self.blocks.parameters(m, self.beta1)
        p = self.beta1/d
        self.cost['scalar_complex_divisions'] += self.n
        self.cost['onehop_current_adjoint_p_scalar_divisions'] += self.n
        return p.conj()[:, None]*(rhs+self.blocks.off_action(a.conj()[:, None]*rhs, adjoint=True))

    def _SU(self, a):
        self.cost['S_columns'] += self.ns
        return self.S@self.U.apply(a)

    def _SUH(self, z):
        self.cost['S_adjoint_columns'] += self.ns
        return self.U.adjoint(self.S.conj().T@z)

    def _H(self, a, m):
        ua = self.U.apply(a)
        return self._P(m, ua-m[:, None]*self.green.apply(ua), self.pair[0])

    def _HH(self, z, m):
        v = self._PH(m, z, self.pair[0])
        return self.U.adjoint(v-self.green.adjoint(m[:, None]*v))

    def value_current(self, m=None, J=None, Q=None, *, include_prior=False):
        """D + rho/2 ||Pc(J-mQ)||²; paid only for source/block products.

        Default arguments are cached states.  Supplying J without Q is rejected
        rather than silently paying another Green action or physical solve.
        """
        m = self.m if m is None else np.asarray(m, float)
        if J is not None and Q is None:
            raise ValueError('Supply cached J and Q together')
        J, Q = (self.J, self.Q) if J is None else (np.asarray(J), np.asarray(Q))
        D = self.S@J-self.Y
        self.cost['S_columns'] += J.shape[1]
        R = self._P(m, J-m[:, None]*Q, self.pair[0])
        self.cost['diagnostic_current_value_evaluations'] += 1
        value = .5*_inner(D, D, self.Od, self.cost)+.5*self.rho*_inner(R, R, self.Os, self.cost)
        return value+(float(self.reg(m)[0]) if include_prior else 0.)

    def value_material(self, m=None, J=None, Q=None):
        """Fixed-current material state loss plus original prior, full Pm(m)."""
        m = self.m if m is None else np.asarray(m, float)
        if J is not None and Q is None:
            raise ValueError('Supply cached J and Q together')
        J, Q = (self.J, self.Q) if J is None else (np.asarray(J), np.asarray(Q))
        R = self._P(m, J-m[:, None]*Q, self.pair[1])
        self.cost['diagnostic_material_value_evaluations'] += 1
        return .5*self.rho*_inner(R, R, self.Os, self.cost)+float(self.reg(m)[0])

    def material_value_gradient(self, m=None, J=None, Q=None):
        """Optional exact derivative; never used as an extra optimization step."""
        m = self.m if m is None else np.asarray(m, float)
        if J is not None and Q is None:
            raise ValueError('Supply cached J and Q together')
        J, Q = (self.J, self.Q) if J is None else (np.asarray(J), np.asarray(Q))
        if self.pair[1] == '1':
            return onehop_material_value_gradient(self.blocks, J, Q, m, self.beta1,
                self.Os, self.rho, self.reg, self.cost)
        b = self.beta0
        d = 1+b*m
        R = b*(J-m[:, None]*Q)/d[:, None]
        dR = -b*(Q+b*J)/d[:, None]**2
        W = R@self.Os
        prior, gradient = self.reg(m)
        self.cost['diagnostic_scalar_material_gradient_evaluations'] += 1
        self.cost['source_metric_products'] += 1
        return (.5*self.rho*float(np.vdot(R, W).real)+float(prior),
            self.rho*np.real(np.sum(dR.conj()*W, axis=1))+np.asarray(gradient, float))

    def _finish_stage(self, reason):
        if not np.all(np.isfinite(self.m)):
            raise FloatingPointError('Nonfinite real material readout')
        self.stage_stop = reason
        self.old_coefficients = {k: a.copy() for k, a in zip(self.frequencies, self.A)}
        state = dict(stage=self.stage, m=self.m.copy(), J=self.J.copy(), Q=self.Q.copy(),
            coefficients=self.A.copy(), beta_current=self.beta0 if self.pair[0] == '0' else self.beta1,
            beta_material=self.beta0 if self.pair[1] == '0' else self.beta1,
            receiver_rank=self.receiver_rank, basis_info=dict(self.basis_info), pair=self.pair,
            global_iteration=self.iterations, iteration=self.iteration, stage_stop=reason)
        self.stage_states.append(state)
        self.history.append(dict(stage=self.stage, event='stage_end', status=reason,
            **self.basis_info, beta_current=state['beta_current'], beta_material=state['beta_material'],
            pair=self.pair, wall=self._elapsed()))
        self.stage_pending = self.stage+1 < len(self.stages)
        self.done = not self.stage_pending
        callback = self.config.get('stage_callback')
        if callback is not None:
            callback(deepcopy(state))

    def step(self):
        """Advance exactly one native PR/current + original material readout."""
        self._check()
        if self.done:
            raise StopIteration(self.status)
        if self.iterations >= int(self.config.get('max_iterations', 10**12)):
            self.status = 'ITERATION_BUDGET'
            raise StopIteration(self.status)
        if self.stage_pending:
            self._enter_stage(self.stage+1)
        while int(self.stages[self.stage].get('maxiter', 500)) <= 0:
            self._finish_stage('MAXITER')
            if self.done:
                raise StopIteration(self.status)
            self._enter_stage(self.stage+1)
        self._check()
        m_before = self.m.copy()
        previous_gradient = None if self.old_grad is None else self.old_grad.copy()
        previous_direction = None if self.old_direction is None else self.old_direction.copy()
        T = self._P(self.m, self.JD-self.m[:, None]*(self.E+self.GJD), self.pair[0])
        D, R = self.SJD-self.Y+self._SU(self.A), T+self._H(self.A, self.m)
        self.cost['scalar_current_operator_builds' if self.pair[0] == '0' else 'onehop_current_operator_builds'] += 1
        self.cost['current_objective_evaluations'] += 1
        data_before = .5*_inner(D, D, self.Od, self.cost)
        state_before = .5*self.rho*_inner(R, R, self.Os, self.cost)
        grad = self._SUH(D@self.Od)+self.rho*self._HH(R@self.Os, self.m)
        self.cost['source_metric_products'] += 2
        self.cost['current_gradient_evaluations'] += 1
        axis = 0 if self.sourcewise else None
        if self.old_grad is None:
            direction, pr = -grad, 0.
        else:
            pr = np.maximum(0., np.sum(np.real(grad.conj()*(grad-self.old_grad)), axis=axis)/
                np.maximum(np.sum(abs(self.old_grad)**2, axis=axis), 1e-300))
            direction = -grad+self.old_direction*pr
        gd = np.sum(np.real(grad.conj()*direction), axis=axis)
        restarted = bool(np.any(gd >= 0))
        if self.sourcewise:
            bad = gd >= 0
            direction[:, bad] = -grad[:, bad]
            self.cost['non_descent_source_restarts'] += int(np.count_nonzero(bad))
        elif gd >= 0:
            direction = -grad
            self.cost['non_descent_global_restarts'] += 1
        DD, RR = self._SU(direction), self._H(direction, self.m)
        numerator = -np.sum(np.real(grad.conj()*direction), axis=axis)
        if self.sourcewise:
            denom = self.wd*np.sum(abs(DD)**2, axis=0)+self.rho*self.ws*np.sum(abs(RR)**2, axis=0)
        else:
            denom = _inner(DD, DD, self.Od, self.cost)+self.rho*_inner(RR, RR, self.Os, self.cost)
        alpha = np.divide(numerator, denom, out=np.zeros_like(numerator), where=denom > 1e-300)
        change = direction*alpha
        if self.sourcewise:
            delta = float(np.mean(np.sum(abs(change)**2, axis=0)/np.maximum(np.sum(abs(self.A)**2, axis=0), 1e-300)))
        else:
            delta = float(np.vdot(change, change).real/max(np.vdot(self.A, self.A).real, 1e-300))
        self.A = self.A+change
        fluctuation = self.U.apply(self.A)
        self.J, self.Q = self.JD+fluctuation, self.E+self.GJD+self.green.apply(fluctuation)
        Dafter = self.SJD-self.Y+self.S@fluctuation
        self.cost['S_columns'] += self.ns
        Rafter = T+self._H(self.A, self.m)
        data_after = .5*_inner(Dafter, Dafter, self.Od, self.cost)
        state_after = .5*self.rho*_inner(Rafter, Rafter, self.Os, self.cost)
        if self.pair[1] == '1':
            self.m, inner = mobius_material_readout(self.blocks, self.J, self.Q, self.m,
                self.beta1, self.Os, self.Ls, self.rho, (self.lo, self.hi), self.reg,
                self.config, self.cost)
        elif self.sourcewise:
            self.m, inner = _baseline_readout(self.J, self.Q, self.m, self.beta0, self.ws,
                self.rho, (self.lo, self.hi), self.reg, self.has_reg, self.config, self.cost)
        else:
            self.m, inner = _baseline_readout(self.J@self.Ls, self.Q@self.Ls, self.m,
                self.beta0, np.ones(self.ns), self.rho, (self.lo, self.hi), self.reg,
                self.has_reg, self.config, self.cost)
        if not self.sourcewise or self.pair[1] == '1':
            self.cost['material_readout_source_whitening_products'] += 2
            self.cost['source_metric_products'] += 2
        self.old_grad, self.old_direction = grad, direction
        data_final = .5*_inner(Dafter, Dafter, self.Od, self.cost)
        def state_value(family):
            e = self._P(self.m, self.J-self.m[:, None]*self.Q, family)
            return .5*self.rho*_inner(e, e, self.Os, self.cost)
        self.iteration += 1
        self.iterations += 1
        row = dict(stage=self.stage, iteration=self.iteration, global_iteration=self.iterations,
            beta_current=self.beta0 if self.pair[0] == '0' else self.beta1,
            beta_material=self.beta0 if self.pair[1] == '0' else self.beta1,
            stage_beta_current=self.beta0, stage_beta_material=self.beta1,
            receiver_rank=self.receiver_rank, **self.basis_info,
            coefficient_relative_change_squared=delta, data_term=data_final,
            state_term_current=state_value(self.pair[0]), state_term_material=state_value(self.pair[1]),
            regularizer=float(self.reg(self.m)[0]),
            current_transform='scalar' if self.pair[0] == '0' else 'onehop',
            material_transform='scalar' if self.pair[1] == '0' else 'onehop', pair=self.pair,
            wall=self._elapsed(), sourcewise_cg=self.sourcewise, current_solver='pr',
            PR_plus=np.asarray(pr).tolist(), line_alpha=np.asarray(alpha).tolist(),
            direction_restarted=restarted, data_before_current=data_before, state_before_current=state_before,
            data_after_current=data_after, state_after_current=state_after,
            current_step_objective_change=data_after+state_after-data_before-state_before,
            readout_material_change=float(np.linalg.norm(self.m-m_before)), **inner)
        self.history.append(row.copy())
        self.cost['PR_CG_steps'] += 1
        self.cost['matrix_free_current_steps'] += 1
        row.update(m=self.m.copy(), m_before=m_before, J=self.J.copy(), Q=self.Q.copy(),
            coefficients=self.A.copy(), gradient=grad.copy(), old_grad=grad.copy(),
            direction=direction.copy(), previous_gradient=previous_gradient,
            previous_direction=previous_direction, alpha=np.asarray(alpha).copy(),
            PR=np.asarray(pr).copy(), stage_end=False, stage_stop=None)
        reason = None
        if self.iteration >= int(self.config.get('miniter', 2)) and delta < float(self.config.get('coefficient_tol', 1e-3)):
            reason = 'COEFFICIENT_TOL'
        elif self.iteration >= int(self.stages[self.stage].get('maxiter', 500)):
            reason = 'MAXITER'
        if reason is not None:
            self._finish_stage(reason)
            row.update(stage_end=True, stage_stop=reason)
        row['operation_counts'] = dict(self.cost)
        callback = self.config.get('iteration_callback')
        if callback is not None:
            callback(deepcopy(row) if self.config.get('deepcopy_iteration_callback', False) else dict(row))
        self._check()
        return row

    def switch_pair(self, pair, reset_memory=False):
        """Change only residual families; optionally reset the inherited PR pair."""
        self._validate_pair(pair)
        self.pair = pair
        if reset_memory:
            self.old_grad = self.old_direction = None
        return self

    def clone(self, pair=None, reset_memory=False):
        """Detached same-state fork using the same supplied physical model."""
        return self.from_state_dict(self.model, self.Y, self.config,
            self.state_dict(), pair=pair, reset_memory=reset_memory)

    def state_dict(self):
        """Detached complete numerical state, with no callback/model references."""
        fields = ('m', 'A', 'J', 'Q', 'Vd', 'JD', 'GJD', 'SJD', 'Od', 'Ld', 'Os', 'Ls',
            'old_grad', 'old_direction', 'old_coefficients', 'receiver_left',
            'receiver_singular_values', 'receiver_rank', 'initial_physical_j',
            'iterations', 'stage', 'iteration', 'stage_pending', 'done', 'status',
            'stage_stop', 'history', 'stage_states', 'frequencies', 'basis_info', 'beta0', 'beta1')
        result = {key: deepcopy(getattr(self, key)) for key in fields}
        result.update(checkpoint_version=self.checkpoint_version, pair=self.pair,
            operation_counts=dict(self.cost), frozen_config=_plain_config(self.config),
            elapsed_wall_seconds=self._elapsed(), shape=self.shape, n=self.n, ns=self.ns,
            observations=self.Y.copy(), incident_field=self.E.copy(), receiver_operator=self.S.copy(),
            green_cache={key: deepcopy(getattr(self.green, key)) for key in ('shape', 'n', 'pad', 'spectrum')},
            onehop_cache={key: deepcopy(getattr(self.blocks, key)) for key in ('indices', 'goff', 'diag_g', 'n', 'nblocks')})
        return result

    def checkpoint(self, path):
        """Write inspectable NPZ with named arrays and a safe recursive JSON tree."""
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        _save_state_npz(path, self.state_dict())
        return path

    @classmethod
    def from_checkpoint(cls, model, Y, config, path, pair=None, reset_memory=False):
        """Restore cached operators and all PR/stage state without repeating SVD/FFT.

        The frozen numerical flags in the checkpoint take precedence.  Runtime
        budget callbacks, prior objects and regularizer must be supplied again.
        No scientific work is rerun; checkpoint I/O remains caller-paid time.
        """
        saved = _load_state_npz(path)
        return cls.from_state_dict(model, Y, config, saved, pair=pair, reset_memory=reset_memory)

    @classmethod
    def from_state_dict(cls, model, Y, config, saved, pair=None, reset_memory=False):
        """In-memory twin of checkpoint restore, useful for same-state controls."""
        if saved.get('checkpoint_version') != cls.checkpoint_version:
            raise ValueError('Unsupported R14 checkpoint version')
        cfg = dict(config)
        cfg.update(deepcopy(saved['frozen_config']))
        result = cls.__new__(cls)
        result._bind(model, Y, cfg)
        if tuple(saved['shape']) != result.shape or not np.array_equal(saved['observations'], result.Y):
            raise ValueError('Checkpoint grid/observations mismatch')
        if not np.array_equal(saved['incident_field'], result.E) or not np.array_equal(saved['receiver_operator'], result.S):
            raise ValueError('Checkpoint public incident/receiver operators mismatch')
        skip = {'checkpoint_version', 'pair', 'operation_counts', 'frozen_config', 'elapsed_wall_seconds',
            'shape', 'n', 'ns', 'observations', 'incident_field', 'receiver_operator', 'green_cache', 'onehop_cache'}
        for key, value in saved.items():
            if key not in skip:
                setattr(result, key, deepcopy(value))
        result.cost = Counter(saved['operation_counts'])
        result.config['_counter'] = result.cost
        result.green = DiskPulseToeplitzGreen.__new__(DiskPulseToeplitzGreen)
        result.green.__dict__.update(deepcopy(saved['green_cache']))
        result.green.cost = result.cost
        result.blocks = OneHopAbsorption.__new__(OneHopAbsorption)
        result.blocks.__dict__.update(deepcopy(saved['onehop_cache']))
        result.blocks.cost = result.cost
        result.U = NativeFourierAPCS(result.shape, result.stages[result.stage]['mf'], result.Vd, result.cost)
        if result.frequencies != [tuple(map(int, f)) for f in result.U.frequencies] or result.A.shape != (result.U.k, result.ns):
            raise ValueError('Checkpoint native Fourier mask/coefficients mismatch')
        result._bind_weights()
        result.elapsed_before = float(saved['elapsed_wall_seconds'])
        result.started = time.perf_counter()
        result.switch_pair(saved['pair'] if pair is None else pair, reset_memory=reset_memory)
        result._check()
        return result
