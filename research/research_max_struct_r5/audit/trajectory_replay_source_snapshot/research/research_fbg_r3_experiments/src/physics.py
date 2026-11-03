"""Thin API over the frozen R4 disk-pulse VIE; no new forward approximation."""
from pathlib import Path
from collections import Counter
import sys
import time
import numpy as np
from scipy.linalg import lu_solve, cho_factor, cho_solve

REPO = Path(__file__).resolve().parents[3]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
from research.research_z_alg_r4_diversity.src.physics import (
    Model as FrozenModel, build_operators, Regularizer, material_basis)


class MaxwellModel:
    def __init__(self, public, shape=None):
        self._physics = FrozenModel(public)
        self.G, self.E, self.S = self._physics.G, self._physics.E0, self._physics.S
        self.n, self.ns, self.nr = self.G.shape[0], self.E.shape[1], self.S.shape[0]
        self.shape = tuple(shape if shape is not None else public.get('grid_shape',
                           (int(np.sqrt(self.n)), int(np.sqrt(self.n)))))
        self.cost = self._physics.counts
        self.counts = self.cost
        self.timings = self._physics.timings
        self.reg = Regularizer(self.shape)
        self._current = None

    @staticmethod
    def pack(y):
        return FrozenModel.pack(y)

    def state(self, m):
        m = np.asarray(m, float).reshape(self.n)
        if self._current is not None and np.array_equal(m, self._current['m']):
            return self._current
        state = self._physics.forward(m)
        state['Y'], state['lu'] = state['y'], self._physics.factor
        state['L'] = self._physics.identity - m[:, None] * self.G
        self._current = state
        return state

    def forward(self, m):
        return self.state(m)['Y']

    def jacobian(self, m):
        return self._physics.jacobian(self.state(m))

    def tangent(self, m, h):
        state = self.state(m)
        hs = np.asarray(h, float).reshape(-1, self.n)
        rhs = (hs[:, :, None] * state['Q'][None]).transpose(1, 0, 2).reshape(self.n, -1)
        with self._physics.operation('CURRENT_TANGENT', rhs=rhs.shape[1], tangent_RHS=rhs.shape[1]):
            z = lu_solve(state['lu'], rhs, check_finite=False)
        return z.reshape(self.n, len(hs), self.ns).transpose(1, 0, 2)

    def regularizer(self, m):
        value, gradient, weights = self.reg.evaluate(m)
        hessian = (self.reg.D.T @ self.reg.D.multiply(weights[:, None])).toarray()
        return value, gradient, hessian

    def regularizer_residual(self, m):
        d = self.reg.D @ m
        q = np.sqrt(d*d + self.reg.epsilon**2)
        a = np.sqrt(2*self.reg.smooth/self.reg.edges)
        b = np.sqrt(2*self.reg.tv/self.reg.edges)
        rr = np.r_[a*d, b*d/np.sqrt(q+self.reg.epsilon)]
        derivative = b*(1/np.sqrt(q+self.reg.epsilon)-
                        .5*d*d/(q*(q+self.reg.epsilon)**1.5))
        jj = np.vstack((a*self.reg.D.toarray(), self.reg.D.multiply(derivative[:, None]).toarray()))
        return rr, jj

    def physical_score(self, m, Y):
        state = self.state(m)
        scale = np.vdot(Y, Y).real
        data = .5*np.vdot(state['Y']-Y, state['Y']-Y).real/scale
        return float(data + self.reg.evaluate(m)[0])

    def data_value_gradient(self, m, Y, beta=None, reference=None):
        state = self.state(m)
        R, Q = state['Y']-Y, state['Q']
        scale = float(np.vdot(Y,Y).real)
        if beta is None:
            value = .5*np.vdot(R,R).real/scale
            with self._physics.operation('DATA_ADJOINT', rhs=self.ns, adjoint_RHS=self.ns):
                a = lu_solve(state['lu'], self.S.conj().T@R, trans=2, check_finite=False)
            gradient = np.real(np.sum(a.conj()*Q, axis=1))/scale
        else:
            ref = self.state(np.full(self.n,.1)) if reference is None else reference
            # Retained original state/factor is explicit; reference lookup can change cache.
            V = Q + float(beta)*state['j']
            factor = (np.vdot(ref['Q'],ref['Q']).real/self.ns)*(1+.1*float(beta))**2/scale
            with self._physics.operation('FEEDBACK_GRAM', Gram=1, Gram_factorizations=1):
                chol = cho_factor(V.conj().T@V, lower=True, check_finite=False)
                W = cho_solve(chol, R.conj().T, check_finite=False).conj().T
            value = .5*factor*np.vdot(R,W).real
            rhs = self.S.conj().T@W-(self.G.conj().T+float(beta)*np.eye(self.n))@(V@(W.conj().T@W))
            with self._physics.operation('FEEDBACK_ADJOINT', rhs=self.ns, adjoint_RHS=self.ns):
                a = lu_solve(state['lu'],rhs,trans=2,check_finite=False)
            gradient = factor*np.real(np.sum(a.conj()*Q,axis=1))
        reg, greg, _ = self.reg.evaluate(m)
        return float(value+reg), gradient+greg


def make_grid(side, layout='uniform', sources=6, receivers=16, k=2*np.pi):
    """Reassemble each physical grid, covering the same 0.96-square region."""
    side = int(side)
    spacing = .96/side
    axis = (np.arange(side)+.5)*spacing-.48
    yy,xx = np.meshgrid(axis,axis,indexing='ij')
    centers = np.c_[xx.ravel(),yy.ravel()]
    angles = 2*np.pi*np.arange(sources)/sources + .071
    if layout == 'jittered':
        angles = angles + .17*np.sin(2.31*np.arange(sources)+.4)
    source_positions = 1.4*np.c_[np.cos(angles),np.sin(angles)]
    recv_angles = 2*np.pi*np.arange(receivers)/receivers+.031
    receiver_positions = 1.1*np.c_[np.cos(recv_angles),np.sin(recv_angles)]
    radius = .45*spacing
    started = time.perf_counter()
    public = build_operators(centers,source_positions,receiver_positions,radius,float(k))
    public.update(E=public['E0'],grid_shape=np.array([side,side]), centers=centers,
                  sources=source_positions,receivers=receiver_positions,k=float(k),radius=radius,
                  green_setup_wall=time.perf_counter()-started)
    return public


def load_model(path):
    with np.load(path,allow_pickle=False) as f:
        public = {k:f[k] for k in f.files}
    if 'operators_path' in public:
        with np.load(Path(path).parent/str(public['operators_path'].item()),allow_pickle=False) as f:
            public = {**{k:f[k] for k in f.files},**public}
    if 'E0' not in public:
        public['E0'] = public['E']
    return MaxwellModel(public), public
