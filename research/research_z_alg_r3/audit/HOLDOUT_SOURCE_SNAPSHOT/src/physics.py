"""Independent pulse TM VIE and counted exact fixed-active resolvent kernel.

No observations or material ground truth are stored in this module.  All setup
choices depend only on the labeled geometry.  A supplied book may implement
``op(kind, rhs=0, **metadata)`` and optionally ``append``; every model also keeps
an in-memory counter, so uninstrumented calls cannot silently become free.
"""
from collections import Counter
from contextlib import contextmanager
from pathlib import Path
import json
import time

import numpy as np
from scipy.linalg import lu_factor, lu_solve
from scipy.special import hankel1, jv


class Model:
    def __init__(self, book=None, centers=None, sources=None, receivers=None,
                 radius=.04, k=2*np.pi):
        if centers is None or sources is None or receivers is None:
            raise ValueError('Explicit centers, sources, receivers required')
        self.book = book
        self.centers = np.asarray(centers, dtype=float)
        self.sources = np.asarray(sources, dtype=float)
        self.receivers = np.asarray(receivers, dtype=float)
        self.n, self.ns, self.nr = len(centers), len(sources), len(receivers)
        self.radius, self.k = float(radius), float(k)
        self.counts = Counter()
        self.timings = Counter()
        self.state = self.factor = None
        with self.operation('GREEN_SETUP', Green_kernel_evaluations=self.n**2+self.n*self.ns+self.nr*self.n):
            d = np.linalg.norm(self.centers[:, None]-self.centers[None, :], axis=2)
            if self.n > 1 and np.min(d+np.eye(self.n)*1e10) <= 2*radius:
                raise ValueError('Overlapping pulse disks')
            np.fill_diagonal(d, 1.)
            self.G = self.external(d)
            np.fill_diagonal(self.G, .5j*np.pi*k*radius*hankel1(1, k*radius)-1)
            self.E0 = .25j*hankel1(0, k*np.linalg.norm(self.centers[:, None]-self.sources[None, :], axis=2))
            self.S = self.external(np.linalg.norm(self.receivers[:, None]-self.centers[None, :], axis=2))
        with self.operation('RECEIVER_SVD', receiver_SVD=1):
            self.S_u, self.S_s, self.S_vh = np.linalg.svd(self.S, full_matrices=True)

    @contextmanager
    def operation(self, kind, rhs=0, **cost):
        self.counts[kind] += 1
        self.counts['actual_complex_RHS'] += rhs
        for name, value in cost.items():
            if isinstance(value, (int, float)) and not isinstance(value, bool):
                self.counts[name] += value
        start = time.perf_counter()
        try:
            if self.book is None:
                yield None
            else:
                with self.book.op(kind, rhs=rhs, **cost) as action:
                    yield action
        finally:
            self.timings[kind] += time.perf_counter()-start

    def external(self, d):
        return .5j*np.pi*self.k*self.radius*jv(1, self.k*self.radius)*hankel1(0, self.k*d)

    def forward(self, m, label='state', jac=False, independent=False):
        m = np.asarray(m, dtype=float)
        if m.shape != (self.n,) or not np.all(np.isfinite(m)):
            raise ValueError('Finite real material vector required')
        self.state = self.factor = None
        with self.operation('INDEPENDENT_FORWARD' if independent else 'FULL_FORWARD',
                            rhs=self.ns, LU=1, new_states=1,
                            Green_matvec_columns=self.ns, label=label):
            A = np.eye(self.n, dtype=complex)-m[:, None]*self.G
            b = m[:, None]*self.E0
            if independent:
                j, factor = np.linalg.solve(A, b), None
            else:
                factor = lu_factor(A, check_finite=False)
                j = lu_solve(factor, b, check_finite=False)
            E, y = self.E0+self.G@j, self.S@j
            singular_values = np.linalg.svd(A, compute_uv=False)
            state = dict(m=m.copy(), j=j, E=E, y=y, L=A,
                         state_relative_residual=float(np.linalg.norm(A@j-b)/max(np.linalg.norm(b), 1e-300)),
                         resolvent_sigma_min=float(singular_values[-1]))
        self.state, self.factor = state, factor
        if jac:
            if independent:
                raise ValueError('Independent validation deliberately provides no retained LU for J')
            self.jacobian(state)
        return state

    def jacobian(self, state=None):
        state = self.state if state is None else state
        if state is None:
            raise ValueError('No state')
        if 'J' in state:
            return state['J']
        if state is not self.state or self.factor is None:
            raise RuntimeError('Jacobian requires the sole live LU at the same state')
        with self.operation('FULL_MATERIAL_J_ADJOINT', rhs=self.nr,
                            adjoint_RHS=self.nr, material_J_builds=1):
            transfer = lu_solve(self.factor, self.S.T, trans=1, check_finite=False).T
            J = np.einsum('ri,is->sri', transfer, state['E']).reshape(self.ns*self.nr, self.n)
            state['Jcomplex'], state['J'] = J, self.pack_jac(J)
        return state['J']

    def pack(self, y):
        a = np.asarray(y).T
        return np.stack([a.real, a.imag], axis=-1).reshape(-1)

    def unpack(self, y):
        a = np.asarray(y, dtype=float).reshape(self.ns, self.nr, 2)
        return (a[..., 0]+1j*a[..., 1]).T

    def pack_jac(self, J):
        return np.stack([J.real, J.imag], axis=1).reshape(2*self.ns*self.nr, self.n)


class ReducedKernel:
    """Exact kernel for m=m0+U delta with fixed coordinate support U.

    U columns must be distinct coordinate vectors.  A generic dense material
    basis would require rank n updates; pretending its update rank is r would
    be mathematically incorrect.  Setup is charged n-source RHS plus r setup
    RHS.  Evaluations are reduced algebraic solves, never full PDE forwards.
    """
    def __init__(self, model, m0, U=None, active=None):
        self.model = model
        self.m0 = np.asarray(m0, dtype=float).copy()
        if U is not None and active is not None:
            raise ValueError('Supply U or active, not both')
        if U is None:
            active = np.arange(model.n) if active is None else np.asarray(active, dtype=int)
            U = np.eye(model.n)[:, active]
        U = np.asarray(U, dtype=float)
        if U.ndim != 2 or U.shape[0] != model.n:
            raise ValueError('U shape must be (n,r)')
        active = np.argmax(np.abs(U), axis=0)
        if len(set(active.tolist())) != len(active) or not np.array_equal(U, np.eye(model.n)[:, active]):
            raise ValueError('Exact rank-r kernel requires distinct unscaled coordinate selector columns')
        self.U, self.active, self.r = U, active, U.shape[1]
        self.M0 = np.eye(model.n, dtype=complex)-self.m0[:, None]*model.G
        with model.operation('FIXED_ACTIVE_SETUP', rhs=model.ns+self.r, LU=1,
                             source_setup_RHS=model.ns, active_setup_RHS=self.r,
                             fixed_active_setups=1):
            factor = lu_factor(self.M0, check_finite=False)
            self.j0 = lu_solve(factor, self.m0[:, None]*model.E0, check_finite=False)
            self.V = lu_solve(factor, U, check_finite=False)
            self.K = U.T@model.G@self.V
            self.L = model.S@self.V
            self.q = U.T@(model.E0+model.G@self.j0)
            self.y0 = model.S@self.j0

    def coefficients(self, delta):
        delta = np.asarray(delta, dtype=float)
        if delta.shape != (self.r,):
            raise ValueError('delta shape differs from active support')
        return np.eye(self.r, dtype=complex)-delta[:, None]*self.K, delta[:, None]*self.q

    def reduced_current(self, delta):
        A, b = self.coefficients(delta)
        with self.model.operation('REDUCED_RATIONAL_EVALUATION', reduced_LU=1,
                                  reduced_complex_RHS=self.model.ns,
                                  reduced_dimension=self.r):
            return np.linalg.solve(A, b)

    def evaluate(self, delta):
        return self.y0+self.L@self.reduced_current(delta)

    __call__ = evaluate

    def state(self, delta):
        w = self.reduced_current(delta)
        j = self.j0+self.V@w
        m = self.m0+self.U@np.asarray(delta, dtype=float)
        E = self.model.E0+self.model.G@j
        return dict(m=m, j=j, E=E, y=self.y0+self.L@w,
                    reduced_w=w, validation='ALGEBRAIC_ONLY_REQUIRES_INDEPENDENT_FORWARD')

    def jacobian(self, delta):
        """Real response derivative, using only the small fixed-active solve."""
        A, b = self.coefficients(delta)
        with self.model.operation('REDUCED_RATIONAL_J', reduced_LU=1,
                                  reduced_complex_RHS=self.model.ns+self.r,
                                  reduced_J_builds=1):
            factor = lu_factor(A, check_finite=False)
            w = lu_solve(factor, b, check_finite=False)
            transfer = lu_solve(factor, self.L.T, trans=1, check_finite=False).T
            field = self.q+self.K@w
            J = np.einsum('ri,is->sri', transfer, field).reshape(self.model.ns*self.model.nr, self.r)
        return np.stack([J.real, J.imag], axis=1).reshape(2*self.model.ns*self.model.nr, self.r)

    def multiaffine_numerators(self, target, halfwidth=None):
        """Coefficients c[mask,real-channel] of cleared response equations.

        Variables are delta_i/halfwidth_i.  Every monomial is square free.
        Coefficients use row multilinearity of an augmented determinant, not
        fits, derivative linearization, or interpolation.  Returned denominator
        coefficients permit explicit rejection of denominator-zero roots.
        """
        target = np.asarray(target, dtype=complex)
        if target.shape != self.y0.shape:
            raise ValueError('Target must have (receivers,sources) shape')
        h = np.ones(self.r) if halfwidth is None else np.asarray(halfwidth, dtype=float)
        if h.shape != (self.r,) or np.any(h <= 0):
            raise ValueError('Positive halfwidth vector required')
        nc = self.model.nr*self.model.ns
        coeff = np.empty((1 << self.r, nc), dtype=complex)
        denominator = np.empty(1 << self.r, dtype=complex)
        residual = (self.y0-target).T.reshape(-1)
        with self.model.operation('MULTIAFFINE_COEFFICIENTS',
                                  algebraic_determinants=(1 << self.r)*(nc+1),
                                  coefficient_arrays=1):
            for mask in range(1 << self.r):
                rows = np.array([(mask >> i) & 1 for i in range(self.r)], dtype=bool)
                A = np.eye(self.r, dtype=complex)
                A[rows] = -h[rows, None]*self.K[rows]
                denominator[mask] = np.linalg.det(A)
                for channel in range(nc):
                    source, receiver = divmod(channel, self.model.nr)
                    aug = np.zeros((self.r+1, self.r+1), dtype=complex)
                    aug[:self.r, :self.r] = A
                    aug[:self.r, -1][rows] = h[rows]*self.q[rows, source]
                    aug[-1, :self.r] = -self.L[receiver]
                    aug[-1, -1] = residual[channel]
                    coeff[mask, channel] = np.linalg.det(aug)
        real = np.stack([coeff.real, coeff.imag], axis=2).reshape(1 << self.r, 2*nc)
        return dict(coefficients=real, complex_coefficients=coeff,
                    denominator_coefficients=denominator,
                    halfwidth=h, active=self.active.copy(),
                    equation='sum_mask coefficients[mask,c]*prod_i_in_mask x_i = 0',
                    variable='delta=halfwidth*x', completeness='REPRESENTATION_ONLY')

    def restrict(self, active, fixed_material, singular_floor=1e-12):
        """Exact Schur restriction with no new full-state or full-column RHS.

        ``active`` contains global material indices retained as variables.
        Complementary indices in this kernel's support are fixed to the full
        vector ``fixed_material``.  Retained variables keep their original m0
        reference.  Singular complementary blocks raise rather than silently
        omit a slice.  The returned kernel supports all original evaluations.
        """
        active = np.asarray(active, dtype=int)
        fixed_material = np.asarray(fixed_material, dtype=float)
        if fixed_material.shape != (self.model.n,) or len(set(active.tolist())) != len(active):
            raise ValueError('Distinct active indices and full fixed material vector required')
        index = {int(cell): i for i, cell in enumerate(self.active)}
        if any(int(cell) not in index for cell in active):
            raise ValueError('Restriction active coordinates must belong to parent support')
        A = np.asarray([index[int(cell)] for cell in active], dtype=int)
        B = np.asarray([i for i in range(self.r) if i not in set(A.tolist())], dtype=int)
        result = object.__new__(ReducedKernel)
        result.model, result.active, result.r = self.model, active.copy(), len(A)
        result.U = self.U[:, A].copy()
        result.m0 = self.m0.copy()
        result.M0 = None  # No full restricted reference matrix is formed/factored.
        result.parent_active = self.active.copy()
        result.restriction = dict(fixed_indices=self.active[B].copy(),
                                  fixed_material=fixed_material[self.active[B]].copy())
        with self.model.operation('SCHUR_SLICE_RESTRICTION', reduced_LU=int(len(B)>0),
                                  reduced_complex_RHS=self.model.ns+len(A) if len(B) else 0,
                                  reduced_dimension=len(B), schur_slice_restrictions=1,
                                  algebraic_SVD=int(len(B)>0)):
            if len(B):
                delta = fixed_material[self.active[B]]-self.m0[self.active[B]]
                TB = np.eye(len(B), dtype=complex)-delta[:, None]*self.K[np.ix_(B, B)]
                singular = np.linalg.svd(TB, compute_uv=False)
                gap = float(singular[-1]/singular[0])
                if gap <= singular_floor:
                    raise np.linalg.LinAlgError('Complementary Schur block singular/ill-conditioned')
                factor = lu_factor(TB, check_finite=False)
                b = lu_solve(factor, delta[:, None]*self.q[B], check_finite=False)
                transfer = lu_solve(factor, delta[:, None]*self.K[np.ix_(B, A)], check_finite=False)
                result.K = self.K[np.ix_(A, A)]+self.K[np.ix_(A, B)]@transfer
                result.q = self.q[A]+self.K[np.ix_(A, B)]@b
                result.L = self.L[:, A]+self.L[:, B]@transfer
                result.y0 = self.y0+self.L[:, B]@b
                result.j0 = self.j0+self.V[:, B]@b
                result.V = self.V[:, A]+self.V[:, B]@transfer
                result.m0[self.active[B]] = fixed_material[self.active[B]]
                result.restriction['relative_singular_gap'] = gap
            else:
                result.K = self.K[np.ix_(A, A)].copy()
                result.q, result.L = self.q[A].copy(), self.L[:, A].copy()
                result.y0, result.j0, result.V = self.y0.copy(), self.j0.copy(), self.V[:, A].copy()
                result.restriction['relative_singular_gap'] = 1.
        return result


def model_for(book=None, stage='E1'):
    contract_path = Path(__file__).resolve().parents[2]/'research_z_alg_asymmetric_validation'/'EXPERIMENT_CONTRACT.json'
    contract = json.loads(contract_path.read_text())
    c = contract['E1']
    centers = c['centers']
    if stage == 'E2':
        centers = [centers[i] for i in contract['E2']['centers_from_E1']]
    elif stage != 'E1':
        raise ValueError('Stage must be E1 or E2')
    return Model(book, centers, c['sources'], c['receivers'],
                 radius=contract['disk_radius'], k=contract['k'])
