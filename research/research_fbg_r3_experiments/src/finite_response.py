"""Exact fixed-outside physical-cell response, with explicit incremental costs.

No material labels or target geometry enter this module.  The kernel is the
finite Maxwell response, not a Born linearization.  Full-system LU ownership
stays with the supplied model; the local ledger counts the work requested here.
"""
from collections import Counter
from contextlib import contextmanager, nullcontext
import time
import numpy as np
from scipy.linalg import lu_solve, qr


class WorkLedger:
    def __init__(self):
        self.counts = Counter()
        self.events = []

    @contextmanager
    def operation(self, name, **counts):
        begin = time.perf_counter()
        self.counts.update(counts)
        try:
            yield
        finally:
            self.events.append(dict(operation=name, wall_seconds=time.perf_counter()-begin,
                                    counts=counts))

    def as_dict(self):
        return dict(counts=dict(self.counts), events=self.events)


def pack_complex(a):
    """Frozen parent convention: source, receiver, real/imaginary."""
    return np.stack([np.asarray(a).T.real, np.asarray(a).T.imag], axis=-1).ravel()


def pack_jacobian(a):
    # a: receiver, source, material.
    a = a.transpose(1, 0, 2)
    return np.stack([a.real, a.imag], axis=2).reshape(-1, a.shape[-1])


def model_operation(model, name, rhs=0, **counts):
    physics = getattr(model, '_physics', None)
    return physics.operation(name, rhs=rhs, **counts) if physics is not None else nullcontext()


def spatial_supports(shape, sizes=(4, 6), stride=1):
    """Public contiguous-cell dictionary; row-major order is the tie breaker."""
    ny, nx = map(int, shape)
    grid = np.arange(ny*nx).reshape(ny, nx)
    blocks = []
    for k in sizes:
        orientations = [(2, 2)] if k == 4 else [(2, 3), (3, 2)] if k == 6 else []
        for h, w in orientations:
            for y in range(0, ny-h+1, stride):
                for x in range(0, nx-w+1, stride):
                    blocks.append(grid[y:y+h, x:x+w].ravel().tolist())
    return blocks


class ExactFiniteResponse:
    """F(d)=F0+O solve(I-diag(d) C,diag(d) E_A), all original data."""
    def __init__(self, model, m0, support, state=None, ledger=None):
        self.model = model
        self.m0 = np.asarray(m0, float).copy()
        self.support = np.asarray(support, int)
        self.k = self.r = len(self.support)
        if self.k not in (4, 6) or len(set(self.support)) != self.k:
            raise ValueError('Support must contain exactly 4 or 6 distinct physical cells')
        if np.any(self.support < 0) or np.any(self.support >= model.n):
            raise ValueError('Support index outside full material grid')
        self.ledger = ledger or WorkLedger()
        if state is None:
            with self.ledger.operation('ANCHOR_FULL_STATE', full_state_requests=1):
                state = model.state(self.m0)
        self.state = state
        self.F0 = np.asarray(state['Y'])
        U = np.eye(model.n, dtype=complex)[:, self.support]
        with self.ledger.operation('FINITE_KERNEL_SETUP', full_support_RHS=self.k,
                                   green_dense_matrix_products=1, receiver_matrix_products=1):
            with model_operation(model,'N2_FINITE_KERNEL_SETUP',rhs=self.k,support_setup_RHS=self.k):
                self.V = lu_solve(state['lu'], U)
                self.C = (model.G@self.V)[self.support, :]
                self.E_A = state['Q'][self.support, :]
                self.O = model.S@self.V
        self.U = U.real

    def material(self, delta):
        m = self.m0.copy()
        m[self.support] += np.asarray(delta, float)
        return m

    def coefficients(self, delta):
        d = np.asarray(delta)
        return np.eye(self.k)-d[:, None]*self.C, d[:, None]*self.E_A

    def evaluate(self, delta, jacobian=False):
        d = np.asarray(delta, float)
        with self.ledger.operation('FINITE_ALL_DATA', small_LU=1,
                                   small_source_RHS=self.E_A.shape[1],
                                   finite_output_products=1):
            A, rhs = self.coefficients(d)
            # One retained small factorization for forward and optional derivative.
            from scipy.linalg import lu_factor
            lu = lu_factor(A)
            eta = lu_solve(lu, rhs)
            out = self.F0+self.O@eta
        if not jacobian:
            return out
        with self.ledger.operation('FINITE_ALL_DATA_JACOBIAN', small_receiver_RHS=self.O.shape[0]):
            Z = lu_solve(lu, self.O.T, trans=1).T
            q = self.E_A+self.C@eta
            J = np.einsum('ri,ip->rpi', Z, q)
        return out, pack_jacobian(J)

    def select_channels(self, count=1):
        s = self.k//2
        if self.E_A.shape[1] < s:
            raise ValueError('INSUFFICIENT_EXISTING_SOURCES')
        with self.ledger.operation('ORIGINAL_CHANNEL_SELECTION', pivoted_QR=1):
            _, R, pivots = qr(self.E_A, mode='economic', pivoting=True)
            diag = np.abs(np.diag(R))
            rank = int(np.sum(diag > max(diag[0], 1e-300)*1e-12))
            if rank < s:
                raise ValueError('SOURCE_SUBKERNEL_RANK_DEFICIENT')
            sources = np.asarray(pivots[:s], int)
            receivers = np.argsort(-np.linalg.norm(self.O, axis=1), kind='stable')[:count]
        return sources, receivers

    def quotient_view(self, sources, receiver):
        """Frozen Z-ALG adapter, exposing only the selected original channels."""
        from types import SimpleNamespace
        return SimpleNamespace(r=self.k, q=self.E_A[:, sources], K=self.C,
                               L=self.O[receiver:receiver+1],
                               y0=self.F0[receiver:receiver+1, sources],
                               model=SimpleNamespace(ns=self.k//2, nr=1,
                                                     operation=self.ledger.operation))


def full_data_jacobian(model, state, ledger):
    """Analytic full-material Jacobian using the already retained physical LU."""
    with ledger.operation('FULL_MATERIAL_JACOBIAN', full_receiver_transpose_RHS=model.S.shape[0]):
        with model_operation(model,'N2_FULL_MATERIAL_JACOBIAN',rhs=model.S.shape[0],receiver_transpose_RHS=model.S.shape[0]):
            Z = lu_solve(state['lu'], model.S.T, trans=1).T
            return pack_jacobian(np.einsum('ri,ip->rpi', Z, state['Q']))
