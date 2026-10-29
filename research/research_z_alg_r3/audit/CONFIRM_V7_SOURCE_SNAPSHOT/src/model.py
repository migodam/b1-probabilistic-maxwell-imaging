"""Counted cached-operator VIE, shared by every R3 method.

The analytic Jacobian uses min(n,nr) complex RHS at a retained trial LU.
All source fields share that LU. No numerical derivatives are used.
"""
from collections import Counter
import numpy as np
from scipy.linalg import lu_factor, lu_solve
from .physics import Model, ReducedKernel


class PublicModel(Model):
    def __init__(self, public):
        self.book = None
        self.counts, self.timings = Counter(), Counter()
        self.state = self.factor = None
        if all(key in public for key in ('G', 'E0', 'S')):
            self.G = np.asarray(public['G'], dtype=complex)
            self.E0 = np.asarray(public['E0'], dtype=complex)
            self.S = np.asarray(public['S'], dtype=complex)
            self.n, self.ns, self.nr = self.G.shape[0], self.E0.shape[1], self.S.shape[0]
            if self.G.shape != (self.n, self.n) or self.E0.shape[0] != self.n or self.S.shape[1] != self.n:
                raise ValueError('Inconsistent cached operator shapes')
            if not all(np.all(np.isfinite(a)) for a in (self.G, self.E0, self.S)):
                raise ValueError('Finite cached operators required')
            with self.operation('CACHED_OPERATOR_IMPORT', cached_operator_imports=1):
                pass
        else:
            geometry = public.get('geometry', public)
            super().__init__(centers=geometry['centers'], sources=geometry['sources'],
                             receivers=geometry['receivers'], radius=float(geometry.get('radius', .04)),
                             k=float(geometry.get('k', 2*np.pi)))

    def forward(self, m, label='trial', jac=False, independent=False):
        m = np.asarray(m, dtype=float)
        if m.shape != (self.n,) or not np.all(np.isfinite(m)):
            raise ValueError('Finite real material vector of shape(n,) required')
        self.state = self.factor = None
        with self.operation('FULL_FORWARD', rhs=self.ns, LU=1, new_states=1,
                            full_factorization_dimension=self.n, label=label):
            A = np.eye(self.n, dtype=complex)-m[:, None]*self.G
            b = m[:, None]*self.E0
            factor = lu_factor(A, check_finite=False)
            j = lu_solve(factor, b, check_finite=False)
            E, y = self.E0+self.G@j, self.S@j
            state = dict(m=m.copy(), j=j, E=E, y=y, L=A,
                         state_relative_residual=float(np.linalg.norm(A@j-b)/max(np.linalg.norm(b), 1e-300)))
        if not np.all(np.isfinite(y)):
            raise FloatingPointError('Nonfinite full forward')
        self.state, self.factor = state, factor
        if jac:
            self.jacobian(state)
        return state

    def jacobian(self, state=None):
        state = self.state if state is None else state
        if state is None:
            raise ValueError('No state')
        if 'J' in state:
            return state['J']
        if state is not self.state or self.factor is None:
            raise RuntimeError('Jacobian requires retained matching full-state LU')
        tangent = self.n < self.nr
        rhs = self.n if tangent else self.nr
        with self.operation('FULL_MATERIAL_J_TANGENT' if tangent else 'FULL_MATERIAL_J_ADJOINT',
                            rhs=rhs, material_J_builds=1,
                            tangent_RHS=rhs if tangent else 0, adjoint_RHS=0 if tangent else rhs):
            if tangent:
                transfer = self.S@lu_solve(self.factor, np.eye(self.n), check_finite=False)
            else:
                transfer = lu_solve(self.factor, self.S.T, trans=1, check_finite=False).T
            J = np.einsum('ri,is->sri', transfer, state['E']).reshape(self.ns*self.nr, self.n)
            state['Jcomplex'], state['J'] = J, self.pack_jac(J)
        return state['J']


class ObservationView:
    """Subset metadata delegates every operation to the common physical ledger."""
    def __init__(self, model, receivers, sources):
        self.parent = model
        self.nr, self.ns, self.n = len(receivers), len(sources), model.n
        self.counts, self.timings = model.counts, model.timings

    def operation(self, *args, **kwargs):
        return self.parent.operation(*args, **kwargs)


def subset_kernel(kernel, receiver, sources):
    """Pure restriction of exact kernel to one receiver and selected sources."""
    result = object.__new__(ReducedKernel)
    result.__dict__.update(kernel.__dict__)
    sources = np.asarray(sources, dtype=int)
    result.model = ObservationView(kernel.model, [receiver], sources)
    result.q = kernel.q[:, sources].copy()
    result.L = kernel.L[receiver:receiver+1].copy()
    result.y0 = kernel.y0[np.ix_([receiver], sources)].copy()
    return result
