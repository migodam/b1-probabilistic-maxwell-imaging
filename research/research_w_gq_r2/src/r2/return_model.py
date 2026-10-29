"""Return-only compression with exact discrete free-target scattering.

The public ports are FE load/field coefficient coordinates on the *whole*
target ROI.  This is neither an angular-spectrum port nor an L2-whitened
operator.  A deterministic QR of 32 Gaussian range actions is an elementary
randomized approximation, with no oversampling or operator-norm certificate.

All inverse actions are supplied by ``r2.physics.Eval.factor``.  Its factors
reserve every RHS before solving; this module must not reserve them twice.
No full restricted free Green matrix C0, target mask port space, dense m by m
inverse, pseudoinverse, Born replacement, or nodal-incident load is formed.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np


def _columns(value, rows, name):
    a = np.asarray(value, dtype=complex)
    if a.ndim == 1:
        a = a[:, None]
    if (a.ndim != 2 or a.shape[0] != rows or not a.shape[1]
            or not np.isfinite(a).all()):
        raise ValueError(f'{name} must have finite shape ({rows}, columns)')
    return a


def ratio(numerator, denominator):
    """No denominator jitter; a nonzero error against zero has no finite ratio."""
    numerator, denominator = float(numerator), float(denominator)
    return numerator/denominator if denominator else (0. if not numerator else None)


def relative_error(value, reference):
    return ratio(np.linalg.norm(np.asarray(value)-reference), np.linalg.norm(reference))


def linear_defect(matrix, solution, rhs):
    """Informational residual; a small backward error is not a field bound."""
    x, b = np.asarray(solution), np.asarray(rhs)
    residual = matrix@x-b
    scale = np.linalg.norm(matrix.data if hasattr(matrix, 'tocsc') else matrix)
    return dict(absolute=float(np.linalg.norm(residual)),
                rhs_relative=ratio(np.linalg.norm(residual), np.linalg.norm(b)),
                backward=ratio(np.linalg.norm(residual),
                               scale*np.linalg.norm(x)+np.linalg.norm(b)))


@dataclass(frozen=True)
class PublicPorts:
    indices: np.ndarray
    nfree: int

    def __post_init__(self):
        raw = np.asarray(self.indices)
        if (raw.ndim != 1 or not np.issubdtype(raw.dtype, np.integer)
                or not len(raw) or np.any(raw < 0) or np.any(raw >= self.nfree)
                or np.any(np.diff(raw) <= 0)):
            raise ValueError('ports must be increasing public-ROI free-DOF indices')
        owned = raw.astype(np.int64, copy=True)
        owned.flags.writeable = False
        object.__setattr__(self, 'indices', owned)

    @property
    def size(self):
        return len(self.indices)

    def inject(self, value):
        a = _columns(value, self.size, 'port coefficients')
        result = np.zeros((self.nfree, a.shape[1]), complex, order='F')
        result[self.indices] = a
        return result

    def restrict(self, value):
        return _columns(value, self.nfree, 'free coefficients')[self.indices].copy()

    def target_block(self, matrix, forcing):
        """Reject target terms outside the frozen public ROI, including columns."""
        if matrix.shape != (self.nfree, self.nfree):
            raise ValueError('target matrix must use the common free FE space')
        force = _columns(forcing, self.nfree, 'analytic target forcing')
        inside = np.zeros(self.nfree, bool)
        inside[self.indices] = True
        entries = matrix.tocoo()
        nonzero = entries.data != 0.
        if (np.any(~inside[entries.row[nonzero]])
                or np.any(~inside[entries.col[nonzero]])
                or np.any(force[~inside] != 0.)):
            raise ValueError('target matrix/analytic forcing leaves public target ROI')
        return matrix[self.indices][:, self.indices].tocsc(), force[self.indices].copy()


def deterministic_qr(actions):
    """Keep the Gaussian column order and fix each complex QR phase."""
    q, r = np.linalg.qr(np.asarray(actions, complex), mode='reduced')
    diagonal = np.diag(r)
    phase = np.ones(len(diagonal), complex)
    nonzero = np.abs(diagonal) != 0.
    phase[nonzero] = diagonal[nonzero]/np.abs(diagonal[nonzero])
    return q*phase[None, :], phase.conj()[:, None]*r


@dataclass(frozen=True)
class ReturnApproximation:
    Q: np.ndarray
    Hcore: np.ndarray
    Ccore: np.ndarray
    raw: dict[str, np.ndarray]
    metadata: dict[str, Any]

    def rank(self, rank):
        if (isinstance(rank, (bool, np.bool_)) or int(rank) != rank
                or not 1 <= rank <= self.Q.shape[1]):
            raise ValueError('rank must fit the predeclared range basis')
        r = int(rank)
        return self.Q[:, :r], self.Hcore[:r, :r], self.Ccore[:r, :r]


def build_return_approximation(free_factor, wall_factor, ports, *,
                               columns=32, seed=28102931, label='E3'):
    """Pay 4*columns RHS: H Omega and H Q, with separate Gw/G0 actions.

    H=P(Gw-G0)P^T.  The separate Cw Q action also supplies the whole-C
    diagnostic in exactly the same Q; it is not a range finder for Cw.
    The same seed is reset for each wall, so all walls get the same Omega.
    """
    if (isinstance(columns, (bool, np.bool_)) or int(columns) != columns
            or not 1 <= columns <= ports.size):
        raise ValueError('Gaussian column count must fit the public port space')
    columns = int(columns)
    rng = np.random.default_rng(seed)
    omega = (rng.standard_normal((ports.size, columns))
             + 1j*rng.standard_normal((ports.size, columns)))/np.sqrt(2.)
    load = ports.inject(omega)
    wall_omega = ports.restrict(wall_factor.solve(
        load, label=label+' Gw Omega', purpose='randomized return range setup'))
    free_omega = ports.restrict(free_factor.solve(
        load, label=label+' G0 Omega', purpose='randomized return range setup'))
    del load
    actions = wall_omega-free_omega
    q, qr_r = deterministic_qr(actions)
    load = ports.inject(q)
    wall_q = ports.restrict(wall_factor.solve(
        load, label=label+' Gw Q', purpose='projected return and same-basis whole-C core'))
    free_q = ports.restrict(free_factor.solve(
        load, label=label+' G0 Q', purpose='projected return core'))
    del load
    hq = wall_q-free_q
    hcore, ccore = q.conj().T@hq, q.conj().T@wall_q
    range_singular = np.linalg.svd(actions, compute_uv=False)
    raw = dict(Omega=omega, Q=q, QR_R=qr_r, Gw_Omega_ports=wall_omega,
               G0_Omega_ports=free_omega, H_Omega=actions,
               Gw_Q_ports=wall_q, G0_Q_ports=free_q, H_Q=hq,
               Hcore=hcore, Ccore=ccore,
               singular_values_H_Omega=range_singular,
               singular_values_H_Q=np.linalg.svd(hq, compute_uv=False),
               singular_values_Hcore=np.linalg.svd(hcore, compute_uv=False),
               singular_values_Cw_Q=np.linalg.svd(wall_q, compute_uv=False),
               singular_values_Ccore=np.linalg.svd(ccore, compute_uv=False))
    metadata = dict(seed=int(seed), range_columns=columns, oversampling=0,
                    basis='deterministic complex QR of H Omega; original Gaussian column order',
                    metric='Euclidean FE load/field coefficient coordinates; no L2 whitening',
                    return_operator='P (Gw-G0) P^T on the whole public target ROI',
                    approximation='Q_r (Q_r^* H Q_r) Q_r^*; nested QR columns',
                    same_basis_whole_C='Q_r (Q_r^* Cw Q_r) Q_r^*; basis fitted to H only',
                    operator_norm_certificate=False, physical_L2_spectrum=False,
                    periodic_angular_spectrum_claim=False, oversampling_guarantee=False,
                    singular_values_scope='sampled action matrices and projected cores only',
                    setup_rhs=4*columns,
                    orthogonality_defect=float(np.linalg.norm(q.conj().T@q-np.eye(columns))),
                    qr_reconstruction_error=relative_error(q@qr_r, actions))
    return ReturnApproximation(q, hcore, ccore, raw, metadata)


@dataclass(frozen=True)
class TargetSelfTemplate:
    q: np.ndarray
    j0: np.ndarray
    KQ: np.ndarray
    VQ: np.ndarray
    Vports: Any
    bXports: np.ndarray
    raw: dict[str, np.ndarray]
    metadata: dict[str, Any]


def build_target_self(ports, assembly, free_target_factor, Q, wall_u, *, label='E3'):
    """Exact (L0-V)^-1 target self action; never project C0 into Q.

    q=bX+V uw; j0=q+V(L0-V)^-1 q and
    KQ=VQ+V(L0-V)^-1 VQ, where injection/restriction are implicit.
    bX is the analytic-incident cut-integral supplied by the backend.
    """
    v, bx = ports.target_block(assembly.Vtarget, assembly.bX)
    q = bx+v@ports.restrict(wall_u)
    free_q = ports.restrict(free_target_factor.solve(
        ports.inject(q), label=label+' exact isolated-target source action',
        purpose='exact free-target self scattering; analytic bX retained'))
    j0 = q+v@free_q
    vq = np.asarray(v@Q)
    free_vq = ports.restrict(free_target_factor.solve(
        ports.inject(vq), label=label+' exact isolated-target basis action',
        purpose='exact free-target self scattering; C0 is not projected'))
    kq = vq+v@free_vq
    raw = dict(q=q, j0=j0, KQ=kq, VQ=vq, bX_ports=bx,
               free_target_q_ports=free_q, free_target_VQ_ports=free_vq)
    metadata = dict(target_self='exact same-space sparse FE factor L0-VX',
                    target_forcing='analytic bX cut integral plus V uw; no nodal incident substitution',
                    target_self_projection=False,
                    setup_rhs=int(Q.shape[1]+q.shape[1]),
                    source_count=int(q.shape[1]))
    return TargetSelfTemplate(q, j0, kq, vq, v, bx, raw, metadata)


@dataclass(frozen=True)
class ClosureSolution:
    current: np.ndarray
    core: np.ndarray
    coefficients: np.ndarray
    reduced_rhs: np.ndarray
    metadata: dict[str, Any]


def _small_closure(base, left, Q, *, rank, method):
    closure = np.eye(rank, dtype=complex)-Q.conj().T@left
    reduced_rhs = Q.conj().T@base
    coefficients = np.linalg.solve(closure, reduced_rhs)
    current = base+left@coefficients
    if not np.isfinite(current).all():
        raise np.linalg.LinAlgError('nonfinite approximate closure current')
    singular = np.linalg.svd(closure, compute_uv=False)
    metadata = dict(method=method, rank=int(rank),
                    reduced_closure='I - Q_r^* left; all source RHS solved together',
                    closure_singular_values=singular.tolist(),
                    closure_condition_2=ratio(singular[0], singular[-1]),
                    closure_inverse_norm_2=ratio(1., singular[-1]),
                    closure_norms='numerical values of this small coordinate matrix; no certified physical resolvent bound',
                    reduced_backward=linear_defect(closure, coefficients, reduced_rhs),
                    current_relation_relative=relative_error(
                        current, base+left@(Q.conj().T@current)),
                    algebra_only=True, full_feedback_validation='separate paid Gw current action required')
    return ClosureSolution(current, closure, coefficients, reduced_rhs, metadata)


def solve_return(approximation, template, rank):
    q, hcore, _ = approximation.rank(rank)
    return _small_closure(template.j0, template.KQ[:, :rank]@hcore, q,
                          rank=rank, method='return_only_exact_target_self')


def solve_whole_same_basis(approximation, template, rank):
    """Diagnostic only: approximate all Cw in the *same H-fitted basis*.

    This deliberately removes the exact isolated-target factor from the
    online closure. It is not the optimal whole-C low-rank method and does
    not establish a general superiority claim for return-only compression.
    """
    q, _, ccore = approximation.rank(rank)
    return _small_closure(template.q, template.VQ[:, :rank]@ccore, q,
                          rank=rank, method='whole_C_same_return_basis_diagnostic')


def quadratic_norm(value, mass):
    a = np.asarray(value)
    energy = float(np.sum((a.conj()*(mass@a)).real))
    # The matrix is a physical L2 Gram matrix.  Only a floating-point-sized
    # negative quadratic-form sum is rounded; a genuine negative fails.
    scale = float(np.linalg.norm(a)*np.linalg.norm(mass@a))
    if energy < -128*np.finfo(float).eps*max(scale, np.finfo(float).tiny):
        raise RuntimeError('ROI L2 mass produced a negative quadratic form')
    return float(np.sqrt(max(0., energy)))


def measure_field(ports, mass_ports, model_u, reference_u, wall_u):
    """L2 is a field diagnostic only; it never changes the QR metric."""
    delta, delta_reference = model_u-wall_u, reference_u-wall_u
    u_port, ref_port = ports.restrict(model_u), ports.restrict(reference_u)
    d_port, dref_port = ports.restrict(delta), ports.restrict(delta_reference)
    return dict(scattered_free_FEcoeff_relative=relative_error(model_u, reference_u),
                target_induced_free_FEcoeff_relative=relative_error(delta, delta_reference),
                ROI_scattered_FEcoeff_relative=relative_error(u_port, ref_port),
                ROI_target_induced_FEcoeff_relative=relative_error(d_port, dref_port),
                ROI_scattered_L2_absolute=quadratic_norm(u_port-ref_port, mass_ports),
                ROI_scattered_L2_relative=ratio(quadratic_norm(u_port-ref_port, mass_ports),
                                               quadratic_norm(ref_port, mass_ports)),
                ROI_target_induced_L2_absolute=quadratic_norm(d_port-dref_port, mass_ports),
                ROI_target_induced_L2_relative=ratio(quadratic_norm(d_port-dref_port, mass_ports),
                                                    quadratic_norm(dref_port, mass_ports)),
                L2_scope='P3 FE field on fixed public target ROI, summed over sources; not an operator spectrum')
