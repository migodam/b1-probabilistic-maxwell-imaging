"""Counted actual two-pulse-cell TM-VIE controls; never imported by the engine.

Uses the inherited Z-ALG Model (physical Hankel Green kernel). Every setup,
forward and adjoint operation enters the root atomic ledger before execution.
These tiny discrete examples are controls, not evidence of a 576-cell theorem.
"""
from contextlib import contextmanager
from pathlib import Path
import importlib.util
import json
import sys
import numpy as np
from scipy.linalg import null_space, cho_factor, cho_solve

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))
from research.p3l_r12_bilateral_fixedpoint.shared.budget import Budget

OUT = Path(__file__).resolve().parent/'two_cell'
SOURCE = REPO/'research/research_z_alg_r3/src/physics.py'
spec = importlib.util.spec_from_file_location('r12_tiny_tm_vie', SOURCE)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
SIGMA = .00119211798514519


class Book:
    def __init__(self):
        self.budget = Budget()
        self.rhs = 0
        self.lu = 0
        self.tickets = []

    @contextmanager
    def op(self, kind, rhs=0, **metadata):
        lu = int(metadata.pop('LU', 0))
        if self.rhs+rhs > 1000: raise RuntimeError('Tiny-control hard cap')
        self.rhs += rhs; self.lu += lu
        with self.budget.operation('R12_AGENT_B_TWO_CELL', LU=lu, RHS=rhs,
                                   order=2, kind=kind, **metadata) as ticket:
            self.tickets.append(ticket)
            yield ticket


def inverse_spd(A, b):
    return cho_solve(cho_factor(A, lower=True, check_finite=False), b, check_finite=False)


def run():
    OUT.mkdir(parents=True, exist_ok=False)
    book = Book()
    model = module.Model(book=book,
        centers=np.array([[-.12, 0.], [.12, 0.]]),
        sources=np.array([[0., -1.], [0., 1.]]),
        receivers=np.array([[0., -1.2], [0., 1.3], [0., 1.8]]), radius=.04)
    m0 = np.array([.65, 1.1])
    root = model.forward(m0, label='ROOT', jac=True)
    y0 = root['y'].copy()
    J0 = root['J']*np.sqrt(2)/SIGMA
    residual0 = model.pack(root['y']-y0)*np.sqrt(2)/SIGMA
    s = np.linalg.svd(J0, compute_uv=False)
    mu = 1.
    A0 = J0.T@J0+mu*np.eye(2)
    analytic_root = np.eye(2)-inverse_spd(A0, J0.T@J0)
    mapped = []
    h = 1e-5
    Jplus, Jminus = [], []
    def get(x, label):
        q = model.forward(x, label=label, jac=True)
        return q['y'].copy(), q['J'].copy()*np.sqrt(2)/SIGMA
    # These paid J states support both FD DT and DJ; none are free probes.
    for i in range(2):
        ei = np.eye(2)[i]
        yp, jp = get(m0+h*ei, f'FD_PLUS_{i}')
        ym, jm = get(m0-h*ei, f'FD_MINUS_{i}')
        Jplus.append(jp); Jminus.append(jm)
        rp = model.pack(yp-y0)*np.sqrt(2)/SIGMA
        rm = model.pack(ym-y0)*np.sqrt(2)/SIGMA
        tp = m0+h*ei-inverse_spd(jp.T@jp+mu*np.eye(2), jp.T@rp)
        tm = m0-h*ei-inverse_spd(jm.T@jm+mu*np.eye(2), jm.T@rm)
        mapped.append((tp-tm)/(2*h))
    fd_root = np.column_stack(mapped)
    # Symmetry produces two material-distinct responses. Full-column-rank real
    # whitened J at each root certifies local isolation, never global uniqueness.
    swapped = model.forward(m0[::-1], label='PERMUTED_ROOT', jac=True)
    swapped_J = swapped['J']*np.sqrt(2)/SIGMA
    pair_defect = float(np.linalg.norm(model.pack(swapped['y']-y0)*np.sqrt(2)/SIGMA))
    # Construct an observed residual orthogonal to columns of actual Maxwell J.
    # It is a declared algebraic data perturbation in this tiny physical control.
    r = null_space(J0.T)[:, 0]
    r /= np.linalg.norm(r)
    y_nonroot = y0-model.unpack(r*SIGMA/np.sqrt(2))
    g = J0.T@r
    step = inverse_spd(A0, g)
    exact_nonroot, omitted_nonroot = [], []
    nonroot_fd = []
    for i in range(2):
        ei = np.eye(2)[i]
        DJ = (Jplus[i]-Jminus[i])/(2*h)
        DA = DJ.T@J0+J0.T@DJ
        exact_nonroot.append(ei-inverse_spd(A0, J0.T@J0@ei+DJ.T@r)+inverse_spd(A0, DA@step))
        omitted_nonroot.append(ei-inverse_spd(A0, J0.T@J0@ei)+inverse_spd(A0, DA@step))
        # Reusing saved FD y/J is algebra, with no additional solver query.
        # The real residual at each FD node is recovered from saved root maps.
        yp, jp = get(m0+h*ei, f'NONROOT_FD_PLUS_{i}')
        ym, jm = get(m0-h*ei, f'NONROOT_FD_MINUS_{i}')
        rp = model.pack(yp-y_nonroot)*np.sqrt(2)/SIGMA
        rm = model.pack(ym-y_nonroot)*np.sqrt(2)/SIGMA
        tp = m0+h*ei-inverse_spd(jp.T@jp+mu*np.eye(2), jp.T@rp)
        tm = m0-h*ei-inverse_spd(jm.T@jm+mu*np.eye(2), jm.T@rm)
        nonroot_fd.append((tp-tm)/(2*h))
    exact_nonroot = np.column_stack(exact_nonroot)
    omitted_nonroot = np.column_stack(omitted_nonroot)
    nonroot_fd = np.column_stack(nonroot_fd)
    results = dict(status='PASS', source=str(SOURCE), source_physics='actual pulse TM VIE Hankel Green',
        roots=dict(materials=[m0.tolist(), m0[::-1].tolist()],
                   whitened_pair_defect=pair_defect,
                   material_RMS_distance=float(np.linalg.norm(m0-m0[::-1])/np.sqrt(2)),
                   whitened_J_singular_values=s.tolist(),
                   second_J_singular_values=np.linalg.svd(swapped_J, compute_uv=False).tolist(),
                   real_J_rank=int(np.linalg.matrix_rank(J0)),
                   local_isolation_condition=bool(np.linalg.matrix_rank(J0)==2),
                   scope='two isolated roots in this real two-cell parameterization; no continuous fiber proof'),
        root_map=dict(mu=mu, eigenvalues=np.linalg.eigvalsh(analytic_root).tolist(),
                      predicted_mu_over_sigma2_plus_mu=(mu/(s*s+mu)).tolist(),
                      derivative_FD_absolute_error=float(np.linalg.norm(fd_root-analytic_root)),
                      exact_unit_modes=0 if np.linalg.matrix_rank(J0)==2 else None),
        nonroot=dict(residual_norm=float(np.linalg.norm(r)), gradient_norm=float(np.linalg.norm(g)),
                     full_DT_FD_error=float(np.linalg.norm(nonroot_fd-exact_nonroot)),
                     omitted_residual_term_error=float(np.linalg.norm(nonroot_fd-omitted_nonroot)),
                     residual_Hessian_term_norm=float(np.linalg.norm(exact_nonroot-omitted_nonroot))),
        cost=dict(conservative_RHS=book.rhs, physical_LU=book.lu,
                  model_counts=dict(model.counts), tickets=book.tickets),
        truth_read=False, posterior_claim=False, continuous_PDE_claim=False)
    if pair_defect > 1e-8 or np.linalg.matrix_rank(J0) != 2:
        results['status'] = 'CONTROL_CONDITION_NOT_MET'
    np.savez_compressed(OUT/'ARRAYS.npz', J=J0, root_DT=analytic_root,
                        root_DT_FD=fd_root, nonroot_DT=exact_nonroot,
                        nonroot_DT_FD=nonroot_fd, omitted_DT=omitted_nonroot,
                        residual=r, materials=np.array([m0,m0[::-1]]), y=y0)
    (OUT/'RESULTS.json').write_text(json.dumps(results, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(status=results['status'], RHS=book.rhs, LU=book.lu,
                         pair_defect=pair_defect, J_singular=s.tolist(),
                         root_FD_error=results['root_map']['derivative_FD_absolute_error'],
                         nonroot_FD_error=results['nonroot']['full_DT_FD_error'])))


if __name__ == '__main__': run()
