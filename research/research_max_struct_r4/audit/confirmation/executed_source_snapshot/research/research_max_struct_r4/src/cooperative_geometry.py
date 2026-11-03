"""Exact receiver-visible curvature using the existing retained Maxwell LU.

This is a dense contraction of the established R4-BASIN Hessian, not a new
forward solver or a new Hessian theorem. The exact TV Hessian differs from IRLS.
"""
import numpy as np
from scipy import linalg as la


def exact_regularizer(model, m):
    reg = model.reg
    d = reg.D @ m
    q = np.sqrt(d*d + reg.epsilon**2)
    weights = (2*reg.smooth + reg.tv*reg.epsilon**2/q**3)/reg.edges
    H = (reg.D.T @ reg.D.multiply(weights[:,None])).toarray()
    value, gradient, _ = reg.evaluate(m)
    return value, gradient, H


def curvature(model, m, Y, indices):
    state = model.state(m)
    ids = np.asarray(indices, int)
    Q = state['Q'][ids]
    residual = state['Y']-Y
    scale = float(np.vdot(Y,Y).real)
    with model._physics.operation('COOPERATIVE_CURVATURE_COLUMNS', rhs=len(ids),
                                  resolvent_RHS=len(ids), curvature_columns=len(ids)):
        columns = la.lu_solve(state['lu'], np.eye(model.n)[:,ids], check_finite=False)
        O = model.S @ columns
        C = model.G[ids] @ columns
    with model._physics.operation('COOPERATIVE_CURVATURE_CONTRACTIONS',
                                  curvature_matrix_dimension=len(ids), dense_curvature_builds=1):
        GN = ((O.conj().T@O)*(Q.conj()@Q.T)).real/scale
        T = (O.T@residual.conj()@Q.T)*C
        feedback = (T+T.T).real/scale
        _, _, regH = exact_regularizer(model, m)
        regularizer = regH[np.ix_(ids,ids)]
        H = GN+feedback+regularizer
    return dict(state=state, indices=ids, O=O, C=C, Q=Q,
                GN=GN, feedback=feedback, regularizer=regularizer, H=H)


def inherited_hvp(model, state, Y, direction):
    """Same analytic derivative as R4-BASIN hessian_geometry.py, adapted API."""
    Q, factor = state['Q'], state['lu']
    with model._physics.operation('HVP_WIRING_BASE_ADJOINT', rhs=model.ns, adjoint_RHS=model.ns):
        adj = la.lu_solve(factor, model.S.conj().T@(state['Y']-Y), trans=2, check_finite=False)
    with model._physics.operation('HVP_WIRING_TANGENT', rhs=2*model.ns,
                                  tangent_RHS=model.ns, incremental_adjoint_RHS=model.ns):
        U = la.lu_solve(factor, direction[:,None]*Q, check_finite=False)
        dQ = model.G@U
        rhs = model.G.conj().T@(direction[:,None]*adj)+model.S.conj().T@(model.S@U)
        da = la.lu_solve(factor, rhs, trans=2, check_finite=False)
        out = np.sum(da.conj()*Q+adj.conj()*dQ,axis=1).real/np.vdot(Y,Y).real
    return out+exact_regularizer(model,state['m'])[2]@direction
