"""A17 common-core current Galerkin exchange, with an explicit full tail.

This adapter operates on the scalar-TM current equation L T=diag(h)Q.
It does not import vector-DDA polarizabilities or certify the full Jacobian.
"""
import numpy as np
from scipy.linalg import lu_factor,lu_solve,qr


def current_schur(L,core,candidate,rhs):
    """Return retained and enriched tangents and exact finite Schur exchange."""
    U,C=np.asarray(core),np.asarray(candidate)
    A=U.conj().T@L@U
    a=lu_factor(A,check_finite=False)
    R=lambda b:U@lu_solve(a,U.conj().T@b,check_finite=False)
    retained=R(rhs)
    corrected=C-R(L@C)
    gamma=C.conj().T@L@corrected
    innovation=C.conj().T@(rhs-L@retained)
    exchange=corrected@np.linalg.solve(gamma,innovation)
    return retained,retained+exchange,{'gamma':gamma,'exchange':exchange,
                                    'corrected_candidates':corrected,'innovation':innovation}


def krylov_spaces(G,Q,core_rank=None,candidate_rank=None):
    """Public current-space bank from fields and propagated fields; no truth."""
    U,_=qr(Q,mode='economic',check_finite=False)
    U=U[:,:core_rank or U.shape[1]]
    raw=np.concatenate((G@Q,G@(G@Q)),axis=1)
    raw=raw-U@(U.conj().T@raw)
    C,_=qr(raw,mode='economic',check_finite=False)
    return U,C[:,:candidate_rank or Q.shape[1]]


def projector_tail(Q,G,exact_tangent,approx_tangent,feedback=None):
    """Common exact anchor: tail effect on projector JVP, not a basin bound."""
    B=np.zeros_like(G) if feedback is None else np.asarray(feedback)
    V=Q
    pinv=np.linalg.solve(V.conj().T@V,V.conj().T)
    tail=(G+B)@(exact_tangent-approx_tangent)
    transverse=tail-V@(pinv@tail)
    N=transverse@pinv
    delta=N+N.conj().T
    return {'projector_tail_fro':float(np.linalg.norm(delta)),
            'upper_bound':float(2*np.linalg.norm(transverse)*np.linalg.norm(pinv,2)),
            'current_tail_fro':float(np.linalg.norm(exact_tangent-approx_tangent))}
