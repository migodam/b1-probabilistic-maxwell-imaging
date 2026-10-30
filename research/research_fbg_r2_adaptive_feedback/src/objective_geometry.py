"""Directional feedback curvature from exact first/second current tangents.

No finite difference, material Hessian, or full n-by-n projector is required.
These diagnostics do not rank test material truth or certify a basin.
"""
import numpy as np
from scipy.linalg import cho_factor,cho_solve


def directional_geometry(Q,R,S,G,B,T1,T2):
    """Raw psi_B'' along one fixed real material direction at fixed B.

Q must already be the transformed anchor V. T1,T2 are physical current
derivatives. T2=2 L^-1 diag(h)G T1. All source Gram solves are small.
"""
    V=np.asarray(Q);V1=(G+B)@T1;V2=(G+B)@T2
    H=V.conj().T@V;cf=cho_factor(H,lower=True,check_finite=False)
    solve=lambda a:cho_solve(cf,a,check_finite=False)
    H1=V1.conj().T@V+V.conj().T@V1
    H2=V2.conj().T@V+V.conj().T@V2+2*V1.conj().T@V1
    R1,R2=S@T1,S@T2;T=solve(R.conj().T)
    linear_positive=2*np.trace(R1@solve(R1.conj().T)).real
    physical_second=2*np.trace(R2@T).real
    weight_cross=-4*np.trace(R1@solve(H1@T)).real
    weight_second=(2*np.trace(R@solve(H1@solve(H1@T)))-np.trace(R@solve(H2@T))).real
    pinv=solve(V.conj().T)
    U=V1-V@(pinv@V1)
    # U* [V (V*V)^-1]=0, giving exact low-rank Frobenius norm.
    projector_norm_squared=2*np.trace(solve(U.conj().T@U)).real
    return {'psi':float(np.trace(R@T).real),
            'psi_second':float(linear_positive+physical_second+weight_cross+weight_second),
            'positive_weighted_J_term':float(linear_positive),
            'physical_residual_second':float(physical_second),
            'signed_Gram_cross':float(weight_cross),'signed_Gram_second':float(weight_second),
            'transformed_projector_jvp_fro':float(np.sqrt(max(projector_norm_squared,0.))),
            'transformed_field_jvp_fro':float(np.linalg.norm(V1)),
            'Gram_condition':float(np.linalg.cond(H)),
            'interpretation':'one exact directional Hessian; Gram terms are signed, not a basin certificate'}
