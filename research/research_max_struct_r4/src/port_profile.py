"""Fixed-eta full-current profile, lambda0 control and lambda1 port augmentation.

Parent-derived envelope gradient; no novelty or convergence claim. All terms
use the same total observed-data normalization and unchanged material prior.
This module reads no files, truth, labels, or scoring artifacts.
"""
import numpy as np
from scipy import linalg as la


def receiver_resolvent(model,state):
    """C=S A^-1 from the retained Maxwell LU; all receiver RHS counted."""
    with model._physics.operation('PORT_PROFILE_RECEIVER_RESOLVENT',rhs=model.nr,
                                  receiver_transpose_RHS=model.nr,receiver_resolvent_builds=1):
        return la.lu_solve(state['lu'],model.S.T,trans=1,check_finite=False).T


def fixed_eta_reference(model,host=.1):
    m=np.full(model.n,float(host))
    state=model.state(m);C=receiver_resolvent(model,state)
    with model._physics.operation('PORT_PROFILE_REFERENCE_GRAM',reference_Gram_builds=1):
        H=C@C.conj().T
        eta=float(np.trace(H).real/model.nr)
    if not np.isfinite(eta)or eta<=0:raise ValueError('Positive finite public-reference eta required')
    return eta,m,C


class PortProfile:
    def __init__(self,model,Y,eta,lam,reference=None):
        self.model=model;self.Y=np.asarray(Y,complex)
        self.eta,self.lam=float(eta),float(lam)
        if self.eta<=0 or self.lam<0:raise ValueError('Fixed eta>0 and lambda>=0 required')
        self.scale=float(np.vdot(self.Y,self.Y).real)
        self.reference=reference
        self.cached_state=None;self.cached_C=None

    def evaluate(self,m,gradient=True):
        model=self.model;state=model.state(m)
        if state is self.cached_state:
            C=self.cached_C
        elif self.reference is not None and np.array_equal(state['m'],self.reference[0]):
            C=self.reference[1]
            model.cost['shared_reference_C_reuses']+=1
            self.cached_state,self.cached_C=state,C
        else:
            C=receiver_resolvent(model,state)
            self.cached_state,self.cached_C=state,C
        e=state['Y']-self.Y
        with model._physics.operation('PORT_PROFILE_SMALL_SYSTEM',profile_Gram_builds=1,
                                      profile_Cholesky=1,profile_small_RHS=model.ns,
                                      profile_system_dimension=model.nr):
            H=C@C.conj().T
            K=(1+self.lam)*H+self.eta*np.eye(model.nr)
            factor=la.cho_factor(K,lower=True,check_finite=False)
            z=la.cho_solve(factor,e,check_finite=False)
            q=-C.conj().T@z
            port=C@q
            fitted=e+port
            regularizer,greg,_=model.reg.evaluate(m)
            data=.5*float(np.vdot(fitted,fitted).real)/self.scale
            port_penalty=.5*self.lam*float(np.vdot(port,port).real)/self.scale
            state_penalty=.5*self.eta*float(np.vdot(q,q).real)/self.scale
            value=data+port_penalty+state_penalty+regularizer
            physical=.5*float(np.vdot(e,e).real)/self.scale+regularizer
        g=None
        if gradient:
            with model._physics.operation('PORT_PROFILE_CURRENT_INCREMENT',rhs=model.ns,
                                          profile_current_RHS=model.ns):
                increment=la.lu_solve(state['lu'],q,check_finite=False)
                jhat=state['j']+increment
                Qhat=model.E+model.G@jhat
            # Shared evaluation path for both fixed arms; lambda0 term is zero
            # but the same adjoint implementation is evaluated and charged.
            with model._physics.operation('PORT_PROFILE_ENVELOPE_ADJOINT',rhs=model.ns,
                                          adjoint_RHS=model.ns,profile_adjoint_RHS=model.ns):
                adj=la.lu_solve(state['lu'],model.S.conj().T@port,trans=2,check_finite=False)
                g=(-self.lam*np.real(np.sum(adj.conj()*state['Q'],axis=1))
                   -self.eta*np.real(np.sum(q.conj()*Qhat,axis=1)))/self.scale+greg
        if not np.isfinite(value)or (g is not None and not np.isfinite(g).all()):
            raise FloatingPointError('Nonfinite port profile value/gradient')
        detail=dict(auxiliary_objective=float(value),auxiliary_data_term=data,
                    auxiliary_port_term=port_penalty,auxiliary_state_term=state_penalty,
                    prior=float(regularizer),physical_objective=float(physical),
                    physical_relative_residual=float(la.norm(e)/np.sqrt(self.scale)),
                    fitted_current_data_relative_residual=float(la.norm(fitted)/np.sqrt(self.scale)),
                    port_change_relative=float(la.norm(port)/np.sqrt(self.scale)),
                    state_defect_fro=float(la.norm(q)),eta=self.eta,lambda_value=self.lam)
        return float(value),g,detail,state
