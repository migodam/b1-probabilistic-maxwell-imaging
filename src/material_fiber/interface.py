"""Small reusable API; every full solve still requires an active charged Target."""
import numpy as np
from .explore import geometry,direction,constrained_fit
from p1_material.proposals import frozen_latent_step

def evaluate(target,q,purpose='external_evaluation'):
    return target.evaluate(q,purpose)

def propose(target,parent,anchor,radius,rng,method='soft_adaptive'):
    target.parent=parent;g=geometry(target,anchor,radius);v,w=direction(target,anchor,g,method,rng)
    if v is None:return None,dict(status='NO_FEASIBLE_DIRECTION')
    q,matched=frozen_latent_step(target.prior,anchor.q,v,radius)
    return q,dict(radius_matched=matched,material_velocity=w,geometry=g,semantics='candidate_only')

def correct(target,parent,measured_point,origin,velocity,radius,seed):
    target.parent=parent;target.jacobian(measured_point)
    return constrained_fit(target,measured_point,origin,velocity,radius,seed)

def transition(target,current,observation,rng,temperature=1.,beta=.02,property_sd=.01):
    """Prior-reversible pCN/reflection: full-likelihood MH, without correction."""
    q=target.prior.local_prior_reversible(current.q,rng,beta,property_sd)
    proposed=target.evaluate(q,'standalone_prior_reversible_MH')
    ratio=-temperature*(target.phi(proposed,observation)-target.phi(current,observation))
    accepted=bool(np.log(rng.uniform())<min(0.,ratio))
    return (proposed if accepted else current),dict(accepted=accepted,log_MH_ratio=float(ratio),temperature=temperature,posterior_convergence='NOT_ESTABLISHED_BY_ONE_TRANSITION')

def archive(target,point,parent,method,origin=None,**metadata):
    return target.record(point,parent,method,origin,**metadata)
