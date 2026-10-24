"""Output3 local material-target representatives; no backend import or job launch.

min mean|decoder(x)-mu|² subject to ORIGINAL Full nonlinear Phi<=threshold.
Coordinates are full66 probit; the original Gaussianized decoder is unchanged.
This optimization output is not a posterior mean/sample and introduces no prior.
A caller supplies the budgeted Target, durable event hook and hard cost guard.
"""
from dataclasses import dataclass
import numpy as np
from scipy.optimize import minimize
from .core import from_probit,to_probit,probit_derivative,BudgetStop

PHI_THRESHOLD=282.873792
GRADIENT_MODES=('source_adjoint8','receiver_J32')


class LocalBudgetStop(RuntimeError):pass


def material_objective_gradient(prior,x,mu):
 x=np.asarray(x,float);mu=np.asarray(mu,complex)
 if x.shape!=(66,) or not np.isfinite(x).all() or not np.isfinite(float(x@x)):
  raise ValueError('finite full66 probit coordinate/energy required')
 q=from_probit(x);chi,_=prior.decode(q)
 if mu.shape!=chi.shape or not np.isfinite(mu).all():raise ValueError('finite material target in frozen full-grid layout required')
 D=prior.derivative(q)*probit_derivative(x)[None,:]
 difference=chi-mu
 value=float(np.vdot(difference,difference).real/len(chi))
 gradient=2*np.real(D.conj().T@difference)/len(chi)
 return value,gradient,chi


@dataclass
class Evaluation:
 x:np.ndarray
 point:object
 phi:float
 gradient_phi:np.ndarray
 material_objective:float
 gradient_material:np.ndarray
 origin:str
 accepted:bool=False


class ExactOracle:
 """Caches at most the caller's bounded trial set; no rejected trial is final."""
 def __init__(self,target,anchor,observation,mu,x0,gradient_mode,require,used,event):
  if gradient_mode not in GRADIENT_MODES:raise ValueError('explicit source8 or receiver32 gradient mode required')
  self.target=target;self.anchor=anchor;self.observation=np.asarray(observation);self.mu=np.asarray(mu,complex)
  self.x0=np.asarray(x0,float);self.mode=gradient_mode;self.require=require;self.used=used;self.event=event;self.cache={}
  q=from_probit(self.x0)
  if np.max(abs(q-anchor.q))>1e-12:raise ValueError('exact online anchor/probit context mismatch; no inward clipping')

 def fee(self,point=None):
  if point is not None and (getattr(point,'J',None) is not None or (getattr(point,'gradient',None) is not None and getattr(point,'gradient_observation',None) is not None and np.array_equal(point.gradient_observation,self.observation))):return 0
  if point is not None and getattr(point,'T',None) is not None:return 0
  return 8 if self.mode=='source_adjoint8' else 32

 def _gradient(self,point,x):
  # Reuse already-paid receiver T at the anchor by matrix composition only.
  if getattr(point,'J',None) is None and getattr(point,'T',None) is not None:self.target.jacobian(point)
  elif self.mode=='receiver_J32' and getattr(point,'J',None) is None:self.target.jacobian(point)
  native=np.asarray(self.target.phi_gradient(point,self.observation),float)
  if native.shape!=(66,) or not np.isfinite(native).all():raise ValueError('finite native likelihood gradient required')
  return native*probit_derivative(x)

 def get(self,x,*,initial=False,reserve_states=0,reserve_rhs=0):
  x=np.asarray(x,float);key=x.tobytes()
  if key in self.cache:return self.cache[key]
  value,grad,chi=material_objective_gradient(self.target.prior,x,self.mu)
  before=dict(self.used())
  if initial:
   self.require(states=0,RHS=self.fee(self.anchor));point=self.anchor;origin='original_qualified_anchor'
  else:
   gradient_fee=self.fee()
   self.require(states=1+reserve_states,RHS=8+gradient_fee+reserve_rhs)
   point=self.target.evaluate(from_probit(x),purpose='material_representative_trial');origin='paid_trial'
  phi=float(self.target.phi(point,self.observation))
  if not np.isfinite(phi):raise ValueError('nonfinite Full likelihood; trial cannot be selected')
  # Full endpoint is already durable in Target before this gradient attempt.
  try:gphi=self._gradient(point,x)
  except Exception as exc:
   self.event(dict(event='GRADIENT_FAILED',endpoint=point.key,x=x,Phi=phi,material_objective=value,error=repr(exc),cost_before=before,cost_after=dict(self.used())))
   raise
  evaluation=Evaluation(x.copy(),point,phi,gphi,value,grad,origin)
  self.cache[key]=evaluation
  self.event(dict(event='EVALUATED',endpoint=point.key,x=x,Phi=phi,material_objective=value,
    gradient_mode=self.mode,gradient_phi_inf=float(np.max(abs(gphi))),cost_before=before,cost_after=dict(self.used()),accepted_iterate=False))
  return evaluation


def kkt_diagnostic(evaluation,threshold,kind):
 slack=threshold-evaluation.phi
 active=bool(abs(slack)<=1e-6*max(1.,threshold))
 g=evaluation.gradient_material if kind=='material_representative' else evaluation.gradient_phi
 c=evaluation.gradient_phi;den=float(c@c)
 multiplier=max(0.,-float(g@c)/den) if kind=='material_representative' and active and den>1e-30 else 0.
 stationarity=g+multiplier*c
 return dict(scope='NUMERICAL_FIRST_ORDER_DIAGNOSTIC_NOT_GLOBAL_OPTIMUM_CERTIFICATE',
  active=active,constraint_slack=slack,feasibility_violation=max(0.,-slack),
  multiplier=multiplier,stationarity_inf=float(np.max(abs(stationarity))),
  complementarity=abs(multiplier*slack),gradient_phi_inf=float(np.max(abs(c))),
  gradient_material_inf=float(np.max(abs(evaluation.gradient_material))))


def solve_representative(target,anchor,observation,mu,*,x0=None,gradient_mode='source_adjoint8',
                        require,used,event=lambda row:None,phi_threshold=PHI_THRESHOLD,
                        kind='material_representative',max_iterations=40,restoration_scales=(.999,.5)):
 """Budgeted SLSQP with strict feasible accepted-iterate selection.

    Cost guard must count actual backend reservations even on failures. No
    Full or gradient callback is invoked without preflight budget headroom.
    Two possible restoration states are reserved inside the same cap. Only
    their actually checked feasible endpoints can become final representatives.
    Ordinary likelihood optimization is an independent same-cap baseline and
    never uses mu in objective/proposals/selection.
 """
 if kind not in ('material_representative','ordinary_likelihood'):raise ValueError('unknown representative/control kind')
 if not np.isfinite(phi_threshold) or phi_threshold<=0:raise ValueError('original positive Phi threshold required')
 if x0 is None:x0=to_probit(anchor.q)  # boundary => error, never clip the prior
 oracle=ExactOracle(target,anchor,observation,mu,x0,gradient_mode,require,used,event)
 initial_phi=float(target.phi(anchor,observation))
 anchor_backward=float(getattr(anchor.state,'residual',0.))
 if not target.prior.legal(anchor.q) or not np.isfinite(initial_phi) or initial_phi>phi_threshold or not np.isfinite(anchor_backward) or anchor_backward>1e-10:
  raise ValueError('only an existing Full-qualified original online anchor may initialize Output3')
 initial=oracle.get(np.asarray(x0,float),initial=True);initial.accepted=True
 best=initial;accepted=[];iterations=0;optimizer=None;stop_reason=None
 objective_cache={}
 def material_terms(x):
  key=np.asarray(x).tobytes()
  if key not in objective_cache:objective_cache[key]=material_objective_gradient(target.prior,x,mu)
  return objective_cache[key]
 reserve_states=2 if kind=='material_representative' else 0
 reserve_rhs=reserve_states*(8+(8 if gradient_mode=='source_adjoint8' else 32))
 def get(x):return oracle.get(x,reserve_states=reserve_states,reserve_rhs=reserve_rhs)
 def objective(x):return material_terms(x)[0] if kind=='material_representative' else get(x).phi
 def jacobian(x):return material_terms(x)[1] if kind=='material_representative' else get(x).gradient_phi
 def accepted_callback(x):
  nonlocal best,iterations
  e=get(x);e.accepted=True;accepted.append(e);iterations+=1
  feasible=e.phi<=phi_threshold
  improve=(e.material_objective<best.material_objective-1e-14*max(1.,best.material_objective)) if kind=='material_representative' else e.phi<best.phi-1e-12
  if feasible and improve:best=e
  event(dict(event='ACCEPTED_ITERATE',endpoint=e.point.key,x=e.x,Phi=e.phi,material_objective=e.material_objective,
    feasible=feasible,selected_best=best is e,iteration=iterations))
 constraints=[]
 if kind=='material_representative':constraints=[dict(type='ineq',fun=lambda x:phi_threshold-get(x).phi,jac=lambda x:-get(x).gradient_phi)]
 try:
  result=minimize(objective,np.asarray(x0,float),jac=jacobian,method='SLSQP',constraints=constraints,
    callback=accepted_callback,options=dict(maxiter=int(max_iterations),ftol=1e-11,disp=False))
  optimizer=dict(success=bool(result.success),status=int(result.status),message=str(result.message),
    iterations=int(result.nit),nfev=int(result.nfev),njev=int(result.njev))
 except LocalBudgetStop as exc:
  stop_reason='LOCAL_HARD_BUDGET_STOP';optimizer=dict(success=False,message=str(exc),iterations=iterations)
 except BudgetStop:raise
 except Exception as exc:
  stop_reason='LOCAL_SOLVER_OR_PHYSICS_EXCEPTION';optimizer=dict(success=False,message=repr(exc),iterations=iterations)
 # Feasibility restoration is separately accepted only after Full evaluation.
 restored=False;restoration_attempts=0
 if kind=='material_representative' and accepted:
  candidate=min(accepted,key=lambda e:e.material_objective)
  if candidate.material_objective<best.material_objective and candidate.phi>phi_threshold:
   for scale in restoration_scales:
    try:
     x=np.asarray(x0)+float(scale)*(candidate.x-np.asarray(x0));restoration_attempts+=1
     e=oracle.get(x)  # reserve released; still within exactly the same hard cap
     if e.phi<=phi_threshold and e.material_objective<best.material_objective-1e-14*max(1.,best.material_objective):
      e.accepted=True;best=e;restored=True
      event(dict(event='ACCEPTED_FEASIBILITY_RESTORATION',scale=float(scale),endpoint=e.point.key,x=e.x,Phi=e.phi,material_objective=e.material_objective))
      break
     event(dict(event='REJECTED_FEASIBILITY_RESTORATION',scale=float(scale),endpoint=e.point.key,Phi=e.phi,material_objective=e.material_objective))
    except LocalBudgetStop:break
    except BudgetStop:raise
    except Exception as exc:event(dict(event='RESTORATION_EXCEPTION',error=repr(exc)));break
 if best.phi>phi_threshold:raise AssertionError('final feasible-anchor invariant')
 result=dict(status='CANDIDATE_REPRESENTATIVE' if kind=='material_representative' else 'ORDINARY_LIKELIHOOD_CONTROL',kind=kind,
  optimizer=optimizer,stop_reason=stop_reason,endpoint=best.point.key,x=best.x.copy(),q=best.point.q.copy(),
  Phi=best.phi,white_residual=float(np.sqrt(2*best.phi)),phi_threshold=phi_threshold,
  material_objective=best.material_objective,anchor_material_objective=initial.material_objective,
  material_objective_improvement=initial.material_objective-best.material_objective,
  anchor_Phi=initial.phi,data_feasible=True,fallback_to_original_anchor=best is initial,
  final_origin='feasibility_restoration' if restored else 'original_anchor' if best is initial else 'accepted_SLSQP_iterate',
  accepted_iterations=iterations,restoration_attempts=restoration_attempts,cached_trial_points=len(oracle.cache)-1,
  actual_cost=dict(used()),kkt=kkt_diagnostic(best,phi_threshold,kind),gradient_mode=gradient_mode,
  posterior_mean=False,posterior_sample=False,new_prior=False,global_minimum_claim=False)
 return result,best.point
