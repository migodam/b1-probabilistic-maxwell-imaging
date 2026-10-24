"""Separate Full-J Gauss-Newton trust-region ordinary control; no mu/GT API."""
import numpy as np
from .core import from_probit,probit_derivative

def exact_psd_trust_region(H,g,radius):
 """Eigen/secular solution of min g's+.5s'Hs, ||s||<=radius, H PSD.

 Only real parameter algebra, no complex-current solve. Any clipped negative
 eigenvalues must be within a dimension-scaled floating Gram roundoff bound.
 """
 H=np.asarray(H,float);g=np.asarray(g,float);H=(H+H.T)/2
 if H.shape!=(len(g),len(g)) or not radius>0 or not np.isfinite(H).all() or not np.isfinite(g).all():raise ValueError('finitePSDsubproblemrequired')
 w,Q=np.linalg.eigh(H);floor=100*np.finfo(float).eps*len(g)*max(1.,float(np.max(abs(w))))
 if w.min() < -floor:raise ValueError('nonPSDGNHessian')
 clipped=int(np.count_nonzero(w<0));w=np.maximum(w,0);c=Q.T@g
 if np.linalg.norm(g)==0:return np.zeros_like(g),dict(multiplier=0.,boundary=False,step_norm=0.,KKT_residual=0.,negative_eigen_roundoff_clips=clipped)
 s0=np.zeros_like(c);positive=w>0;s0[positive]=-c[positive]/w[positive]
 interior=not np.any(c[~positive]!=0) and np.linalg.norm(s0)<=radius
 if interior:multiplier=0.;s=Q@s0
 else:
  lo=0.;hi=max(1.,float(np.linalg.norm(g)/radius),float(w.max()))
  norm=lambda v:float(np.linalg.norm(c/(w+v)))
  while norm(hi)>radius:hi*=2
  for _ in range(240):
   mid=(lo+hi)/2
   if norm(mid)>radius:lo=mid
   else:hi=mid
   if hi==lo or abs(norm(hi)-radius)<=1e-13*max(radius,1e-12):break
  multiplier=hi;s=-Q@(c/(w+multiplier))
 residual=float(np.linalg.norm(g+H@s+multiplier*s));pred=-float(g@s+.5*s@H@s)
 return s,dict(multiplier=multiplier,boundary=not interior,step_norm=float(np.linalg.norm(s)),radius=radius,
  KKT_residual=residual,predicted_reduction=pred,negative_eigen_roundoff_clips=clipped,H_eigen_min=float(w.min()),H_eigen_max=float(w.max()))

def stationarity(r,J,g,H,x,objective,phi):
 w,Q=np.linalg.eigh((H+H.T)/2);cut=np.finfo(float).eps*len(g)*max(1.,float(np.max(abs(w))));active=w>cut;c=Q.T@g
 decrement2=float(np.sum(c[active]**2/w[active]))
 return dict(scope='exactFullJgradient;GNcurvaturedecrement,notfullnonlinearHessian/globalcertificate',gradient_total_inf=float(np.max(abs(g))),gradient_data_inf=float(np.max(abs(J.T@r))),
  scaled_GN_gradient_inf=float(np.max(abs(g)/np.sqrt(np.maximum(np.diag(H),1)))),GN_Newton_decrement_squared=decrement2,GN_Newton_decrement=float(np.sqrt(max(0.,decrement2))),
  GN_pseudoinverse_rank=int(active.sum()),GN_pseudoinverse_eigen_cutoff=cut,gradient_in_discarded_subspace_norm=float(np.linalg.norm(c[~active])),objective=objective,Phi=phi,prior_energy=float(.5*x@x))

def solve_control(target,anchor,x0,y,*,regularization,require,used,event,save_derivative):
 """Same16state/256RHS cap supplied externally; accepted step gets fresh J32.

 Preflight40RHS before each trial guarantees an accepted endpoint's J and KKT
 can always be bought. Rejected endpoints cost only Full8. Cached initialJ
 is supplied by the runner and separately disclosed as previously paid.
 """
 if regularization not in (0,1):raise ValueError('two preregistered objectives only')
 x=np.asarray(x0,float).copy();point=anchor;radius=.5;accepted=0;rejected=0;iterations=0;stop='HARD_LOCAL_BUDGET';history=[]
 if anchor.J is None:raise ValueError('paid cached initial nativeJrequired')
 def terms(p,v):
  r=target.white(p.state.data-y);J=p.J*probit_derivative(v)[None,:]
  phi=float(.5*r@r);objective=phi+.5*regularization*float(v@v);g=J.T@r+regularization*v;H=J.T@J+regularization*np.eye(len(v))
  return r,J,phi,objective,g,H
 r,J,phi,objective,g,H=terms(point,x);initial_kkt=stationarity(r,J,g,H,x,objective,phi);save_derivative(point,x,J,g,regularization,'initial_cached')
 for iterations in range(1,65):
  step,sub=exact_psd_trust_region(H,g,radius);pred=-float(g@step+.5*step@H@step)
  if pred<=np.finfo(float).eps*max(1.,abs(objective)) or np.linalg.norm(step)<=1e-12:stop='NO_NUMERICAL_MODEL_DESCENT';break
  try:require(states=1,RHS=40)
  except Exception as exc:
   if exc.__class__.__name__ not in ('LocalBudgetStop','ControlBudgetStop'):raise
   stop='HARD_LOCAL_BUDGET_BEFORE_FULL_PLUS_ACCEPTED_J';break
  candidate_x=x+step;candidate=target.evaluate(from_probit(candidate_x),purpose=f'ordinary_fullJ_TR_lambda{regularization}_trial{iterations}')
  trial_phi=target.phi(candidate,y);trial_objective=trial_phi+.5*regularization*float(candidate_x@candidate_x)
  actual=objective-trial_objective;rho=float(actual/pred);accept=bool(np.isfinite(rho) and rho>.1 and actual>0)
  oldradius=radius
  if not np.isfinite(rho) or rho<.25:radius*=.25
  elif rho>.75 and np.linalg.norm(step)>=.99*oldradius:radius=min(2.,2*oldradius)
  record=dict(iteration=iterations,lambda_=regularization,current_endpoint=point.key,candidate_endpoint=candidate.key,current_Phi=phi,trial_Phi=trial_phi,
   current_objective=objective,trial_objective=trial_objective,predicted_reduction=pred,actual_reduction=float(actual),ratio=rho,accepted=accept,
   radius_before=oldradius,radius_after=radius,step_norm=float(np.linalg.norm(step)),subproblem=sub,cost=dict(used()))
  if accept:
   require(states=0,RHS=32);target.jacobian(candidate)
   x=candidate_x;point=candidate;r,J,phi,objective,g,H=terms(point,x);accepted+=1
   save_derivative(point,x,J,g,regularization,'fresh_accepted_J32');record['cost']=dict(used())
  else:rejected+=1
  history.append(record);event(record)
  if radius<1e-12:stop='NUMERICAL_RADIUS_FLOOR';break
 else:stop='FIXED64NUMERICAL_ITERATION_CAP'
 kkt=stationarity(r,J,g,H,x,objective,phi)
 return dict(status='ORDINARY_FULLJ_TR_CONTROL',lambda_=regularization,x=x,q=point.q,endpoint=point.key,Phi=phi,objective=objective,kkt=kkt,initial_kkt=initial_kkt,
  accepted_updates=accepted,rejected_trials=rejected,completed_trials=len(history),stop_reason=stop,radius_final=radius,actual_cost=dict(used()),
  original_anchor_fallback=accepted==0,GT_read=False,mu_used=False,posterior_reference=False,global_optimum_claim=False),point
