"""Complex8-RHS global GMRES math only; zero Maxwell queries or paid runs."""
from pathlib import Path
import json,time,resource,sys
import numpy as np
from scipy.sparse.linalg import LinearOperator,gmres
from p3r4.decision_refresh_global_gmres import DecisionRefreshGlobalGMRES,resolve_global
from p3r4.interval_chain import IntervalState,approximate_state
wall,cpu=time.monotonic(),time.process_time();rng=np.random.default_rng(20261025162)
n,sources,receivers=14,8,7;G=.025*(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)));S=.1*(rng.normal(size=(receivers,n))+1j*rng.normal(size=(receivers,n)))
e=rng.normal(size=(n,sources))+1j*rng.normal(size=(n,sources));chi0=np.full(n,.4+.1j);L0=np.eye(n)-chi0[:,None]*G;R0=np.linalg.solve(L0,np.eye(n));j0=R0@(chi0[:,None]*e);E0=e+G@j0;sigma=.03
y=S@j0+.03*(rng.normal(size=(receivers,sources))+1j*rng.normal(size=(receivers,sources)))
small=[0]
def make():
 cache=DecisionRefreshGlobalGMRES(chi0,j0,E0,G,S,e,R0,sigma,n,True)
 def solve(H,b):
  assert H.shape[0]<=5 and H.shape[1]<=4 and b.shape==(H.shape[0],1)
  small[0]+=1;return np.linalg.lstsq(H,b,rcond=None)[0],1
 cache.configure_reduced_solver(solve);return cache
cache=make();chi=chi0+.05*(rng.normal(size=n)+1j*rng.normal(size=n));x=rng.normal(size=66)
candidate=approximate_state(cache,x,x,chi,y);initial=candidate.point['current'].copy();dummy=IntervalState(np.zeros(66),np.zeros(66),chi0,{},True)
cache.retain_endpoint_pair(dummy,candidate);errors=[]
def K(v):
 a=v.reshape(n,sources);return (R0@(a-chi[:,None]*(G@a))).ravel()
op=LinearOperator((n*sources,n*sources),matvec=K,dtype=complex);rhs=(R0@(chi[:,None]*e)).ravel()
for iteration in range(1,5):
 candidate.point=cache.refine_state(candidate,y);candidate.refinements+=1
 reference,info=gmres(op,rhs,x0=initial.ravel(),rtol=0.,atol=0.,restart=iteration,maxiter=1)
 error=float(np.linalg.norm(candidate.point['current'].ravel()-reference));assert error<1e-10;errors.append(error)
assert cache.counts['global_gmres_reduced_LS_complex_RHS']==4
# Sameworld cachedsolution may be reused, but no uncharged new least-squares.
copy=approximate_state(cache,x,x,chi,y);before=small[0];copy.point=cache.refine_state(copy,y);assert small[0]==before
# Force wrong reducedcharge to be rejected before any claim of a free solve.
bad=make();bad.configure_reduced_solver(lambda H,b:(np.linalg.lstsq(H,b,rcond=None)[0],0));s=approximate_state(bad,x,x,chi,y)
try:bad.refine_state(s,y)
except AssertionError:pass
else:raise AssertionError('unchargedLS accepted')
# Independent mock exact values for both endpoints, including accuracy-first.
def exact(state):
 current=np.linalg.solve(np.eye(n)-state.chi[:,None]*G,state.chi[:,None]*e);field=e+G@current;pred=S@current;phi=float(np.linalg.norm(pred-y)**2/sigma**2)
 return IntervalState(state.x.copy(),state.q.copy(),state.chi.copy(),dict(current=current,field=field,prediction=pred,residual=state.chi[:,None]*field-current,phi_lower=phi,phi_upper=phi,approximate_phi=phi,inverse_upper=None),True)
misses=0;accuracy_violations=0;fallbacks=0;refined=0;cache=make();state=IntervalState(np.zeros(66),np.zeros(66),chi0,{});state=exact(state)
for step in range(40):
 proposed=np.sqrt(1-.002**2)*state.x+.002*rng.normal(size=66);material=chi0+.5*(proposed[:n]+1j*proposed[n:2*n]);candidate=approximate_state(cache,proposed,proposed,material,y);cur_phi=exact(state).lower;new_phi=exact(candidate).lower;logu=float(np.log(rng.uniform()))
 def paid(s,role):return exact(s),sources
 answer=resolve_global(cache,state,candidate,y,logu,paid,fixed_accuracy=True);refined+=answer['candidate'].refinements;fallbacks+=answer['Full_calls']
 for row in answer['history']:
  assert row['current_lower']<=cur_phi+1e-10<=row['current_upper']+2e-10 and row['candidate_lower']<=new_phi+1e-10<=row['candidate_upper']+2e-10
 misses+=answer['accepted']!=(logu<min(0.,cur_phi-new_phi));accuracy_violations+=not(answer['history'][-1]['current_accuracy'] and answer['history'][-1]['candidate_accuracy'])
 state=answer['candidate'] if answer['accepted'] else answer['current'];assert len(cache._basis_cache)<=2
assert misses==0 and accuracy_violations==0
result=dict(status='PASS',scope='arbitrary14x14complexsystem,8sources; no Maxwell/GT/paidprocess',global_steps=4,global_iteration_current_errors_vs_independent_SciPy=errors,
 toy_chain_steps=40,decision_mismatches=misses,fixed_accuracy_violations=accuracy_violations,candidate_refinements=refined,mockFullfallbacks=fallbacks,
 bounded_two_endpoint_basis=True,source_columns=8,synthetic_reduced_complex_RHS=small[0],zerocharge_rejected=True,physics_RHS=0,physics_queries=0,training_updates=0,
 wall_seconds=time.monotonic()-wall,cpu_seconds=time.process_time()-cpu,peak_rss_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)*(1 if sys.platform=='darwin' else 1024))
Path('research/p3l_r4_conditional_transport/AUDIT/global_gmres/GLOBAL_GMRES_MATH_CHECK.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
