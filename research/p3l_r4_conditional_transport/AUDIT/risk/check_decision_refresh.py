"""Small arbitrary complex linear algebra tests; not Maxwell simulation."""
import json,time,resource,sys
from pathlib import Path
import numpy as np
from p3r4.decision_refresh import *
wall,cpu=time.monotonic(),time.process_time();rng=np.random.default_rng(20261025043)
n,sources,receivers=18,3,7
G=.008*(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))
S=.1*(rng.normal(size=(receivers,n))+1j*rng.normal(size=(receivers,n)))
e=rng.normal(size=(n,sources))+1j*rng.normal(size=(n,sources));chi0=.5+.1j+np.zeros(n)
L0=np.eye(n)-chi0[:,None]*G;j0=np.linalg.solve(L0,chi0[:,None]*e);E0=e+G@j0
charged=[0]
def acquire(rhs):charged[0]+=rhs.shape[1];return np.linalg.solve(L0,rhs),rhs.shape[1]
R0=charged_inverse(acquire,n)
# The new defect formula must remain valid for a genuinely approximate R0.
R0=R0+1e-5*(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))
cache=DecisionRefresh(chi0,j0,E0,G,S,e,R0,.01,n,True)
observation=S@j0+.01*(rng.normal(size=(receivers,sources))+1j*rng.normal(size=(receivers,sources)))
base_phi=float(np.linalg.norm(S@j0-observation)**2/.01**2)
misses=0;mismatches=0;resolved=0;unknown=0;max_inverse_ratio=0.;preconditioner_charged=[0]
def precondition(rhs):preconditioner_charged[0]+=1;return np.linalg.solve(L0,rhs),1
for candidate in range(100):
 chi=chi0+.002*(rng.normal(size=n)+1j*rng.normal(size=n));logu=float(np.log(rng.uniform()))
 L=np.eye(n)-chi[:,None]*G;truth=np.linalg.solve(L,chi[:,None]*e);phi=float(np.linalg.norm(S@truth-observation)**2/.01**2)
 constants=cache.bound_constants(chi)
 exactnorm=np.linalg.norm(np.linalg.solve(L,np.eye(n)),2)
 if constants['inverse_upper'] is not None:
  assert exactnorm<=constants['inverse_upper']*(1+1e-12)
  max_inverse_ratio=max(max_inverse_ratio,float(exactnorm/constants['inverse_upper']))
 rr=cache.richardson_decision(chi,observation,base_phi,logu)
 for h in rr['history']:
  assert h['phi_lower']<=phi+1e-10 and phi<=h['phi_upper']+1e-10
 gm=cache.ordinary_krylov(chi,observation,base_phi,logu,paid_preconditioner=precondition)
 for arm in (rr,gm):
  p=arm['point'];misses+=not(p['phi_lower']<=phi+1e-10 and phi<=p['phi_upper']+1e-10)
  if arm['decision']=='UNKNOWN':unknown+=1
  else:
   resolved+=1;mismatches+=(arm['decision']=='ACCEPT')!=(logu<min(0.,base_phi-phi))
 assert rr['refinements']<=4 and max(gm['iterations_per_source'])<=8
far=chi0+100;unqualified=cache.richardson_decision(far,observation,base_phi,-1.)
assert unqualified['point']['eta_upper']>=1 and unqualified['decision']=='UNKNOWN'
try:charged_inverse(lambda rhs:(R0,n-1),n)
except AssertionError:pass
else:raise AssertionError('undercharged setup accepted')
try:cache.ordinary_krylov(chi,observation,base_phi,logu,paid_preconditioner=lambda rhs:(rhs,0))
except AssertionError:pass
else:raise AssertionError('uncharged preconditioner accepted')
assert misses==0 and mismatches==0 and cache.left_inverse_defect_F>0
result=dict(status='PASS',scope='arbitrary18x18 complex algebra only; not a Maxwell state/pilot',proposals=100,
             inverse_defect_F=cache.left_inverse_defect_F,bound_misses=misses,resolved_decisions=resolved,unknown=unknown,decision_mismatches=mismatches,
             maximum_true_inverse_norm_over_bound=max_inverse_ratio,eta_ge1_fallback='PASS',setup_undercharge_rejected='PASS',preconditioner_undercharge_rejected='PASS',
             synthetic_setup_callback_columns=charged[0],synthetic_preconditioner_callback_columns=preconditioner_charged[0],
             physics_RHS=0,physics_states=0,training_updates=0,wall_seconds=time.monotonic()-wall,cpu_seconds=time.process_time()-cpu,
             peak_rss_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)*(1 if sys.platform=='darwin' else 1024))
Path('research/p3l_r4_conditional_transport/AUDIT/risk/DECISION_REFRESH_MATH_CHECK.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
