"""Independent toy Full/interval chains; zero Maxwell work."""
from pathlib import Path
import json,time,resource,sys
import numpy as np
from p3r4.interval_chain import *
from p3r4.decision_refresh import DecisionRefresh
wall,cpu=time.monotonic(),time.process_time();rng=np.random.default_rng(20261025044)
n,sources,m=12,3,5;G=.008*(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)));S=.1*(rng.normal(size=(m,n))+1j*rng.normal(size=(m,n)))
e=rng.normal(size=(n,sources))+1j*rng.normal(size=(n,sources));chi0=np.full(n,.4+.1j);L0=np.eye(n)-chi0[:,None]*G;R0=np.linalg.solve(L0,np.eye(n));j0=R0@(chi0[:,None]*e);E0=e+G@j0;sig=.01
y=S@j0+.01*(rng.normal(size=(m,sources))+1j*rng.normal(size=(m,sources)))
def material(x):return chi0+.1*(x[:n]+1j*x[n:2*n])
def exact(x,q,chi):
 j=np.linalg.solve(np.eye(n)-chi[:,None]*G,chi[:,None]*e);E=e+G@j;pred=S@j;phi=float(np.linalg.norm(pred-y)**2/sig**2)
 return IntervalState(x.copy(),q.copy(),chi.copy(),dict(current=j,field=E,prediction=pred,residual=chi[:,None]*E-j,phi_lower=phi,phi_upper=phi,approximate_phi=phi,inverse_upper=None),True)
def run_interval():
 cache=DecisionRefresh(chi0,j0,E0,G,S,e,R0,sig,n,True);stream=np.random.default_rng(777);x=np.zeros(66);current=exact(x,x,chi0);paths=[x.copy()];flags=[];intervals=[];fullcalls=[0];approxaccepted=0
 def paid(state,role):fullcalls[0]+=1;return exact(state.x,state.q,state.chi),sources
 for step in range(512):
  x=np.sqrt(1-.002**2)*current.x+.002*stream.normal(size=66);logu=float(np.log(stream.uniform()));cand=approximate_state(cache,x,x,material(x),y)
  answer=resolve_step(cache,current,cand,y,logu,paid);current,cand=answer['current'],answer['candidate'];intervals.append(([current.lower,current.upper],[cand.lower,cand.upper]))
  if answer['accepted']:current=cand;approxaccepted+=not current.exact
  paths.append(current.x.copy());flags.append(answer['accepted'])
 # Explicit eta>=1 case must invoke paid callback, not pretend an exact cache.
 x=current.x+.1;far=approximate_state(cache,x,x,chi0+100,y);before=fullcalls[0]
 answer=resolve_step(cache,current,far,y,-1.,paid);assert fullcalls[0]>before
 return np.array(paths),np.array(flags),np.array(intervals),fullcalls[0],approxaccepted
interval=run_interval()
# Separate stream and no inspection of interval outputs in this function.
def run_full():
 stream=np.random.default_rng(777);x=np.zeros(66);current=exact(x,x,chi0);paths=[x.copy()];flags=[];phis=[]
 for step in range(512):
  x=np.sqrt(1-.002**2)*current.x+.002*stream.normal(size=66);logu=float(np.log(stream.uniform()));cand=exact(x,x,material(x));phis.append([current.lower,cand.lower]);accept=logu<min(0.,current.lower-cand.lower)
  if accept:current=cand
  paths.append(current.x.copy());flags.append(accept)
 return np.array(paths),np.array(flags),np.array(phis)
full=run_full();assert np.array_equal(interval[0],full[0]);assert np.array_equal(interval[1],full[1]);assert interval[4]>0
bound=interval[2];phi=full[2];assert np.all(bound[:,:,0]<=phi+1e-10) and np.all(phi<=bound[:,:,1]+1e-10)
result=dict(status='PASS',scope='512step independent arbitrary12x12 complex toy chains; no Maxwell/GT data',steps=512,
            exact_path_equality=True,accept_flag_mismatches=0,potential_bound_misses=0,accepted_approximate_states=interval[4],
            Full_callback_calls_including_forced_eta_fallback=interval[3],forced_eta_fallback='PASS',
            physics_RHS=0,physics_states=0,training_updates=0,wall_seconds=time.monotonic()-wall,cpu_seconds=time.process_time()-cpu,
            peak_rss_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)*(1 if sys.platform=='darwin' else 1024))
Path('research/p3l_r4_conditional_transport/AUDIT/risk/INTERVAL_CHAIN_MATH_CHECK.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
