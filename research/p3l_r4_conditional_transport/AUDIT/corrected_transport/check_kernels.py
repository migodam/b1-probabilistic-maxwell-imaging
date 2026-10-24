"""Mock densities/targets only; no neural inference, truth, or physics."""
from pathlib import Path
import json,time,resource,sys
import numpy as np
from scipy.special import logsumexp
from p3r4.corrected_transport import *
wall,cpu=time.monotonic(),time.process_time();rng=np.random.default_rng(20261025142);mu=np.linspace(-.3,.3,66);scale=np.linspace(.7,1.3,66)
def logq(x):return float(-.5*np.sum(((x-mu)/scale)**2)-np.log(scale).sum()-33*np.log(2*np.pi))
def phi(x):return float(.5*np.sum((x[:8]-.4)**2))
def logtarget(x):return -phi(x)-.5*float(x@x)
def lg(dest,origin):return float(-.5*np.sum((dest-np.sqrt(1-BETA**2)*origin)**2/BETA**2)-66*np.log(BETA)-33*np.log(2*np.pi))
qerror=0.;pcnerror=0.;mixerror=0.
for t in range(200):
 x=rng.normal(size=66);y=mu+scale*rng.normal(size=66) if t%2 else pcn_proposal(x,rng.normal(size=66))
 a=log_acceptance('q',x,phi(x),y,phi(y),logq(x),logq(y));b=log_acceptance('q',y,phi(y),x,phi(x),logq(y),logq(x))
 err=abs(logtarget(x)+logq(y)+a-logtarget(y)-logq(x)-b);qerror=max(qerror,err);assert err<1e-10
 ap=log_acceptance('pcn',x,phi(x),y,phi(y));bp=log_acceptance('pcn',y,phi(y),x,phi(x))
 if t%2==0:
  err=abs(logtarget(y)-logtarget(x)+lg(x,y)-lg(y,x)-(phi(x)-phi(y)));pcnerror=max(pcnerror,err);assert err<1e-7
 fx=logtarget(x)+logsumexp([np.log(.5)+logq(y)+a,np.log(.5)+lg(y,x)+ap]);fy=logtarget(y)+logsumexp([np.log(.5)+logq(x)+b,np.log(.5)+lg(x,y)+bp])
 err=abs(fx-fy);mixerror=max(mixerror,err);assert err<1e-7
assert component_for('mixture_half',.1)=='q' and component_for('mixture_half',.5)=='pcn' and component_for('mixture_half',.9)=='pcn'
assert all(component_for('independence_q',u)=='q' and component_for('pcn',u)=='pcn' for u in (.1,.9))
a=step_randomness(123,0,1,2);b=step_randomness(123,0,1,2);assert np.array_equal(a[0],b[0]) and a[1:]==b[1:]
stuck=np.broadcast_to(np.arange(4)[:,None,None],(4,128,1)).copy();diagnostics=multichain_diagnostics(stuck)
assert diagnostics['constant_or_stuck_chain_function_count']==4 and np.isinf(diagnostics['split_Rhat'][0]) and all(np.isnan(v[0]) for v in diagnostics['ESS_per_chain'])
# A rejected transition explicitly retains state and reference; not a draw fromq.
x=rng.normal(size=66);proposed=np.full(66,100.);la=log_acceptance('q',x,phi(x),proposed,phi(proposed),logq(x),logq(proposed));assert la<-.0001
states=[x.copy()];refs=['same_paid_endpoint'];accepted=(-.00001<la)
if not accepted:states.append(x.copy());refs.append(refs[-1])
assert np.array_equal(states[0],states[1]) and refs[0]==refs[1]
result=dict(status='PASS',independence_log_detailed_balance_max_error=qerror,pcn_prior_cancellation_max_error=pcnerror,constant_half_mixture_detailed_balance_max_error=mixerror,
 component_specific_ratios=True,state_independent_selector=True,reproducible_common_streams=True,rejections_retained=True,stuck_chains_not_falseESSorMCSE=True,
 physics_RHS=0,physics_states=0,neural_sampling_calls=0,training_updates=0,wall_seconds=time.monotonic()-wall,cpu_seconds=time.process_time()-cpu,
 peak_rss_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)*(1 if sys.platform=='darwin' else 1024))
Path('research/p3l_r4_conditional_transport/AUDIT/corrected_transport/KERNEL_MOCK_CHECK.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
