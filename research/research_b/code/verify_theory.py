"""Deterministic algebra/probability checks. Zero Maxwell solver calls."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.stats import norm, truncnorm
from acquisition_core import (bayes_risk, discrete_utility, effective_information,
    posterior_update, shared_gain_loglik, certified_interval_winner,
    gaussian_bank_utility, gaussian_equal_covariance_tv)

rng = np.random.default_rng(26101402)
results = {'scope': 'algebra_and_finite_probability_only', 'maxwell_solver_calls': 0, 'checks': {}}

def record(name, result, details):
    results['checks'][name] = {'passed': bool(result), 'details': details}
    assert result, name

# Exact nonlinear current secant identity (generic complex algebra, NOT a Maxwell experiment).
n = 8
g = (rng.normal(size=(n,n)) + 1j*rng.normal(size=(n,n))) * .02
a1, a2 = np.diag(rng.uniform(.1,.4,n)), np.diag(rng.uniform(.1,.4,n))
e = rng.normal(size=n)+1j*rng.normal(size=n)
l1, l2 = np.eye(n)-a1@g, np.eye(n)-a2@g
j1, j2 = np.linalg.solve(l1, a1@e), np.linalg.solve(l2,a2@e)
predicted = np.linalg.solve(l1,(a1-a2)@(e+g@j2))
err = float(np.linalg.norm(j1-j2-predicted))
record('exact_current_secant', err < 1e-12, {'absolute_error': err})

# Individually uninformative calibration + signal become informative together.
b0,n0 = np.array([[1.]]),np.array([[1.]])
ba,na = np.array([[0.]]),np.array([[1.]])
separate = effective_information(b0,n0)+effective_information(ba,na)
stacked = effective_information(np.vstack([b0,ba]),np.vstack([n0,na]))
record('stack_then_eliminate_shared_nuisance', abs(separate.item()) < 1e-12 and
       abs(stacked.item()-.5)<1e-12, {'separately_projected': separate.item(),'joint':stacked.item()})

# Shared compensation quadratic conflict = matrix parallel sum.
x,y = rng.normal(size=(4,4)),rng.normal(size=(4,4))
h0,ha = x.T@x+np.eye(4),y.T@y+np.eye(4)
v0,va = rng.normal(size=4),rng.normal(size=4)
v = np.linalg.solve(h0+ha,h0@v0+ha@va)
actual = (v-v0)@h0@(v-v0)+(v-va)@ha@(v-va)
formula = (v0-va)@np.linalg.solve(np.linalg.inv(h0)+np.linalg.inv(ha),v0-va)
record('compensation_conflict_parallel_sum', abs(actual-formula)<1e-10,
       {'joint_penalty':float(actual),'formula':float(formula)})

# Proper truncated shared-gain marginal, and conditional increment identity.
f0,fnew = np.array([1.,.3]),np.array([-.7, .8])
y0,ynew = np.array([1.12,.29]),np.array([-.83,.93])
sigma=.12
joint = shared_gain_loglik(np.r_[y0,ynew],np.r_[f0,fnew],sigma)
old = shared_gain_loglik(y0,f0,sigma)
prior = lambda gain: truncnorm.pdf(gain,-2.,2.,loc=1.,scale=.15)
rawold = quad(lambda gain: np.prod(norm.pdf(y0,loc=gain*f0,scale=sigma))*prior(gain),.7,1.3,
              epsabs=1e-11,epsrel=1e-11)[0]
rawjoint = quad(lambda gain: np.prod(norm.pdf(np.r_[y0,ynew],loc=gain*np.r_[f0,fnew],scale=sigma))*prior(gain),.7,1.3,
                epsabs=1e-11,epsrel=1e-11)[0]
conditional = np.log(rawjoint/rawold)
wrong = shared_gain_loglik(ynew,fnew,sigma)
record('shared_gain_conditional_likelihood', abs(joint-np.log(rawjoint))<1e-9 and
       abs(joint-old-conditional)<1e-9 and abs(wrong-conditional)>1e-3,
       {'correct_log_increment':float(conditional),'independent_gain_wrong_log_increment':float(wrong)})

# Entropy and expected classification risk can choose different actions.
w = np.array([.49,.49,.02])
loss = 1-np.eye(3)
# A reveals the rare state perfectly, but does not distinguish the two common states.
# B weakly discriminates the common states, and does not directly reveal the rare state.
ca=np.array([[1.,0.],[1.,0.],[0.,1.]])
cb=np.array([[.55,.45],[.45,.55],[.5,.5]])
ua=discrete_utility(w,ca,loss,np.arange(3),3)
ub=discrete_utility(w,cb,loss,np.arange(3),3)
record('event_eig_vs_task_risk_counterexample', ua['event_eig_nats']>ub['event_eig_nats'] and
       ua['risk_reduction']<ub['risk_reduction'],
       {'A_eig':ua['event_eig_nats'],'B_eig':ub['event_eig_nats'],
        'A_risk_gain':ua['risk_reduction'],'B_risk_gain':ub['risk_reduction']})
record('posterior_probability_martingale', np.allclose(ua['posterior_mean_weights'],w) and
       np.allclose(ub['posterior_mean_weights'],w), {'initial_weights': w.tolist(),
       'expected_updated_weights':ub['posterior_mean_weights'].tolist()})

# Fixed-prior expected optimal loss is Lipschitz in joint predictive TV.
max_ratio=0.
for _ in range(100):
    ww=rng.dirichlet(np.ones(5)); c=rng.dirichlet(np.ones(4),size=5)
    d=rng.dirichlet(np.ones(4),size=5); ll=rng.uniform(0,1,(7,5))
    rr1=discrete_utility(ww,c,ll)['expected_future_risk']
    rr2=discrete_utility(ww,d,ll)['expected_future_risk']
    tv=.5*np.sum(ww[:,None]*abs(c-d)); bound=ll.max()*tv
    max_ratio=max(max_ratio,abs(rr1-rr2)/bound)
    assert abs(rr1-rr2)<=bound+1e-12
record('risk_perturbation_bound_100_trials', max_ratio<=1, {'max_observed_ratio':float(max_ratio)})

record('constant_likelihood_no_update', np.allclose(posterior_update(w,np.zeros(3)),w),{})
record('utility_interval_abstention', certified_interval_winner([.2,.19],[.02,.02]) is None and
       certified_interval_winner([.3,.1],[.01,.01])==0,{})

# Two Gaussian worlds: finite-bank MC risk versus the known binary error.
m=10000
u=rng.uniform(size=m); z=rng.normal(size=(m,1))
stats=gaussian_bank_utility(np.array([.5,.5]),np.array([[-1.],[1.]]),np.eye(1),1-np.eye(2),
                            np.arange(2),2,u,z)
exact=float(norm.cdf(-1))
record('gaussian_predictive_risk', abs(stats.expected_future_risk-exact)<.01,
       {'mc_risk':stats.expected_future_risk,'exact_risk':exact,
        'mc_standard_error':stats.mc_standard_error,'hoeffding_halfwidth':stats.hoeffding_halfwidth})

# No assumption that each noisy observation improves true-event probability.
post=posterior_update(np.array([.5,.5]),np.array([-2.,0.]))
record('realized_true_probability_can_decrease',post[0]<.5,{'true_state_0_probability':float(post[0])})

# Gaussian hard-pair numbers are illustrations only, NOT composite event error estimates.
distances=[.9001699997509817,2.4679439214721497,.107504]
results['hard_pair_illustrations']=[{'distance':d,'equal_prior_simple_error':float(norm.cdf(-d/2))}
                                    for d in distances]
results['all_passed']=all(v['passed'] for v in results['checks'].values())
path=Path(__file__).resolve().parents[1]/'checks'/'THEORY_CHECKS.json'
path.write_text(json.dumps(results,indent=2),encoding='utf-8')
print(json.dumps(results,indent=2))
