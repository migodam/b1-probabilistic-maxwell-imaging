"""Small shifted-prior linear Gaussian control; no physics backend calls."""
from pathlib import Path
from types import SimpleNamespace
import json
import sys

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from scipy.linalg import solve, cho_solve
from scipy.stats import multivariate_normal
from p3l.core import Point, plain
from p3l.conditional import (make_frame, compose, conditional_reference,
    reference_logpdf, importance_correct, evaluate_slice)


class AnalyticPrior:
    def __init__(self,rng):
        self.A = rng.normal(size=(5,64))/8
        self.config = dict(field_mean=.13)
        self.temperature = .37
        self.host = .1+.04j
        self.mask = np.array([1,1,0,1,1],bool)
        self.real_interval = np.array([.8,1.8])
        self.imag_interval = np.array([.03,.18])
    def decode(self,q):
        field = self.config['field_mean']+self.A@q[:64]
        t = np.clip((field/self.temperature+1)/2,0,1)
        f = t**3*(10-15*t+6*t*t)*self.mask
        delta = (.8+q[64])+1j*(.03+.15*q[65])
        return self.host+delta*f,dict(field=field,fraction=f)


class AnalyticTarget:
    def __init__(self,rng):
        self.prior = AnalyticPrior(rng)
        self.D = rng.normal(size=(9,64))/4
        self.property_D = rng.normal(size=(9,2))
        self.calls = 0
        self.Jcalls = 0
    def evaluate(self,q,purpose):
        self.calls += 1
        return Point(q.copy(),SimpleNamespace(data=self.D@q[:64]+self.property_D@q[64:]),str(self.calls))
    def white(self,r):return r
    def phi(self,p,y):return .5*float(np.sum((p.state.data-y)**2))
    def jacobian(self,p):
        if p.J is None:
            self.Jcalls += 1
            p.J = np.c_[self.D,self.property_D]
        return p.J


def main():
    rng = np.random.default_rng(2026102207)
    target = AnalyticTarget(rng)
    q0 = np.r_[rng.normal(size=64),.3,.8]
    y = rng.normal(size=9)
    frame = make_frame(q0,target.D,3)
    for kind in ('material_generalized','posterior_material_risk'):
        f = make_frame(q0,target.D,3,kind=kind,material_derivative=target.prior.A)
        assert np.linalg.norm(np.c_[f['Vw'],f['Vs']].T@np.c_[f['Vw'],f['Vs']]-np.eye(64)) < 1e-12
    a,u = np.array([.2,-.4,.1]),np.array([.4,.9])
    events = []
    out = conditional_reference(target,y,frame,a,u,max_evaluations=6,callback=events.append)
    Vs = frame['Vs']
    base = q0[:64]+frame['Vw']@a
    Db = target.D@Vs
    H = Db.T@Db+np.eye(Vs.shape[1])
    expected_mean = -solve(H,Db.T@(target.D@base+target.property_D@u-y)+Vs.T@base)
    expected_cov = np.linalg.inv(H)
    ref = out['reference']
    mean_error = float(np.max(np.abs(expected_mean-ref['mean_b'])))
    covariance_error = float(np.max(np.abs(expected_cov-ref['covariance_b'])))
    assert mean_error < 1e-9 and covariance_error < 1e-10
    assert out['counts']['target_evaluations'] == target.calls <= 6
    assert out['counts']['jacobian_evaluations'] == target.Jcalls
    probes = rng.normal(size=(3,Vs.shape[1]))
    logpdf_error = float(np.max(np.abs(reference_logpdf(ref,probes)-multivariate_normal.logpdf(probes,expected_mean,expected_cov))))
    assert logpdf_error < 1e-10
    corrected = importance_correct(target,y,frame,a,u,ref,samples=8,rng=rng,callback=events.append)
    logweight_spread = float(np.ptp(corrected['logweights']))
    normalizer_error = abs(corrected['logZ_estimate']-ref['logZ_GN'])
    assert logweight_spread < 1e-8 and normalizer_error < 1e-8
    # Independent field mean/variance under the supplied exact Gaussian.
    mm = out['material_moments']
    expected_field = .13+target.prior.A@(base+Vs@expected_mean)
    expected_var = np.diag(target.prior.A@Vs@expected_cov@Vs.T@target.prior.A.T)
    assert np.max(np.abs(mm['field_mean']-expected_field)) < 1e-10
    assert np.max(np.abs(mm['field_sd']**2-expected_var)) < 1e-10
    sl = evaluate_slice(target,y,frame,a,u,b=np.zeros(Vs.shape[1]))
    assert sl['target_subtype'] == 'SLICE' and not sl['original_weak_marginal']
    zero = importance_correct(target,y,frame,a,u,ref,samples=2,rng=rng,
                              target_family='LOCAL_C',region=lambda q:False,callback=events.append)
    assert zero['status'] == 'NO_POSITIVE_WEIGHTS' and zero['material_mean'] is None
    # Callback events are directly writable through core.plain without infinities.
    json.dumps(plain(events),allow_nan=False)
    result = dict(status='PASS_ANALYTIC_CONDITIONAL_GAUSSIAN_ONLY',physics_calls=0,
        gaussian_mean_error=mean_error,gaussian_covariance_error=covariance_error,
        normalized_logpdf_error=logpdf_error,constant_logweight_spread=logweight_spread,
        logZ_error=normalizer_error,mode_state_calls=out['counts']['target_evaluations'],
        shifted_prior_center_preserved=True,zero_weight_local_case='PASS',
        nondefault_actual_halfwidth=.37,nonzero_actual_field_mean=.13)
    Path(__file__).with_name('results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))


if __name__ == '__main__':main()
