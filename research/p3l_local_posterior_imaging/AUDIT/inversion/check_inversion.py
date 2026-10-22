"""Tiny exact linear/box optimization check; no physics backend instantiated."""
from pathlib import Path
from types import SimpleNamespace
import json
import shutil
import sys

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from p3l.core import Point
from p3l.inversion import reconstruct


class AnalyticTarget:
    def __init__(self, directory):
        self.book = SimpleNamespace(jobdir=directory)
        self.calls = 0
        self.Jcalls = 0
        self.Gcalls = 0
    def evaluate(self,q,purpose):
        self.calls += 1
        assert np.all(q[64:] >= 0) and np.all(q[64:] <= 1)
        return Point(q.copy(),SimpleNamespace(data=2*q),str(self.calls))
    def white(self,data): return data
    def phi(self,p,y): return float(np.sum((p.state.data-y)**2)/2)
    def jacobian(self,p):
        if p.J is None:
            self.Jcalls += 1
            p.J = 2*np.eye(66)
        return p.J
    def phi_gradient(self,p,y):
        if p.gradient is None:
            self.Gcalls += 1
            p.gradient = 2*(p.state.data-y)
        return p.gradient
    def persist_cache(self,p,path): np.savez(path,q=p.q,prediction=p.state.data)


def main():
    directory = Path(__file__).resolve().parent
    jobs = directory/'analytic_runs'
    if jobs.exists():
        shutil.rmtree(jobs)
    y = np.r_[np.linspace(-2.,3.,64),6.,-2.]
    expected = np.r_[2*y[:64]/5,1.,0.]
    results = {}
    for method in ('full_TRF','coarse_TRF'):
        target = AnalyticTarget(jobs/method)
        p,r = reconstruct(target,y,'analytic_only',method,max_evaluations=48,polish=True)
        error = float(np.max(np.abs(p.q-expected)))
        assert error < 1e-5, (method,error,r)
        assert r['fresh_forward_states'] == target.calls <= 48
        assert r['receiver_J32_calls'] == target.Jcalls
        assert r['fresh_source_gradient_calls'] == target.Gcalls
        assert r['kkt']['scaled_KKT'] < 1e-5
        assert r['failure'] is None
        results[method] = dict(max_parameter_error=error, fresh_states=target.calls,
                              scaled_KKT=r['kkt']['scaled_KKT'],phases=r['phases'])
    # A one-state cap must preserve/persist the finite initial endpoint.
    target = AnalyticTarget(jobs/'cap')
    p,r = reconstruct(target,y,'analytic_cap','full_TRF',max_evaluations=1)
    assert target.calls == 1 and np.allclose(p.q,np.r_[np.zeros(64),.5,.5])
    results['one_state_cap_preserves_start'] = True
    summary = dict(status='PASS_ANALYTIC_ONLY',physics_calls=0,results=results)
    (directory/'results.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary))


if __name__ == '__main__':main()
