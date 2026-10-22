"""Independent SMC probability controls; no data reads or physics calls."""
from pathlib import Path
import json
import sys
import time

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from scipy.stats import multivariate_normal
from p3l.tempered import (tempered_smc,choose_beta,conditional_ess,systematic_resample)


def main():
    start = time.monotonic()
    rng = np.random.default_rng(2026102217)
    y = np.array([.7,-1.1,.35])
    sigma = .7
    variance = sigma*sigma
    exact_mean = y/(1+variance)
    exact_cov = variance/(1+variance)
    exact_evidence = multivariate_normal.logpdf(y,np.zeros(3),(1+variance)*np.eye(3))
    all_runs = []
    for repeat in range(8):
        counts = dict(prior=0,likelihood=0,gradient=0,precision=0)
        events = []
        def prior(x):
            counts['prior'] += 1
            return -.5*(x@x)-1.5*np.log(2*np.pi)
        def likelihood(x):
            counts['likelihood'] += 1
            return -.5*np.sum((x-y)**2)/variance-1.5*np.log(2*np.pi*variance)
        def gradient(x):
            counts['gradient'] += 1
            return -(x-y)/variance
        def precision(x,beta):
            counts['precision'] += 1
            return np.diag((1+beta/variance)*(1+.05*np.tanh(x)**2))
        out = tempered_smc(lambda r,n:r.normal(size=(n,3)),prior,likelihood,
            rng=np.random.default_rng(2026102230+repeat),N=1024,mutation_steps=3,
            mutation='mala',step_size=.8,gaussian_dimension=3,
            gradloglikelihood=gradient,precision=precision,
            callback=events.append if repeat == 0 else None)
        assert out['beta'] == 1 and out['status'] == 'COMPLETED_FINITE_N_UNQUALIFIED'
        assert abs(out['weights'].sum()-1) < 1e-12
        assert out['counters']['likelihood_evaluations'] == counts['likelihood']
        assert out['counters']['likelihood_gradient_evaluations'] == counts['gradient']
        assert out['counters']['precision_evaluations'] == counts['precision']
        assert out['counters']['logprior_evaluations'] == counts['prior']
        assert abs(sum(h['log_evidence_increment'] for h in out['history'])-out['log_evidence']) < 1e-12
        roots = np.arange(1024)
        for history in out['history']:
            roots = roots[history['parent_indices']]
            assert np.array_equal(roots,history['ancestor_ids'])
            assert history['conditional_ess'] >= .8*1024-1e-8
        assert np.array_equal(out['ancestor_ids'],roots)
        weighted_mean = np.sum(out['weights'][:,None]*out['particles'],axis=0)
        weighted_variance = np.sum(out['weights'][:,None]*(out['particles']-weighted_mean)**2,axis=0)
        all_runs.append(dict(mean=weighted_mean.tolist(),variance=weighted_variance.tolist(),
            log_evidence=out['log_evidence'],stages=len(out['history']),
            weight_ess=out['weight_ess'],counts=out['counters']))
        if repeat == 0:
            assert sum(e['event'] == 'LIKELIHOOD_EVALUATION_ATTEMPT' for e in events) == counts['likelihood']
    pooled_mean = np.mean([r['mean'] for r in all_runs],axis=0)
    pooled_variance = np.mean([r['variance'] for r in all_runs],axis=0)
    evidence_mean = np.mean([r['log_evidence'] for r in all_runs])
    assert np.max(np.abs(pooled_mean-exact_mean)) < .05
    assert np.max(np.abs(pooled_variance-exact_cov)) < .05
    assert abs(evidence_mean-exact_evidence) < .08

    # Prior invariance controls at a constant likelihood: evidence exactly one.
    invariance = {}
    for mutation in ('random_walk','pcn','mala'):
        out = tempered_smc(lambda r,n:r.normal(size=(n,3)),lambda x:-.5*(x@x),lambda x:0.,
            rng=np.random.default_rng(2026102250+len(invariance)),N=8192,mutation_steps=2,
            mutation=mutation,proposal_scale=.7,step_size=.6,gaussian_dimension=3,
            gradloglikelihood=lambda x:np.zeros(3),
            precision=lambda x,beta:np.diag(1+.2*np.tanh(x)**2))
        mean,diagvar = out['particles'].mean(axis=0),out['particles'].var(axis=0)
        assert np.max(np.abs(mean)) < .06 and np.max(np.abs(diagvar-1)) < .08
        assert abs(out['log_evidence']) < 1e-12 and len(out['history']) == 1
        invariance[mutation] = dict(mean=mean.tolist(),variance=diagvar.tolist())
    def native_draw(r,n):return np.c_[r.normal(size=(n,2)),r.uniform(size=(n,2))]
    def native_prior(x):
        if np.any(x[2:] < 0) or np.any(x[2:] > 1):return -np.inf
        return -.5*(x[:2]@x[:2])
    native = tempered_smc(native_draw,native_prior,lambda x:0.,rng=rng,N=8192,
        mutation='pcn',proposal_scale=.7,gaussian_dimension=2,mutation_steps=2)
    assert np.max(np.abs(native['particles'][:,2:].mean(axis=0)-.5)) < .03
    assert np.max(np.abs(native['particles'][:,2:].var(axis=0)-1/12)) < .02
    assert np.all((native['particles'][:,2:] >= 0)&(native['particles'][:,2:] <= 1))
    assert np.all(systematic_resample(np.array([0.,0.,0.,1.]),rng) == 3)
    beta_choice = choose_beta(0.,np.full(4,-np.log(4)),np.array([0.,-1e20,-2e20,-3e20]))
    assert beta_choice['status'] == 'TEMPERING_STALLED_BELOW_BETA_FLOOR'
    zero = tempered_smc(lambda r,n:r.normal(size=(n,1)),lambda x:-.5*(x@x),lambda x:-np.inf,
        rng=rng,N=8,mutation_steps=0)
    assert zero['beta'] == 0 and zero['status'] == 'TEMPERING_STALLED_BELOW_BETA_FLOOR'
    truncated = tempered_smc(lambda r,n:r.normal(size=(n,1)),lambda x:-.5*(x@x),
        lambda x:-10*(x[0]-2)**2,rng=rng,N=32,mutation_steps=0,max_stages=1)
    assert truncated['status'] == 'MAXIMUM_TEMPERING_STAGES_REACHED' and truncated['beta'] < 1
    report = dict(status='PASS_ANALYTIC_SMC_CONTROLS_ONLY',physics_calls=0,
        seed=2026102217,repeat_count=8,particles_per_run=1024,
        exact_mean=exact_mean.tolist(),pooled_mean=pooled_mean.tolist(),
        exact_variance=exact_cov,pooled_variance=pooled_variance.tolist(),
        exact_log_evidence=float(exact_evidence),mean_estimated_log_evidence=float(evidence_mean),
        runs=all_runs,prior_invariance=invariance,native_uniform_property_control='PASS',
        systematic_genealogy_counts_and_callbacks='PASS',beta_floor_and_stage_failures='PASS',
        wall_seconds=time.monotonic()-start,
        caution='Adaptive finite-N bias and missed modes remain; no Maxwell or posterior-reference qualification.')
    Path(__file__).with_name('results.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('runs','prior_invariance')}))


if __name__ == '__main__':main()
