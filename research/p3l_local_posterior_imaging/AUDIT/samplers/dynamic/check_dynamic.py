"""Normalized variable-metric Gaussian controls, never a physics target."""
from pathlib import Path
import json
import sys
import time

ROOT = Path(__file__).resolve().parents[5]
sys.path.insert(0,str(ROOT/'src'))
import numpy as np
import scipy
from scipy.linalg import cho_solve, solve_triangular
from scipy.stats import multivariate_normal
from p3l.kernels import ChainState, initialize, PositionDependentMALA
from p3l.diagnostics import material_diagnostics


def main():
    start = time.monotonic()
    counts = dict(target=0,gradient=0,geometry=0)
    direction = np.sin(np.arange(66)+1.)/np.sqrt(66)
    def target(x):
        counts['target'] += 1
        return -.5*float(x@x)
    def gradient(x):
        counts['gradient'] += 1
        return -x
    def precision(x):
        counts['geometry'] += 1
        return np.diag(1+.15*np.tanh(x)**2)+.03*np.outer(direction,direction)

    # Reproduce one proposal and evaluate both normalized densities independently.
    rng = np.random.default_rng(2026102211)
    q = rng.normal(size=66)
    state,_ = initialize(target,q)
    h,seed = .6,4881
    H = precision(q)
    L = np.linalg.cholesky(H)
    covariance = np.linalg.inv(H)
    drift = -.5*h*covariance@q
    proposal = q+drift+np.sqrt(h)*solve_triangular(L.T,np.random.default_rng(seed).normal(size=66),lower=False)
    Hproposal = precision(proposal)
    covariance_proposal = np.linalg.inv(Hproposal)
    reverse_drift = -.5*h*covariance_proposal@proposal
    logforward = multivariate_normal.logpdf(proposal,q+drift,h*covariance)
    logreverse = multivariate_normal.logpdf(q,proposal+reverse_drift,h*covariance_proposal)
    expected_ratio = target(proposal)-target(q)+logreverse-logforward
    kernel = PositionDependentMALA(target,gradient,h,precision)
    transition = kernel.step(state,np.random.default_rng(seed))
    md = transition.metadata
    errors = dict(logforward=abs(md['log_proposal_forward']-logforward),
                  logreverse=abs(md['log_proposal_reverse']-logreverse),
                  mh_ratio=abs(md['log_mh_ratio']-expected_ratio))
    assert max(errors.values()) < 1e-10
    assert md['geometry_evaluations'] == 2
    assert abs(md['logdet_precision_proposed']-md['logdet_precision_current']) > 1e-4

    # Exact support rejection reuses both current gradient and geometry.
    exact_point = q.copy()
    support_calls = dict(target=0,gradient=0,geometry=0)
    def point_target(x):
        support_calls['target'] += 1
        return 0. if np.array_equal(x,exact_point) else -np.inf
    def point_gradient(x):
        support_calls['gradient'] += 1
        return -x
    def point_precision(x):
        support_calls['geometry'] += 1
        return np.eye(66)
    reject = PositionDependentMALA(point_target,point_gradient,.2,point_precision)
    state = ChainState(exact_point,0.)
    first = reject.step(state,rng)
    second = reject.step(first.state,rng)
    assert first.metadata['repeated_state'] and second.metadata['repeated_state']
    assert first.metadata['geometry_evaluations'] == 1 and second.metadata['geometry_evaluations'] == 0
    assert first.metadata['gradient_evaluations'] == 1 and second.metadata['gradient_evaluations'] == 0
    assert support_calls == dict(target=2,gradient=1,geometry=1)

    # Four independent overdispersed chains; every rejection is retained.
    arrays, chain_records = [], []
    callback_events = []
    for chain in range(4):
        chain_rng = np.random.default_rng(2026102220+chain)
        state,_ = initialize(target,chain_rng.normal(size=66)*1.7)
        kernel = PositionDependentMALA(target,gradient,.6,precision,
            callback=callback_events.append if chain == 0 else None)
        before = dict(counts)
        accumulated = dict(target=0,gradient=0,geometry=0)
        draws, accepted = [], 0
        for i in range(10000):
            out = kernel.step(state,chain_rng)
            if not out.metadata['accepted']:
                assert np.array_equal(out.state.q,state.q)
                assert out.state.gradient is not None
            state = out.state
            accepted += int(out.metadata['accepted'])
            for key in accumulated:
                accumulated[key] += out.metadata[key+'_evaluations']
            if i >= 2000:
                draws.append(state.q.copy())
        for key in accumulated:
            assert accumulated[key] == counts[key]-before[key]
        assert len(kernel._precision_cache) <= 2
        arrays.append(draws)
        chain_records.append(dict(acceptance=accepted/10000,invocations=accumulated))
    diagnostic = material_diagnostics(np.asarray(arrays))
    mean_error = float(np.max(np.abs(diagnostic['mean'])))
    variance_error = float(np.max(np.abs(diagnostic['variance']-1)))
    max_rhat = float(np.max(diagnostic['rhat']))
    assert mean_error < .08 and variance_error < .12 and max_rhat < 1.02
    assert sum(e['geometry_evaluations'] for e in callback_events) == chain_records[0]['invocations']['geometry']
    result = dict(status='PASS_ANALYTIC_VARIABLE_METRIC_CONTROLS_ONLY',physics_calls=0,
        seed=2026102211,numpy=np.__version__,scipy=scipy.__version__,
        normalized_density_errors=errors,chains=chain_records,
        production_draws_per_chain=8000,chain_count=4,max_mean_error=mean_error,
        max_variance_error=variance_error,max_rhat=max_rhat,
        min_bulk_ess=float(np.min(diagnostic['bulk_ess'])),
        max_mean_mcse=float(np.max(diagnostic['mean_mcse'])),
        cache_and_rejection_control='PASS',separate_invocation_and_callback_counts='PASS',
        wall_seconds=time.monotonic()-start,
        fixed_control_tolerances=dict(max_mean_error=.08,max_variance_error=.12,max_rhat=1.02),
        caution='No Maxwell geometry, posterior, speedup, global mixing, ideal diffusion or scientific acceptance claim.')
    Path(__file__).with_name('results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))


if __name__ == '__main__':main()
