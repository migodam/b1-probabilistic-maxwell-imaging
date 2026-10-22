"""Independent analytic controls only: never imports any physics module."""
import json
import os
from pathlib import Path
import sys
import time

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "src"))
import numpy as np
import scipy
from scipy.integrate import quad
from scipy.stats import multivariate_normal
from p3l.kernels import (ChainState, PCN, PropertyReflection, EllipticalSlice,
                        PreconditionedMALA, initialize, reflect_unit_box, run_chain)
from p3l.diagnostics import (rank_normalized_split_rhat, effective_sample_size,
                             material_diagnostics, crosscheck_chains)


def main():
    start = time.monotonic()
    rng = np.random.default_rng(2026102203)
    checks = {}
    mean = np.array([.45, -.6, .2, -.15])
    C = np.array([[.7,.2,0,0],[.2,.6,0,0],[0,0,.8,.1],[0,0,.1,.9]])
    precision = np.linalg.inv(C)
    full_cov = np.eye(64)
    full_cov[:4,:4] = C
    calls = {"target": 0, "gradient": 0}
    def target(q):
        calls["target"] += 1
        if np.any(q[64:] < 0) or np.any(q[64:] > 1):
            return -np.inf
        d = q[:4]-mean
        return -.5*(d@precision@d+q[4:64]@q[4:64])+2*q[64]-3*q[65]
    def grad(q):
        calls["gradient"] += 1
        out = -q[:64].copy()
        out[:4] = -precision@(q[:4]-mean)
        return out

    # Independent quadrature of the bounded exponential laws.
    u_mean, u_variance = [], []
    for t in (2., -3.):
        Z = quad(lambda u: np.exp(t*u), 0., 1.)[0]
        m = quad(lambda u: u*np.exp(t*u), 0., 1.)[0]/Z
        v = quad(lambda u: (u-m)**2*np.exp(t*u), 0., 1.)[0]/Z
        u_mean.append(m)
        u_variance.append(v)
    exact_mean = np.r_[mean, u_mean]
    exact_variance = np.r_[np.diag(C), u_variance]
    results = {}
    chains_by_method = {}
    for name in ("pcn", "elliptical_slice", "mala"):
        chains, accepted, used = [], [], []
        for chain in range(4):
            q = np.r_[rng.normal(size=64), rng.uniform(size=2)]
            state, init = initialize(target, q)
            counts = dict(target=init["target_evaluations"], gradient=0)
            if name == "pcn":
                kernels = [PCN(target, beta=.7, property_sd=.5)]
            elif name == "elliptical_slice":
                kernels = [EllipticalSlice(target), PropertyReflection(target, .5)]
            else:
                kernels = [PreconditionedMALA(target, grad, .5, full_cov),
                           PropertyReflection(target, .5)]
            draws, acc = [], []
            before = dict(calls)
            for i in range(10000):
                for kernel in kernels:
                    transition = kernel.step(state, rng)
                    state = transition.state
                    md = transition.metadata
                    counts["target"] += md["target_evaluations"]
                    counts["gradient"] += md["gradient_evaluations"]
                    if not md["accepted"]:
                        assert md["repeated_state"]
                    acc.append(md["accepted"])
                if i >= 2000:
                    draws.append(state.q[np.r_[0:4,64:66]])
            assert calls["target"]-before["target"] == counts["target"]-1
            assert calls["gradient"]-before["gradient"] == counts["gradient"]
            chains.append(draws)
            accepted.append(float(np.mean(acc)))
            used.append(counts)
        values = np.asarray(chains)
        chains_by_method[name] = values
        d = material_diagnostics(values)
        errors = np.abs(d["mean"]-exact_mean)
        # Stochastic control tolerance fixed here, not a production criterion.
        ok = bool(np.all(errors < .05) and np.all(np.abs(d["variance"]-exact_variance) < .08)
                  and np.all(d["rhat"] < 1.02))
        assert ok, (name, errors, d["rhat"])
        results[name] = {"mean": d["mean"].tolist(), "exact_mean": exact_mean.tolist(),
            "variance": d["variance"].tolist(), "exact_variance": exact_variance.tolist(),
            "mean_abs_error": errors.tolist(), "rhat": d["rhat"].tolist(),
            "bulk_ess": d["bulk_ess"].tolist(), "mean_mcse": d["mean_mcse"].tolist(),
            "acceptance": accepted, "invocation_counts": used,
            "fixed_check_tolerances": {"mean": .05, "variance": .08, "rhat": 1.02}, "passed": ok}
    checks["analytic_gaussian_times_bounded_exponential"] = results

    # Independently calculate normalized proposal densities for MALA's MH.
    q = np.r_[rng.normal(size=64), [.4,.8]]
    state, _ = initialize(target, q)
    mala = PreconditionedMALA(target, grad, .5, full_cov)
    # Reproduce the proposed point with an independent RNG initialized alike.
    seed = 722
    xi = np.random.default_rng(seed).normal(size=64)
    g0 = grad(q)
    qp = q.copy()
    qp[:64] += .25*(full_cov@g0)+np.sqrt(.5)*(np.linalg.cholesky(full_cov)@xi)
    gp = grad(qp)
    reverse = multivariate_normal.logpdf(q[:64], qp[:64]+.25*full_cov@gp, .5*full_cov)
    forward = multivariate_normal.logpdf(qp[:64], q[:64]+.25*full_cov@g0, .5*full_cov)
    expected = target(qp)-target(q)+reverse-forward
    got = mala.step(state, np.random.default_rng(seed)).metadata["log_mh_ratio"]
    assert abs(got-expected) < 1e-10
    checks["mala_normalized_proposal_logratio_error"] = abs(got-expected)

    # Canonical pCN cancellation, repeated rejections and reflection crossings.
    prior = lambda q: -.5*(q[:64]@q[:64])
    s = ChainState(q, prior(q))
    for _ in range(100):
        out = PCN(prior, .9, .5).step(s, rng)
        assert out.metadata["accepted"] and abs(out.metadata["log_mh_ratio"]) < 1e-11
        s = out.state
    reflected = reflect_unit_box(np.array([-7.2,-2.,-.2,0.,1.,1.2,8.8]))
    assert np.allclose(reflected, [.8,0.,.2,0.,1.,.8,.8])
    point_target = lambda q: 0. if np.array_equal(q, np.r_[np.zeros(64), [.5,.5]]) else -np.inf
    r = run_chain(point_target, np.r_[np.zeros(64), [.5,.5]], [PCN(point_target)], 20, rng)
    assert r["states"].shape == (20,66) and np.all(r["states"] == r["states"][0])
    assert all(not t["accepted"] for t in r["transitions"])
    checks["pcn_prior_cancellation_reflection_and_repeated_rejections"] = True

    # Transformed 66D MALA never invokes the canonical property-box check.
    transformed = lambda x: -.5*(x@x)
    transformed_grad = lambda x: -x
    transformed_state, _ = initialize(transformed, rng.normal(size=66))
    transformed_step = PreconditionedMALA(transformed, transformed_grad, .1, dimension=66)
    t = transformed_step.step(transformed_state, rng)
    assert t.metadata["gradient_evaluations"] == 2
    checks["transformed_66_mala_callable_contract"] = True

    independent = rng.normal(size=(4,16000))
    iid_ess = float(effective_sample_size(independent, "mean"))
    assert .75*independent.size < iid_ess < 1.3*independent.size
    ar = np.zeros_like(independent)
    rho = .8
    for i in range(1, ar.shape[1]):
        ar[:,i] = rho*ar[:,i-1]+np.sqrt(1-rho*rho)*independent[:,i]
    ar_ess = float(effective_sample_size(ar, "mean"))
    theoretical = ar.size*(1-rho)/(1+rho)
    assert .75 < ar_ess/theoretical < 1.3
    shifted = independent.copy()
    shifted[0] += 2.
    scale = independent.copy()
    scale[0] *= 4.
    assert rank_normalized_split_rhat(shifted) > 1.1
    assert rank_normalized_split_rhat(scale) > 1.1  # folded Rhat detects scale.
    assert np.isnan(rank_normalized_split_rhat(np.ones((4,100))))
    stuck = np.repeat(np.arange(4.)[:,None],100,axis=1)
    assert np.isinf(rank_normalized_split_rhat(stuck))
    assert np.isnan(material_diagnostics(stuck)["mean_mcse"])
    try:
        material_diagnostics(np.full((4,100), np.nan))
    except ValueError:
        pass
    else:
        raise AssertionError("nonfinite input silently accepted")
    complex_diag = material_diagnostics(independent[:,:1000,None]*(1+2j))
    assert complex_diag["mean"].shape == (1,2)
    cc = crosscheck_chains(chains_by_method["pcn"], chains_by_method["mala"], independent=True, tolerance=.1)
    assert np.all(cc["agreement_within_frozen_tolerance"])
    checks["diagnostics"] = {"iid_mean_ess": iid_ess, "iid_count": independent.size,
        "ar1_rho": rho, "ar1_estimated_ess": ar_ess, "ar1_theoretical_ess": theoretical,
        "shifted_rhat": float(rank_normalized_split_rhat(shifted)),
        "scale_rhat": float(rank_normalized_split_rhat(scale)),
        "crossmechanism_mean_difference": cc["difference"].tolist(),
        "constant_stuck_nonfinite_and_complex_handling": "PASS"}
    report = dict(status="PASS_ANALYTIC_CONTROLS_ONLY", seed=2026102203,
        numpy=np.__version__, scipy=scipy.__version__, wall_seconds=time.monotonic()-start,
        new_physics_calls=0, target_family="ANALYTIC_CONTROL_NOT_P1_POSTERIOR",
        checks=checks,
        cautions=["Finite moment checks do not prove detailed balance or convergence.",
                  "No Maxwell imports, likelihoods, truths, or solver calls.",
                  "No claim about P1 posterior, material risk, equal-cost speed, or qualified teacher.",
                  "Diagnostic thresholds here are control-specific, not frozen production gates."])
    path = Path(__file__).with_name("results.json")
    path.write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps({"status": report["status"], "wall_seconds": report["wall_seconds"], "report": str(path)}))


if __name__ == "__main__":
    main()
