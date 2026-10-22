#!/usr/bin/env python3
"""Reproducible analytical probability controls for Research P3-L.

No Maxwell/PDE solver is imported or called. No network or paid computation.
Outputs: results.json, two PNG/PDF scientific figures, summary.txt.
Run: python run_math_checks.py
"""
from __future__ import annotations

import json
import math
import os
import platform
from pathlib import Path

OUT = Path(__file__).resolve().parent
os.environ.setdefault("MPLCONFIGDIR", str(OUT / ".mpl_cache"))

import numpy as np
import scipy
from scipy.integrate import quad
from scipy.special import ndtr
from scipy.stats import norm
from numpy.polynomial import Polynomial
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

SEED = 20261019
ATOL = 2e-10
SQRT_2PI = math.sqrt(2 * math.pi)


def qintegral(fun, lo=-np.inf, hi=np.inf):
    value, error = quad(fun, lo, hi, epsabs=2e-12, epsrel=2e-12, limit=300)
    return float(value), float(error)


def check_close(actual, expected, atol=ATOL, rtol=ATOL):
    if not np.isclose(actual, expected, atol=atol, rtol=rtol):
        raise AssertionError((actual, expected, abs(actual - expected)))


def profile_check():
    """U=.5*a^2+.5*exp(2*k*a)*(b-a^2)^2 on R^2."""
    k = 0.8
    z_exact = 2 * math.pi * math.exp(0.5 * k * k)

    def b_integrated(a):
        return SQRT_2PI * np.exp(-0.5 * a * a - k * a)

    z_quad, z_error = qintegral(b_integrated)
    ea_num, _ = qintegral(lambda a: a * b_integrated(a))
    ea2_num, _ = qintegral(lambda a: a * a * b_integrated(a))
    ea = ea_num / z_quad
    ea2 = ea2_num / z_quad
    # Conditional b|a is Gaussian centered at a^2. Verify its moment using
    # the standard-normal conditional coordinate, not a fitted posterior.
    ez, _ = qintegral(lambda t: t * norm.pdf(t))
    ez2, _ = qintegral(lambda t: t * t * norm.pdf(t))
    eb = ea2
    eb2_num, _ = qintegral(
        lambda a: (
            a**4 * SQRT_2PI * np.exp(-0.5 * a * a - k * a)
            + SQRT_2PI * np.exp(-0.5 * a * a - 3 * k * a)
        )
    )
    eb2 = eb2_num / z_quad
    eb2_exact = (3 + 6 * k * k + k**4) + math.exp(4 * k * k)
    direct_b_checks = []
    for a in [-3.0, -1.0, 0.0, 1.0, 3.0]:
        precision = math.exp(2 * k * a)
        integ, err = qintegral(
            lambda b: math.exp(-0.5 * a * a - 0.5 * precision * (b - a * a)**2)
        )
        exact = float(b_integrated(a))
        check_close(integ, exact)
        direct_b_checks.append({
            "a": a, "original_b_quadrature": integ,
            "analytic_b_integral": exact,
            "absolute_error": abs(integ-exact), "quad_error_estimate": err,
        })
    check_close(z_quad, z_exact)
    check_close(ea, -k)
    check_close(ea2-ea*ea, 1.0)
    check_close(eb, 1+k*k)
    check_close(eb2, eb2_exact)
    return {
        "definition": "U(a,b)=0.5*a^2+0.5*exp(2*k*a)*(b-a^2)^2; a,b in R",
        "k": k,
        "joint_normalizer_exact": z_exact,
        "joint_normalizer_quadrature": z_quad,
        "joint_normalizer_quad_error_estimate": z_error,
        "profile_a_distribution": {"mean": 0.0, "variance": 1.0},
        "integrated_a_distribution": {"mean": -k, "variance": 1.0},
        "quadrature_Ea": ea,
        "quadrature_Var_a": ea2-ea*ea,
        "quadrature_Eb": eb,
        "exact_Eb": 1+k*k,
        "profile_ridge_Eb": 1.0,
        "quadrature_Eb2": eb2,
        "exact_Eb2": eb2_exact,
        "conditional_standard_normal_checks": {"Ez": ez, "Ez2": ez2},
        "direct_original_b_integrals": direct_b_checks,
        "interpretation": (
            "Profile density omits the strong conditional-volume factor "
            "exp(-k*a); exact weak marginal is N(-k,1)."
        ),
    }


def boundary_case(name, g, H, lower, upper):
    """Integrate exp(-g*b-H*b^2/2) on an interval, retaining linear term."""
    mu = -g/H
    sigma = H**-0.5
    alpha, beta = (lower-mu)/sigma, (upper-mu)/sigma
    mass = float(ndtr(beta)-ndtr(alpha))
    pa, pb = float(norm.pdf(alpha)), float(norm.pdf(beta))
    z = math.sqrt(2*math.pi/H) * math.exp(g*g/(2*H)) * mass
    mean = mu+sigma*(pa-pb)/mass
    a_phi = 0.0 if not np.isfinite(alpha) else alpha*pa
    b_phi = 0.0 if not np.isfinite(beta) else beta*pb
    variance = sigma*sigma * (
        1+(a_phi-b_phi)/mass-((pa-pb)/mass)**2
    )
    kernel = lambda b: np.exp(-g*b-0.5*H*b*b)
    zq, ze = qintegral(kernel, lower, upper)
    m1q, _ = qintegral(lambda b: b*kernel(b), lower, upper)
    m2q, _ = qintegral(lambda b: b*b*kernel(b), lower, upper)
    meanq = m1q/zq
    varq = m2q/zq-meanq*meanq
    check_close(zq, z)
    check_close(meanq, mean)
    check_close(varq, variance)
    # Integrating a Hessian-only Gaussian at a boundary omits g*b.
    z_no_linear = math.sqrt(2*math.pi/H) * (
        ndtr(upper*math.sqrt(H))-ndtr(lower*math.sqrt(H))
    )
    mode = min(max(mu, lower), upper)
    return {
        "name": name, "g": g, "H": H,
        "lower": lower if np.isfinite(lower) else "-infinity",
        "upper": upper if np.isfinite(upper) else "+infinity",
        "unconstrained_mean": mu, "unconstrained_sd": sigma,
        "unconstrained_Gaussian_mass_inside_interval": mass,
        "constrained_mode": mode, "normalizer_exact": z,
        "normalizer_quadrature": zq,
        "normalizer_quad_error_estimate": ze,
        "mean_exact": mean, "mean_quadrature": meanq,
        "variance_exact": variance, "variance_quadrature": varq,
        "normalizer_if_linear_term_is_omitted": float(z_no_linear),
        "relative_normalizer_error_if_linear_term_is_omitted": float(z_no_linear/z-1),
    }


def boundary_checks():
    cases = [
        ("Lower active boundary", 1.25, 3.0, 0.0, np.inf),
        ("Upper active boundary", -2.0, 1.5, -np.inf, 0.0),
        ("Finite positive interval", 0.7, 2.4, 0.0, 1.0),
        ("Two finite boundaries", -0.6, 4.0, -0.8, 0.4),
    ]
    return [boundary_case(*case) for case in cases]


def toy_forward(x):
    return 0.15*x*x+np.sin(x)


def toy_phi(x, observation=0.25, sigma=0.6):
    return 0.5*((toy_forward(x)-observation)/sigma)**2


def phi_taylor_interval(x, degree, observation=0.25, sigma=0.6):
    """Global deterministic sine-Taylor remainder; exact polynomial term."""
    if degree % 2 != 1:
        raise ValueError("Degree must be odd.")
    xt = np.asarray(x, dtype=float)
    sin_approx = np.zeros_like(xt)
    for j in range((degree+1)//2):
        sin_approx += (-1)**j*xt**(2*j+1)/math.factorial(2*j+1)
    approximation = 0.15*xt*xt+sin_approx
    # Taylor polynomial through degree+1 has the same coefficients;
    # Lagrange remainder uses the bounded (degree+2)th sine derivative.
    forward_error = np.abs(xt)**(degree+2)/math.factorial(degree+2)
    abs_r = np.abs((approximation-observation)/sigma)
    radius = forward_error/sigma
    lower = 0.5*np.maximum(0.0, abs_r-radius)**2
    upper = 0.5*(abs_r+radius)**2
    return lower, upper


def certified_mh_check():
    rng = np.random.default_rng(SEED)
    n = 50000
    proposal_sd = 0.8
    x = 0.0
    pair_x = np.empty(n)
    pair_y = np.empty(n)
    log_u = np.empty(n)
    log_prior_ratio = np.empty(n)
    exact_ratio = np.empty(n)
    exact_accept = np.empty(n, dtype=bool)
    for i in range(n):
        y = x+proposal_sd*rng.normal()
        lu = np.log(rng.random())
        lp = -0.5*y*y+0.5*x*x
        ratio = lp+toy_phi(x)-toy_phi(y)
        accept = lu < ratio
        pair_x[i], pair_y[i], log_u[i] = x, y, lu
        log_prior_ratio[i], exact_ratio[i], exact_accept[i] = lp, ratio, accept
        if accept:
            x = y

    degrees = [1, 3, 5, 7, 9, 13]
    intervals = {}
    all_contain = True
    max_violation = 0.0
    for degree in degrees:
        lx, ux = phi_taylor_interval(pair_x, degree)
        ly, uy = phi_taylor_interval(pair_y, degree)
        lo, hi = log_prior_ratio+lx-uy, log_prior_ratio+ux-ly
        violation = np.maximum(lo-exact_ratio, exact_ratio-hi)
        max_violation = max(max_violation, float(np.max(violation)))
        all_contain &= bool(np.all(violation < 1e-10))
        intervals[degree] = (lo, hi)

    initial_reports = []
    for initial_degree in [1, 3, 5, 7]:
        lo0, hi0 = intervals[initial_degree]
        ambiguous0 = (log_u >= lo0) & (log_u <= hi0)
        width = np.maximum(0.0, hi0-lo0)
        conditional_refine_prob = (
            np.exp(np.minimum(0.0, hi0))-np.exp(np.minimum(0.0, lo0))
        )
        conditional_refine_prob = np.clip(conditional_refine_prob, 0, 1)
        if np.any(conditional_refine_prob > np.minimum(1,width)+1e-12):
            raise AssertionError("Conditional refinement bound violated.")
        cumulative_lo = np.full(n, -np.inf)
        cumulative_hi = np.full(n, np.inf)
        decisions = np.zeros(n, dtype=bool)
        resolved = np.zeros(n, dtype=bool)
        levels_used = np.zeros(n, dtype=int)
        refinements = []
        for degree in degrees:
            if degree < initial_degree:
                continue
            lo, hi = intervals[degree]
            cumulative_lo = np.maximum(cumulative_lo, lo)
            cumulative_hi = np.minimum(cumulative_hi, hi)
            newly_accept = (~resolved) & (log_u < cumulative_lo)
            newly_reject = (~resolved) & (log_u > cumulative_hi)
            newly_resolved = newly_accept | newly_reject
            decisions[newly_accept] = True
            levels_used[newly_resolved] = degree
            resolved |= newly_resolved
            refinements.append({
                "taylor_degree": degree,
                "unresolved_after_level": int(np.count_nonzero(~resolved)),
            })
        fallback_count = int(np.count_nonzero(~resolved))
        decisions[~resolved] = exact_accept[~resolved]
        levels_used[~resolved] = -1
        mismatches = int(np.count_nonzero(decisions != exact_accept))
        if mismatches:
            raise AssertionError(f"{mismatches} certified/exact disagreements")
        initial_reports.append({
            "initial_degree": initial_degree,
            "empirical_initial_refine_fraction": float(np.mean(ambiguous0)),
            "average_exact_conditional_refine_probability": float(np.mean(conditional_refine_prob)),
            "martingale_scale_for_refine_fraction": float(
                math.sqrt(np.sum(conditional_refine_prob*(1-conditional_refine_prob)))/n
            ),
            "average_min_one_log_ratio_interval_width": float(np.mean(np.minimum(1,width))),
            "mean_log_ratio_interval_width": float(np.mean(width)),
            "initial_interval_width_p95": float(np.quantile(width, .95)),
            "exact_decision_mismatches": mismatches,
            "full_analytic_fallback_count": fallback_count,
            "full_analytic_fallback_fraction": fallback_count/n,
            "progressive_refinement": refinements,
        })
    return {
        "label": "ANALYTIC_TOY_CONTROL_NOT_MAXWELL",
        "forward": "f(x)=0.15*x^2+sin(x)",
        "observation": 0.25, "noise_sd": 0.6,
        "prior": "N(0,1)", "proposal": "symmetric Gaussian random walk",
        "proposal_sd": proposal_sd, "seed": SEED, "n_steps": n,
        "exact_chain_acceptance_fraction": float(np.mean(exact_accept)),
        "all_log_ratio_intervals_contain_exact_value_to_float_tolerance": all_contain,
        "max_float_containment_violation": max_violation,
        "initial_degrees": initial_reports,
        "important_limit": (
            "Deterministic Taylor bounds are exact analytically; floating arithmetic "
            "checked to 1e-10. This checks the acceptance logic, not a Maxwell "
            "certificate, sampling convergence, speedup, or posterior calibration."
        ),
    }


def smootherstep(f, scale=0.5):
    t = np.clip((np.asarray(f)/scale+1)/2, 0, 1)
    return 6*t**5-15*t**4+10*t**3


def standardized_truncated_moments(lower, upper, max_order):
    vals = np.zeros(max_order+1)
    vals[0] = ndtr(upper)-ndtr(lower)
    if max_order >= 1:
        vals[1] = norm.pdf(lower)-norm.pdf(upper)
    for j in range(2, max_order+1):
        vals[j] = (
            lower**(j-1)*norm.pdf(lower)-upper**(j-1)*norm.pdf(upper)
            +(j-1)*vals[j-2]
        )
    return vals


def gaussian_smootherstep_moment(mu, sd, power, scale=0.5):
    t = Polynomial([0.5, 0.5/scale])
    hp = (6*t**5-15*t**4+10*t**3)**power
    a, b = (-scale-mu)/sd, (scale-mu)/sd
    zm = standardized_truncated_moments(a,b,hp.degree())
    raw = np.zeros(hp.degree()+1)
    for k in range(hp.degree()+1):
        raw[k] = sum(
            math.comb(k,j)*mu**(k-j)*sd**j*zm[j]
            for j in range(k+1)
        )
    exact_from_moments = float(np.dot(hp.coef,raw)+norm.sf(b))
    interior, error = qintegral(
        lambda f: float(smootherstep(f,scale))**power*norm.pdf(f,loc=mu,scale=sd),
        -scale, scale
    )
    tail, tail_error = qintegral(
        lambda f: norm.pdf(f,loc=mu,scale=sd), scale, np.inf
    )
    quadrature = interior+tail
    check_close(exact_from_moments, quadrature, atol=3e-10)
    return {
        "polynomial_truncated_Gaussian_moment_formula": exact_from_moments,
        "quadrature": quadrature,
        "absolute_difference": abs(exact_from_moments-quadrature),
        "quad_error_estimate": error+tail_error,
    }


def prior_checks():
    parameters = [(0.0,.15),(0.0,.5),(0.0,1.0),(.25,.35),(-.4,.6)]
    reports=[]
    for mu,sd in parameters:
        m1=gaussian_smootherstep_moment(mu,sd,1)
        m2=gaussian_smootherstep_moment(mu,sd,2)
        if mu == 0:
            check_close(m1["quadrature"],.5)
        reports.append({
            "Gaussian_field_mean":mu, "Gaussian_field_sd":sd,
            "saturation_scale":.5, "E_phi":m1, "E_phi_squared":m2,
            "Var_phi":m2["quadrature"]-m1["quadrature"]**2,
        })
    chi_b=.1+.04j
    amp_mean=1.3+.105j
    chi_roi=chi_b+.5*amp_mean
    chi_average=chi_b+(256/576)*.5*amp_mean
    return {
        "label":"P1_DECODER_ANALYTIC_CONTROL_NO_MAXWELL",
        "decoder":"phi=1_ROI*h(f/.5), h=6*t^5-15*t^4+10*t^3, t=clip((f/.5+1)/2,0,1)",
        "five_parameter_checks":reports,
        "exact_native_prior": {
            "reason": "For every zero-mean Gaussian marginal f, h(-s)=1-h(s), hence E_phi=.5 in ROI independently of the local field variance.",
            "independence_required":"Uniform material amplitudes are independent of the Gaussian field, as in P1.",
            "E_phi_inside_ROI":.5,
            "E_phi_outside_ROI":0.0,
            "E_a_uniform":1.3,
            "E_b_uniform":.105,
            "E_chi_inside_ROI":{"real":chi_roi.real,"imag":chi_roi.imag},
            "E_chi_outside_ROI":{"real":chi_b.real,"imag":chi_b.imag},
            "E_epsilon_r_inside_ROI":{"real":1+chi_roi.real,"imag":chi_roi.imag},
            "E_epsilon_r_outside_ROI":{"real":1+chi_b.real,"imag":chi_b.imag},
            "ROI_cells":256,"grid_cells":576,
            "spatial_average_of_prior_mean_chi":{"real":chi_average.real,"imag":chi_average.imag},
            "limit":"This is the unconditional prior mean. It is not a posterior mean and does not use scattering data.",
        },
    }


def configure_plots():
    plt.rcParams.update({
        "font.family":"DejaVu Sans",
        "font.size":11,
        "axes.titlesize":13,
        "axes.labelsize":11,
        "legend.fontsize":9.5,
        "figure.facecolor":"white",
        "axes.facecolor":"#fbfcfe",
        "axes.spines.top":False,
        "axes.spines.right":False,
        "axes.grid":True,
        "grid.color":"#dce3eb",
        "grid.alpha":.6,
        "savefig.facecolor":"white",
        "pdf.fonttype":42,
        "ps.fonttype":42,
    })


def figures(results):
    configure_plots()
    blue, orange, green = "#2559a7", "#c77823", "#14856f"
    fig,axs=plt.subplots(1,2,figsize=(11.6,4.4),layout="constrained")
    a=np.linspace(-4,3,600)
    k=results["profile_marginal"]["k"]
    axs[0].plot(a,norm.pdf(a,loc=-k),color=blue,lw=2.7,label=r"Integrated marginal: $N(-0.8,1)$")
    axs[0].plot(a,norm.pdf(a),color=orange,lw=2.4,ls="--",label=r"Normalized profile: $N(0,1)$")
    axs[0].axvline(-k,color=blue,lw=1.2,alpha=.6)
    axs[0].axvline(0,color=orange,lw=1.2,alpha=.6)
    axs[0].set(xlabel="Weak coordinate a",ylabel="Normalized density",
               title="A. Strong-direction volume shifts weak probability",ylim=(0,.54))
    axs[0].legend(loc="upper left",frameon=True,facecolor="white",fontsize=9)
    axs[0].text(.98,.97,
        r"$E[a]=-0.8,\quad E[b]=1.64$"+"\n"
        "Profile ridge: E[b] = 1.00",
        transform=axs[0].transAxes,ha="right",va="top",fontsize=9,
        bbox={"boxstyle":"round,pad=.45","fc":"white","ec":"#ced6df","alpha":.96})
    bcase=results["boundary_gaussian"][0]
    b=np.linspace(-1.9,2.0,650)
    mu,sd=bcase["unconstrained_mean"],bcase["unconstrained_sd"]
    mass=bcase["unconstrained_Gaussian_mass_inside_interval"]
    target=np.where(b>=0,norm.pdf(b,mu,sd)/mass,np.nan)
    axs[1].axvspan(-1.9,0,color="#e1e4e8",alpha=.6,label="Excluded by b >= 0")
    axs[1].plot(b,norm.pdf(b,mu,sd),color=orange,ls="--",lw=2.2,label="Untruncated Gaussian")
    axs[1].plot(b,target,color=green,lw=2.7,label="Correct truncated density")
    axs[1].axvline(0,color="#4d5663",lw=1)
    axs[1].axvline(bcase["mean_exact"],color=green,ls=":",lw=1.8)
    axs[1].set(xlabel="Strong coordinate b",ylabel="Normalized density",
               title="B. Active boundary retains a nonzero linear term",
               xlim=(-1.9,2.0),ylim=(0,2.35))
    axs[1].legend(loc="upper left",frameon=True,facecolor="white",fontsize=9)
    axs[1].text(.98,.97,
        f"Mode = 0; mean = {bcase['mean_exact']:.4f}\n"
        f"Mass retained = {mass:.4f}",
        transform=axs[1].transAxes,ha="right",va="top",fontsize=9,
        bbox={"boxstyle":"round,pad=.45","fc":"white","ec":"#ced6df","alpha":.96})
    fig.suptitle("Analytical probability controls: profile volume and boundary mass",fontsize=15)
    for ext in ["png","pdf"]:
        fig.savefig(OUT/f"01_profile_boundary.{ext}",dpi=190)
    plt.close(fig)

    fig,axs=plt.subplots(1,2,figsize=(11.6,4.5),layout="constrained")
    checks=results["certified_MH"]["initial_degrees"]
    degrees=[r["initial_degree"] for r in checks]
    empirical=np.array([r["empirical_initial_refine_fraction"] for r in checks])
    predicted=np.array([r["average_exact_conditional_refine_probability"] for r in checks])
    bound=np.array([r["average_min_one_log_ratio_interval_width"] for r in checks])
    axs[0].plot(degrees,bound,"^-",color=orange,lw=2,label="Mean interval-width bound")
    axs[0].plot(degrees,predicted,"o-",color=blue,lw=2.2,label="Exact conditional refinement probability")
    axs[0].plot(degrees,empirical,"x",color="#121a24",ms=9,mew=1.6,label="Observed refinement fraction")
    axs[0].set(xlabel="Initial sine Taylor degree",ylabel="Fraction requiring refinement",
               title="C. Certified MH agrees with the exact toy reference",
               xticks=degrees,ylim=(0,max(bound)*1.45))
    axs[0].legend(loc="upper right",frameon=True,facecolor="white")
    axs[0].text(.98,.60,
        "50,000 shared proposals and uniforms\n"
        "0 exact / certified decision mismatches\n"
        r"$P(\mathrm{refine}\mid x,x')\leq\min(1,\ell_+-\ell_-)$",
        transform=axs[0].transAxes,ha="right",va="top",fontsize=10,
        bbox={"boxstyle":"round,pad=.45","fc":"white","ec":"#ced6df","alpha":.96})
    pr=results["P1_smootherstep"]["five_parameter_checks"]
    xs=np.arange(len(pr))
    e1=np.array([r["E_phi"]["polynomial_truncated_Gaussian_moment_formula"] for r in pr])
    e2=np.array([r["E_phi_squared"]["polynomial_truncated_Gaussian_moment_formula"] for r in pr])
    q1=np.array([r["E_phi"]["quadrature"] for r in pr])
    q2=np.array([r["E_phi_squared"]["quadrature"] for r in pr])
    axs[1].bar(xs-.18,e1,width=.33,color=blue,alpha=.85,label=r"Analytic $E[\phi]$")
    axs[1].bar(xs+.18,e2,width=.33,color=green,alpha=.85,label=r"Analytic $E[\phi^2]$")
    axs[1].plot(xs-.18,q1,"kx",ms=5.5,label="Independent quadrature")
    axs[1].plot(xs+.18,q2,"kx",ms=5.5)
    axs[1].axhline(.5,color=orange,ls=":",lw=1.5,label="Symmetric-field mean = 0.5")
    labels=[f"mu={r['Gaussian_field_mean']:g}\nsd={r['Gaussian_field_sd']:g}" for r in pr]
    axs[1].set(xlabel="Gaussian field f (saturation scale = 0.5)",ylabel="Fraction moment",
               title="D. P1 smootherstep moments match quadrature",
               xticks=xs,xticklabels=labels,ylim=(0,1.02))
    axs[1].legend(loc="upper left",frameon=True,facecolor="white",fontsize=9)
    fig.suptitle("Acceptance certification and P1 prior moments — zero Maxwell evaluations",fontsize=15)
    for ext in ["png","pdf"]:
        fig.savefig(OUT/f"02_certification_prior.{ext}",dpi=190)
    plt.close(fig)


def main():
    results={
        "research":"P3-L independent mathematical controls",
        "status":"ALL_CHECKS_PASSED",
        "Maxwell_forward_calls":0,
        "Maxwell_adjoint_calls":0,
        "paid_physics_calls":0,
        "seed":SEED,
        "software":{
            "python":platform.python_version(),"numpy":np.__version__,
            "scipy":scipy.__version__,"matplotlib":matplotlib.__version__,
        },
        "profile_marginal":profile_check(),
        "boundary_gaussian":boundary_checks(),
        "certified_MH":certified_mh_check(),
        "P1_smootherstep":prior_checks(),
        "scope_limit":"Mathematical controls only. No Maxwell posterior, empirical imaging improvement, or solver acceleration claim.",
    }
    (OUT/"results.json").write_text(json.dumps(results,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    figures(results)
    first=results["certified_MH"]["initial_degrees"][0]
    prior=results["P1_smootherstep"]["exact_native_prior"]
    summary=(
        "ALL CHECKS PASSED; Maxwell calls = 0.\n"
        f"Profile example: E[a]={results['profile_marginal']['quadrature_Ea']:.12g}, "
        f"E[b]={results['profile_marginal']['quadrature_Eb']:.12g}; profile E[b]=1.\n"
        f"Boundary cases: {len(results['boundary_gaussian'])}; all exact moments match quadrature.\n"
        "Certified MH: 50,000 shared proposals/uniforms; mismatch count = 0.\n"
        f"Initial degree 1 refinement: {first['empirical_initial_refine_fraction']:.6f}; "
        f"conditional prediction {first['average_exact_conditional_refine_probability']:.6f}; "
        f"width bound {first['average_min_one_log_ratio_interval_width']:.6f}.\n"
        "P1 smootherstep: five Gaussian marginal controls, E[phi] and E[phi^2] match quadrature.\n"
        f"P1 exact unconditional prior mean inside ROI: {prior['E_chi_inside_ROI']}.\n"
        f"P1 exact unconditional prior mean outside ROI: {prior['E_chi_outside_ROI']}.\n"
        "Figures are mathematical control figures, not Maxwell experiments.\n"
    )
    (OUT/"summary.txt").write_text(summary)
    print(summary)


if __name__=="__main__":
    main()
