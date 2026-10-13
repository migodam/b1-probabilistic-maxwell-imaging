"""Small algebraic checks for B1 next-stage theory. No Maxwell solver calls.
Run: python verify_theory.py --output checks.json
These tests validate finite-dimensional identities and counterexamples, not
Maxwell approximation certificates, sampler mixing, or experimental accuracy.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.special import logsumexp
from scipy.integrate import quad
from scipy.stats import norm


def relative_error(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.linalg.norm(a - b) / max(np.linalg.norm(b), 1e-14))


def main() -> dict:
    rng = np.random.default_rng(20261013)
    checks = {}
    maxerr = 0.0
    for _ in range(40):
        m = 32
        f = rng.normal(size=m)
        y = rng.normal(size=m)
        sigma = 10 ** rng.uniform(-1.5, .2)
        sg = rng.uniform(0., .4)
        r = y - f
        cov = sigma**2 * np.eye(m) + sg**2 * np.outer(f, f)
        sign, ld = np.linalg.slogdet(cov)
        assert sign > 0
        dense = -.5 * (r @ np.linalg.solve(cov, r) + ld + m*np.log(2*np.pi))
        a = sg**2/sigma**2
        den = 1 + a*(f@f)
        fast = -.5 * ((r@r-a*(r@f)**2/den)/sigma**2 + m*np.log(sigma**2)+np.log(den)+m*np.log(2*np.pi))
        maxerr = max(maxerr, abs(float(fast-dense)))
    checks['gain_marginalization'] = {'max_absolute_loglik_error':maxerr, 'pass':maxerr<1e-7}

    n, q, m, ns = 24, 4, 7, 3
    G = (rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))/100
    X = np.diag(rng.uniform(.5,1.5,n)+.1j)
    U = np.eye(n)[:,[1,4,8,19]]
    D = np.diag([.15+.01j,-.2+.02j,.4+.03j,-.1+.01j])
    H = U@D@U.T
    Einc = rng.normal(size=(n,ns))+1j*rng.normal(size=(n,ns))
    S = rng.normal(size=(m,n))+1j*rng.normal(size=(m,n))
    L0 = np.eye(n)-X@G
    L1 = np.eye(n)-(X+H)@G
    jb = np.linalg.solve(L0,X@Einc)
    eb = Einc+G@jb
    j1 = np.linalg.solve(L1,(X+H)@Einc)
    direct = S@(j1-jb)
    finite = S@np.linalg.solve(L1,H@eb)
    RU = np.linalg.solve(L0,U)
    local = S@RU@D@np.linalg.solve(np.eye(q)-U.T@G@RU@D,U.T@eb)
    checks['finite_material_secant'] = {'direct_identity_rel_error':relative_error(finite,direct),'localized_feedback_rel_error':relative_error(local,direct),'pass':relative_error(local,direct)<1e-12}

    jhat = j1+.001*(rng.normal(size=j1.shape)+1j*rng.normal(size=j1.shape))
    residual = (X+H)@Einc-L1@jhat
    Z = np.linalg.solve(L1.conj().T,S.conj().T)
    output_error = S@(j1-jhat)
    checks['dual_residual_identity'] = {'relative_error':relative_error(Z.conj().T@residual,output_error),'pass':relative_error(Z.conj().T@residual,output_error)<1e-11}

    N = 150
    prior = rng.dirichlet(np.ones(N))
    phi = rng.uniform(0,6,N)
    tilt = rng.normal(0,.3,N)
    logpa = np.log(prior)-phi
    pa = np.exp(logpa-logsumexp(logpa))
    logpf = logpa-tilt
    pf = np.exp(logpf-logsumexp(logpf))
    masks = [np.arange(N)%3==i for i in range(3)]
    pA = np.array([pa[x].sum() for x in masks])
    pF = np.array([pf[x].sum() for x in masks])
    k = np.array([np.sum(pa[x]*np.exp(-tilt[x]))/pa[x].sum() for x in masks])
    exact = pA*k/(pA@k)
    lo = np.array([tilt[x].min() for x in masks])
    hi = np.array([tilt[x].max() for x in masks])
    lowmass,highmass = pA*np.exp(-hi),pA*np.exp(-lo)
    lower=np.array([lowmass[i]/(lowmass[i]+highmass.sum()-highmass[i]) for i in range(3)])
    upper=np.array([highmass[i]/(highmass[i]+lowmass.sum()-lowmass[i]) for i in range(3)])
    # The log-odds derivative along phi+t*tilt.
    eps=1e-5
    def probs(t):
        lp=logpa-t*tilt
        p=np.exp(lp-logsumexp(lp))
        return np.array([p[x].sum() for x in masks])
    pp,pm=probs(eps),probs(-eps)
    fd=((np.log(pp[0])-np.log(pp[1]))-(np.log(pm[0])-np.log(pm[1])))/(2*eps)
    expect=np.array([np.sum(pa[x]*tilt[x])/pa[x].sum() for x in masks])
    tv=float(.5*np.sum(abs(pf-pa)))
    tvbound=float(np.tanh((tilt.max()-tilt.min())/4))
    checks['event_probability_bounds']={'reweighting_max_error':float(max(abs(exact-pF))), 'exact_probabilities':pF.tolist(), 'lower_bounds':lower.tolist(), 'upper_bounds':upper.tolist(), 'log_odds_derivative_error':float(abs(fd+expect[0]-expect[1])), 'TV':tv, 'TV_bound':tvbound,'pass':bool(np.all(lower<=pF) and np.all(pF<=upper) and tv<=tvbound and max(abs(exact-pF))<1e-12)}

    r=np.geomspace(.01,.15,40)
    raw=np.sqrt(r**6+r**10)
    prof=r**5
    checks['compensation_counterexample']={'model':'F(r,theta)=(theta+r^3,r^5); null F(0,theta0)=(theta0,0)', 'raw_log_log_slope':float(np.polyfit(np.log(r),np.log(raw),1)[0]),'profiled_log_log_slope':float(np.polyfit(np.log(r),np.log(prof),1)[0]),'interpretation':'Algebraic example, NOT a Maxwell fifth-power claim','pass':abs(np.polyfit(np.log(r),np.log(prof),1)[0]-5)<1e-10}

    width=1e-4
    height=1e4
    samples=rng.random(1024)
    likelihood=np.where(samples<width,height,1.)
    estimated=float(likelihood[samples<width].sum()/likelihood.sum())
    exact=width*height/(1-width+width*height)
    checks['rare_mode_counterexample']={'prior_mass':width,'peak_likelihood':height,'true_event_mass':exact,'no_hit_probability':float((1-width)**1024),'hits':int((samples<width).sum()),'estimated_event_mass':estimated,'ESS':float(likelihood.sum()**2/(likelihood@likelihood)), 'pass':bool((samples<width).sum()==0 and exact>.49)}

    native_var=3.229090594956384e-6
    refined_var=.000733252515315478
    checks['B1_saved_summary_comparison']={'source_commit':'db5aacc9b7f1cc46118a216657df9af8d3d9d1f7','source_files':['results/PRIMARY_SUMMARY.json','results/refinement/SUMMARY.json'],'case':'distinguishable, H=0, Re(epsilon)','native_variance':native_var,'refined_variance':refined_var,'variance_ratio':refined_var/native_var,'native_sd':float(np.sqrt(native_var)),'refined_sd':float(np.sqrt(refined_var)),'sd_ratio':float(np.sqrt(refined_var/native_var)), 'status':'Arithmetic from retrieved saved values; neither is a certified exact reference'}
    native_mean=2.500024354572197
    native_q05=2.4435757578381154
    native_q95=2.5564493291570263
    variance_lower_bound=.05*((native_mean-native_q05)**2+(native_q95-native_mean)**2)
    checks['B1_saved_moment_quantile_coherence']={
        'source':'PRIMARY_SUMMARY.json, distinguishable, H=0, Re(epsilon)',
        'native_mean':native_mean,'native_q05':native_q05,'native_q95':native_q95,
        'reported_variance':native_var,'necessary_variance_lower_bound':variance_lower_bound,
        'lower_bound_divided_by_reported_variance':variance_lower_bound/native_var,
        'can_all_be_exact_summaries_of_one_distribution':False,
        'interpretation':'Detection of incompatible approximate summaries; not a likelihood normalization failure. Test PASS means the counterexample was detected.',
        'pass':bool(native_q05<native_mean<native_q95 and native_var<variance_lower_bound/90.)}

    means_a=np.array([-1.5,2.])
    means_b=np.array([-1.3,3.5])
    weights_a=np.array([.7,.3])
    weights_b=np.array([.4,.6])
    coupling=np.array([[.4,.3],[0.,.3]])
    sigma=.7
    def density(x, means, weights):
        return float(weights@norm.pdf((x-means)/sigma)/sigma)
    true_bayes_error, quad_error=quad(lambda x:.5*min(density(x,means_a,weights_a),density(x,means_b,weights_b)),-15,15,epsabs=1e-11,limit=300)
    distances=abs(means_a[:,None]-means_b[None,:])/sigma
    lower_bound=float(np.sum(coupling*norm.cdf(-distances/2)))
    checks['composite_gaussian_coupling']={'bayes_error_numeric':true_bayes_error,'quad_error_estimate':quad_error,'coupling_lower_bound':lower_bound,'scope':'one-dimensional Gaussian-mixture numerical example, not Maxwell','pass':bool(true_bayes_error+1e-9>=lower_bound)}

    m,pdim,rank=19,13,3
    J=rng.normal(size=(m,pdim))
    Y=rng.normal(size=(m,rank))
    Gamma=np.eye(rank)+.05*rng.normal(size=(rank,rank))
    Nmat=rng.normal(size=(rank,pdim))
    D=Y@np.linalg.solve(Gamma,Nmat)
    step=.02*rng.normal(size=pdim)
    r0=rng.normal(size=m)
    r_a=r0+J@step
    r_b=r0+(J+D)@step
    direct_tilt=.5*(r_b@r_b-r_a@r_a)
    signature_tilt=(Y.T@r0+Y.T@J@step)@np.linalg.solve(Gamma,Nmat@step)+.5*np.linalg.norm(D@step)**2
    kappas=[.6,.9]
    v0=rng.normal(size=m);v0/=np.linalg.norm(v0)
    v1=rng.normal(size=m);v1/=np.linalg.norm(v1)
    R0=.5*kappas[0]*np.linalg.norm(step)**2*v0
    R1=.5*kappas[1]*np.linalg.norm(step)**2*v1
    nonlinear_tilt=.5*(np.linalg.norm(r_b+R1)**2-np.linalg.norm(r_a+R0)**2)
    remainder=r_b@R1-r_a@R0+.5*(R1@R1-R0@R0)
    remainder_bound=.5*(kappas[1]*np.linalg.norm(r_b)+kappas[0]*np.linalg.norm(r_a))*np.linalg.norm(step)**2+.125*sum(k*k for k in kappas)*np.linalg.norm(step)**4
    checks['A17_affine_tilt_and_nonlinear_remainder']={'affine_absolute_error':float(abs(direct_tilt-signature_tilt)),'remainder_identity_error':float(abs(nonlinear_tilt-direct_tilt-remainder)),'actual_remainder_abs':float(abs(remainder)),'remainder_bound':float(remainder_bound),'scope':'random matrices and explicit bounded quadratic remainder','pass':bool(abs(direct_tilt-signature_tilt)<1e-10 and abs(remainder)<=remainder_bound and abs(nonlinear_tilt-direct_tilt-remainder)<1e-10)}

    return {'date':'2026-10-13','seed':20261013,'new_physical_evaluations':0,'neural_training':False,'scope':'algebra and saved-number arithmetic only','checks':checks, 'all_asserted_checks_pass':all(v.get('pass',True) for v in checks.values())}

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,default=Path('checks.json'))
    args=parser.parse_args()
    result=main()
    args.output.write_text(json.dumps(result,indent=2,default=lambda x: x.item() if isinstance(x,np.generic) else str(x))+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2,default=lambda x: x.item() if isinstance(x,np.generic) else str(x)))
    if not result['all_asserted_checks_pass']:
        raise SystemExit('A numerical identity check failed.')
