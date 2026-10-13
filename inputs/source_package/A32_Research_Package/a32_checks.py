"""A32: deterministic algebra checks and fresh, analytic 2-D TM examples.
No training, no user truth arrays, no historical evaluation-set access.
Run: OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python a32_checks.py
NumPy/SciPy required. All phasors use exp(-i omega t). Every analytic
material-frequency call is counted, including controls and repeated checks.
"""
from __future__ import annotations
import json, math, platform, time
from pathlib import Path
import numpy as np
import scipy
from scipy.special import jv, yv, jvp, yvp, hankel1, h1vp, ndtr, logsumexp

SEED=321012
rng=np.random.default_rng(SEED)
checks={}
physical_calls=0
physical_limit=80
LEDGER=Path(__file__).with_name('A32_RESOURCE_LEDGER.json')
if LEDGER.exists():
    resource_ledger=json.loads(LEDGER.read_text())
else:
    resource_ledger={'cap':80,'charged_calls':0,'history':[]}

max_boundary_residual=0.0

def check(name, passed, **data):
    checks[name]={'passed':bool(passed),**data}

def normalize(logw):
    return np.exp(logw-logsumexp(logw))

def realify(z):
    return np.concatenate((np.asarray(z).real,np.asarray(z).imag),axis=-1)

def cylinder(eps:complex, a:float, k:float, *, R:float=.35, Robs:float=2., nmax:int=16):
    """Exact separated Helmholtz/TM E_z solution, truncated Fourier series.
    mu=1 everywhere; circular shell eps with vacuum cavity a and exterior.
    Boundary continuity: E_z and radial derivative. Plane wave exp(i*k*x).
    Returns eight complex samples of the scattered E_z on Robs.
    A scalar infinite-cylinder reduction, NOT a 3-D vector experiment.
    """
    global physical_calls,max_boundary_residual
    if resource_ledger['charged_calls']>=physical_limit: raise RuntimeError('Cumulative analytic-call budget exhausted; use a declared NEW output directory for independent reproduction')
    physical_calls+=1
    resource_ledger['charged_calls']+=1
    resource_ledger['history'].append({'eps':[complex(eps).real,complex(eps).imag],'a':a,'k':k,'nmax':nmax})
    LEDGER.write_text(json.dumps(resource_ledger,indent=2)+'\n')
    if not 0<=a<R: raise ValueError('Require 0<=a<R')
    ks=k*np.sqrt(complex(eps)); ang=np.arange(8)*2*np.pi/8
    out=np.zeros(8,dtype=complex)
    for n in range(-nmax,nmax+1):
        if a==0:
            A=np.array([[-hankel1(n,k*R),jv(n,ks*R)],[-k*h1vp(n,k*R),ks*jvp(n,ks*R)]],complex)
            b=np.array([jv(n,k*R),k*jvp(n,k*R)],complex)
        else:
            A=np.array([
                [-hankel1(n,k*R),jv(n,ks*R),yv(n,ks*R),0],
                [-k*h1vp(n,k*R),ks*jvp(n,ks*R),ks*yvp(n,ks*R),0],
                [0,jv(n,ks*a),yv(n,ks*a),-jv(n,k*a)],
                [0,ks*jvp(n,ks*a),ks*yvp(n,ks*a),-k*jvp(n,k*a)]],complex)
            b=np.array([jv(n,k*R),k*jvp(n,k*R),0,0],complex)
        # Equilibration reduces cylindrical-wave dynamic range. It is an exact
        # invertible change of variables, not a physical regularization.
        scale=np.maximum(np.max(np.abs(A),axis=0),1e-300)
        Ae=A/scale
        row=np.maximum(np.max(np.abs(Ae),axis=1),1e-300)
        Ae=Ae/row[:,None]; be=b/row
        z=np.linalg.solve(Ae,be); sol=z/scale
        rel=np.linalg.norm(A@sol-b)/(np.linalg.norm(A)*np.linalg.norm(sol)+np.linalg.norm(b)+1e-300)
        max_boundary_residual=max(max_boundary_residual,float(rel))
        out+=(1j**n)*sol[0]*hankel1(n,k*Robs)*np.exp(1j*n*ang)
    return out

def algebra():
    # Exact paired signature of a FROZEN likelihood change, arbitrary finite prior.
    m,p,r=7,9,2
    J=rng.normal(size=(m,p))*.25; Y=rng.normal(size=(m,r))*.2
    Gamma=np.array([[1.3,.2],[-.1,1.1]]); N=rng.normal(size=(r,p))*.25
    T=np.linalg.solve(Gamma,N); J1=J+Y@T; residual=rng.normal(size=m)
    S=rng.normal(size=(150,p)); prior=normalize(rng.normal(size=150))
    phi0=.5*np.sum((residual+S@J.T)**2,axis=1)
    phi1=.5*np.sum((residual+S@J1.T)**2,axis=1)
    a=S@N.T; b=S@J.T@Y; c=a@np.linalg.inv(Gamma).T
    h=np.sum((residual@Y+b)*c,axis=1)+.5*np.sum((c@Y.T)**2,axis=1)
    p0=normalize(np.log(prior)-phi0); p1=normalize(np.log(prior)-phi1)
    p1_signature=normalize(np.log(p0)-h)
    check('paired_likelihood_signature',np.max(np.abs(h-(phi1-phi0)))<1e-11,
          max_potential_error=float(np.max(np.abs(h-(phi1-phi0)))),max_weight_error=float(np.max(np.abs(p1-p1_signature))),
          ambient_material_dimension=p,signature_dimension=2*r)
    # Gaussian mean/covariance changes live in C0[J0^T Y,N^T].
    Lambda=np.diag(1+rng.random(p)); H0=J.T@J+Lambda; H1=J1.T@J1+Lambda
    C0=np.linalg.inv(H0); C1=np.linalg.inv(H1)
    m0=-C0@J.T@residual; m1=-C1@J1.T@residual
    Q=C0@np.column_stack((J.T@Y,N.T)); U=np.linalg.svd(Q,full_matrices=False)[0]
    projection=np.eye(p)-U@U.T
    err=max(np.linalg.norm(projection@(m1-m0)),np.linalg.norm(projection@(C1-C0)))
    check('paired_gaussian_support',err<1e-11,outside_support_norm=float(err))
    # Non-Gaussian conditional invariance, event independent of signature.
    pu=np.array([.15,.35,.50]); ph=.3
    joint=np.column_stack((pu*(1-ph),pu*ph)); w=np.exp(np.array([-2.,.5,2.]))
    post=joint*w[:,None];post/=post.sum()
    tv=.5*np.abs(post-joint).sum()
    check('topology_independent_of_signature',abs(post[:,1].sum()-ph)<1e-14,
          material_posterior_TV=float(tv),hole_probability_before=ph,hole_probability_after=float(post[:,1].sum()))
    # Event depends on signature: exact conditional covariance identity and bound.
    mh=np.array([.05,.4,.95]); pH=float(pu@mh); Ew=float(pu@w)
    pH1=float(pu@(mh*w)/Ew); cov=float(pu@((mh-pH)*(w-Ew)))
    var_m=float(pu@(mh-pH)**2); var_w=float(pu@(w-Ew)**2)
    tau=var_m/(pH*(1-pH)); bound=np.sqrt(pH*(1-pH)*tau)*np.sqrt(var_w)/Ew
    check('topology_signature_coupling',abs((pH1-pH)-cov/Ew)<1e-14 and abs(pH1-pH)<=bound+1e-14,
          before=pH,after=pH1,tau=float(tau),absolute_change=abs(pH1-pH),coupling_bound=float(bound))
    # Sharp oscillation / TV bound achieved by a two-state construction.
    omega=1.6; q=np.array([1/(1+np.exp(omega/2)),1/(1+np.exp(-omega/2))])
    pp=normalize(np.log(q)+np.array([omega,0.])); tv=.5*np.abs(pp-q).sum(); bound=np.tanh(omega/4)
    check('sharp_tilt_TV',abs(tv-bound)<1e-14,TV=float(tv),bound=float(bound))
    # A constant field offset is not a constant potential offset.
    f=np.array([-1.,1.]); yy=.4; offset=2.
    old=normalize(-.5*(f-yy)**2); new=normalize(-.5*(f+offset-yy)**2)
    check('common_field_offset_changes_odds',abs(old[1]-new[1])>.1,
          J_defect=0,probability_before=float(old[1]),probability_after=float(new[1]))
    # Exact finite-field secant, generic complex VIE algebra only.
    n=10; G=(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))*.02
    e_inc=rng.normal(size=(n,2))+1j*rng.normal(size=(n,2))
    x0=.4+.1j+.2*rng.random(n); dx=.1*(rng.normal(size=n)+1j*rng.normal(size=n))
    L0=np.eye(n)-x0[:,None]*G; L1=np.eye(n)-(x0+dx)[:,None]*G
    j0=np.linalg.solve(L0,x0[:,None]*e_inc); j1=np.linalg.solve(L1,(x0+dx)[:,None]*e_inc)
    predicted=np.linalg.solve(L1,dx[:,None]*(e_inc+G@j0))
    err=np.linalg.norm(j1-j0-predicted)/max(1.,np.linalg.norm(j1-j0))
    check('nonlinear_finite_secant',err<1e-12,error=float(err),scope='Generic VIE matrix, not a physical dyadic Maxwell experiment')
    # Equal best residual, different Bayesian evidence (proper normalized priors).
    sigma=.1; A=2.; evidence_ratio=np.sqrt(2*np.pi)*sigma/(2*A)*(2*ndtr(A/sigma)-1)
    check('occam_equal_best_residual',evidence_ratio<.1,
          both_minimum_residuals=0,BF_flexible_vs_fixed=float(evidence_ratio),posterior_flexible_equal_model_prior=float(evidence_ratio/(1+evidence_ratio)))
    # Discrete covariance prior: likelihood sees only u, yet v marginal updates.
    prior_joint=np.array([[.45,.05],[.05,.45]])
    ll=np.array([.1,.9]); postj=prior_joint*ll[:,None];postj/=postj.sum()
    old_cond=prior_joint/prior_joint.sum(axis=1)[:,None]
    new_cond=postj/postj.sum(axis=1)[:,None]
    check('fiber_prior_conditional_not_marginal',np.max(np.abs(old_cond-new_cond))<1e-14,
          marginal_v1_prior=float(prior_joint[:,1].sum()),marginal_v1_posterior=float(postj[:,1].sum()))
    # Variable mobility without Ito correction: stationary density is pi/M.
    grid=np.linspace(-3,3,1201); target=normalize(-grid**2/2); mob=1+.5*grid**2
    wrong=target/mob;wrong/=wrong.sum()
    check('missing_divergence_bias',abs(target@grid**2-wrong@grid**2)>.1,
          intended_second_moment=float(target@grid**2),uncorrected_second_moment=float(wrong@grid**2))
    # Finite-support witness of TV >= target mass excluded by the approximation.
    target_mass=np.array([.3,.4,.1,.2]); approximate_mass=np.array([.5,.5,0.,0.])
    excluded=float(target_mass[approximate_mass==0].sum())
    tv=float(.5*np.abs(target_mass-approximate_mass).sum())
    check('support_exclusion_bound',tv+1e-14>=excluded and abs(excluded-.3)<1e-14,
          example_full_posterior_mass_outside_support=excluded,TV_lower_bound=excluded,actual_TV=tv,
          scope='Computed finite-support witness; not a physical measurement or proof of general discretization coverage')

def original_a30_reaggregation():
    # Manually transcribed PARENT_STATISTICS.csv rows from pinned original source.
    # This checks aggregation only; it does not rerun fits, physics, or read labels.
    gains={
    'full_hull':[.15039768691308428,.11065399695233184,.18331733822196034,.10294346715217474,.1388359359854376,.15056176596829468,.1345267865058784,.11501839317095551,.14269357340508015,.12634054040009918,.06925013401461053,.15016010828019474],
    'full_L_mean':[-.025124372506541787,-.03894838287368815,.011473174523050654,-.12644326923110155,-.04557907994635282,.0006974522572794427,-.04491750630737762,.03854445244385085,-.05140996789004735,-.02792355562920762,-.08246439778672388,-.005749041532960068],
    'sparse_64':[.40694577171182134,.16155566084281534,.3540659784556771,.6520125088274469,.29096731508192125,.352323182619962,.33598250389437423,.39709818239079125,.43170737241343693,.35229375010385106,.3382788290683077,.2095237983024259],
    'single_chart_64':[.4127697490748703,.2047036591960957,.2791075069452215,-.022478549432477398,.28581544251094004,.3348060931975827,.20354639920396572,.2185083981208151,.2844317885723488,.10456220902431769,.07137590596052683,.2517422161586392]}
    out={}
    for name,vals in gains.items():
        v=np.array(vals); rr=np.random.default_rng(291011)
        boot=v[rr.integers(0,12,size=(2000,12))].mean(axis=1)
        out[name]={'parent_mean_gain':float(v.mean()),'bootstrap_ci95':np.quantile(boot,[.025,.975]).tolist(),'parent_gains':vals}
    return out

def physical_demo():
    # Independent boundary-condition controls and Fourier truncation check.
    vacuum=cylinder(1.,.12,4.)
    solid_vacuum=cylinder(1.,0.,4.)
    s16=cylinder(2.5+.08j,.16,7.,nmax=16)
    s20=cylinder(2.5+.08j,.16,7.,nmax=20)
    trunc=np.linalg.norm(s16-s20)/max(np.linalg.norm(s20),1e-30)
    check('analytic_TM_controls',max(np.linalg.norm(vacuum),np.linalg.norm(solid_vacuum))<1e-10 and trunc<1e-10,
          vacuum_scatter_norm=float(np.linalg.norm(vacuum)),solid_vacuum_norm=float(np.linalg.norm(solid_vacuum)),relative_n16_n20_error=float(trunc))
    eps=2.5+.08j; k=4.; base=cylinder(eps,0.,k)
    sigma=.1*np.linalg.norm(base)/np.sqrt(8)  # fixed 10% RMS complex noise
    rows=[]
    for a in [.0035,.0105,.035,.105,.1925]:
        f=cylinder(eps,a,k);d=np.sqrt(2)*np.linalg.norm(f-base)/sigma
        rows.append({'a_over_R':a/.35,'whitened_distance':float(d),'optimal_equal_prior_error':float(ndtr(-d/2))})
    check('small_cavity_indistinguishability',rows[0]['optimal_equal_prior_error']>.49,
          radius_results=rows,scope='Fixed pairs, known permittivity, scalar 2-D TM, Gaussian known noise; not composite topology evidence')
    # Explicit finite PRIOR, not a quadrature masquerading as a continuous prior.
    eps_re=np.linspace(1.8,3.2,9); holes=[0.,.16]; freqs=[4.]
    responses=[]; labels=[]; params=[]
    for H,a in enumerate(holes):
        for er in eps_re:
            z=np.concatenate([cylinder(er+.08j,a,kk) for kk in freqs])
            responses.append(realify(z));labels.append(H);params.append([H,float(er),a])
    responses=np.array(responses); labels=np.array(labels); n=len(labels)
    # Sigma is experiment-fixed from the SAME reference, not candidate-dependent.
    noise_scale=.30*np.linalg.norm(responses[4])/np.sqrt(responses.shape[1])
    obs_dim=responses.shape[1]
    true_index=13; obs=responses[true_index]+noise_scale*rng.normal(size=obs_dim)
    loglik=-.5*np.sum((responses-obs)**2,axis=1)/noise_scale**2
    w=normalize(loglik); p_hole=float(np.clip(w@labels,0.,1.))
    # Proper common real multiplicative calibration gain, marginalized analytically.
    gain_sd=.15; r=obs[None,:]-responses
    ff=np.sum(responses**2,axis=1); rf=np.sum(r*responses,axis=1)
    qform=np.sum(r*r,axis=1)/noise_scale**2-gain_sd**2*rf**2/(noise_scale**4*(1+gain_sd**2*ff/noise_scale**2))
    logdet=np.log1p(gain_sd**2*ff/noise_scale**2)
    wg=normalize(-.5*(qform+logdet))
    # Prior-predictive calibration test within the declared 18-state model.
    trials=2000; idx=rng.integers(0,n,size=trials)
    yy=responses[idx]+noise_scale*rng.normal(size=(trials,obs_dim))
    L=-.5*np.sum((yy[:,None,:]-responses[None,:,:])**2,axis=2)/noise_scale**2
    weights=np.exp(L-logsumexp(L,axis=1)[:,None]); probs=np.clip(weights@labels,0.,1.); truth=labels[idx]
    brier=float(np.mean((probs-truth)**2)); bins=[]
    for l,h in zip(np.linspace(0,1,6)[:-1],np.linspace(0,1,6)[1:]):
        sel=(probs>=l)&(probs<=h if h==1 else probs<h)
        if sel.any():
            bins.append({'range':[float(l),float(h)],'count':int(sel.sum()),'mean_probability':float(probs[sel].mean()),'observed_frequency':float(truth[sel].mean()),'SE_under_independence':float(np.sqrt(np.mean(probs[sel]*(1-probs[sel]))/sel.sum()))})
    check('finite_topology_probability_normalization',abs(w.sum()-1)<1e-13 and 0<=p_hole<=1 and sum(b['count'] for b in bins)==trials,
          p_hole_known_calibration=p_hole,p_hole_gain_marginalized=float(np.clip(wg@labels,0.,1.)),finite_prior_states=n,
          true_state=params[true_index],noise_real_sd=float(noise_scale),Brier_prior_predictive=brier,calibration_trials=trials,
          calibration_bins=bins,scope='Exact finite-prior model plus finite Monte Carlo audit; not full continuous-material posterior calibration')
    return {'finite_prior_parameters':params,'one_observation_weights':w.tolist(),'one_observation_gain_marginalized_weights':wg.tolist(),
            'boundary_backward_residual_max':max_boundary_residual}

def main():
    t=time.perf_counter();algebra();a30=original_a30_reaggregation();physical=physical_demo()
    out={'schema':'A32-v1','seed':SEED,'status':'PASS' if all(x['passed'] for x in checks.values()) else 'FAIL',
         'check_count':len(checks),'checks':checks,'a30_parent_reaggregation':a30,'analytic_TM_details':physical,
         'resources':{'analytic_material_frequency_evaluations':physical_calls,'analytic_call_cap':physical_limit,'cumulative_analytic_calls':resource_ledger['charged_calls'],
                      'generic_matrix_examples':'small CPU only','training_updates':0,'GPU_calls':0,'locked_arrays_read':0,
                      'historical_truth_arrays_read':0,'elapsed_seconds':time.perf_counter()-t},
         'environment':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},
         'limitations':['Floating-point checks, not formal verified numerics.','Analytic cylinder is 2-D scalar TM, not 3-D vector Maxwell.',
                        'A30 aggregation checks transcribed parent scores, not a replay of the original optimizers.',
                        'Initial two-frequency example was too easy for calibration; this one-frequency 30-percent-noise check is exploratory, not independent validation.',
                        'Finite-prior posterior correctness does not establish continuous-prior coverage, algorithmic superiority, or acceleration.']}
    Path(__file__).with_name('A32_CHECKS.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':out['status'],'checks':len(checks),'resources':out['resources'],
                      'failed':[k for k,v in checks.items() if not v['passed']],
                      'a30_means':{k:v['parent_mean_gain'] for k,v in a30.items()}},indent=2))
    if out['status']!='PASS':raise SystemExit(1)
if __name__=='__main__':main()
