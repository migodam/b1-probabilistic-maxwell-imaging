"""Analytic/synthetic saved-array controls, no physical matrices or RHS."""
from pathlib import Path
import json,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[4];sys.path.insert(0,str(ROOT/'src'))
from p3l.decision_bounds import stability_from_saved_inverse,dual_corrected_prediction,same_uniform_decision
OUT=Path(__file__).resolve().parent
rng=np.random.default_rng(2026102251);checks={}
def check(name,ok,**values):
    checks[name]=dict(passed=bool(ok),**values)
    if not ok:raise AssertionError(name)
def cx(shape):return rng.normal(size=shape)+1j*rng.normal(size=shape)
L0=np.diag([2.,3.]).astype(complex);R0=np.diag([.5,1/3]);Lnew=np.diag([2.1,2.9])
s=stability_from_saved_inverse(L0,R0,Lnew)
check('analytic_diagonal_beta',0<s['beta']<=2.1 and s['kappa']<1,beta=s['beta'],kappa=s['kappa'],rho0=s['rho0_upper'])
bad=stability_from_saved_inverse(np.eye(2),np.eye(2),3*np.eye(2))
check('stable_matrix_outside_sufficient_region',bad['status']=='BOUND_UNAVAILABLE' and bad['beta'] is None and bad['kappa']>=1,true_sigma_min=3.)
bad_inverse=stability_from_saved_inverse(np.eye(2),3*np.eye(2),np.eye(2))
check('bad_saved_inverse_not_singular_matrix',bad_inverse['status']=='BOUND_UNAVAILABLE' and bad_inverse['rho0_upper']>=1)
# Small synthetic complex nonnormal matrices: arbitrary projected state, saved
# approximate dual. Full tiny solve is test-only and never passed to evaluator.
trials=[];decision_cases=0;tight_cases=0
for case in range(20):
    n=6;p=3;m=2
    L0=np.diag(np.arange(3,9))+0.03*cx((n,n))
    R0=np.linalg.inv(L0)
    L=L0+0.01*cx((n,n))
    stability=stability_from_saved_inverse(L0,R0,L)
    rhs=cx((n,m));S=cx((p,n));V=np.linalg.qr(cx((n,2)))[0]
    Z0=np.linalg.solve(L0.conj().T,S.conj().T)+.01*cx((n,p))
    y=cx((p,m));sigma=.2
    answer=dual_corrected_prediction(L,rhs,S,V,Z0,stability,observation=y,sigma=sigma)
    truth=S@np.linalg.solve(L,rhs)
    error=float(np.linalg.norm(truth-answer['prediction']))
    phi=float(np.sum(abs(truth-y)**2)/sigma**2)
    lo,hi=answer['likelihood_interval']
    check(f'enclosure_{case}',error<=answer['eps_complex'] and lo<=phi<=hi and stability['beta']<=np.linalg.svd(L,compute_uv=False)[-1],output_error=error,eps_complex=answer['eps_complex'],true_phi=phi,interval=[lo,hi],beta=stability['beta'])
    # Construct a paired state interval without giving true values to chooser.
    current=(phi+.1,phi+.1);candidate=(lo,hi);c=-.2
    for logu in [-.001,-.07,-1.,-10.,-100.]:
        decision=same_uniform_decision(current,candidate,c,logu)
        exact='ACCEPT' if logu<c+(phi+.1)-phi else 'REJECT'
        if decision['action']!='REFINE':
            if decision['action']!=exact:raise AssertionError('decision mismatch')
            decision_cases+=1
        tightened=same_uniform_decision(current,(phi,phi),c,logu)
        if tightened['action']!=exact or tightened['logu']!=logu:raise AssertionError('tightened decision mismatch')
        tight_cases+=1
    trials.append(answer)
check('shared_uniform_decisions_match_when_resolved',decision_cases>0,resolved_cases=decision_cases)
check('same_uniform_tightened_decisions_match',tight_cases==100,resolved_cases=tight_cases)
# Exact old dual gives vanishing algebraic error but still explicit roundoff.
L=np.diag([2.,3.]);rhs=np.ones((2,1));S=np.array([[1.,1.]])
a=dual_corrected_prediction(L,rhs,S,np.array([[1.],[0.]]),np.array([[.5],[1/3]]),2.,observation=np.zeros((1,1)),sigma=1.)
exact=np.array([[5/6]])
check('dual_correction_exact_in_real_arithmetic',np.linalg.norm(a['prediction']-exact)<1e-15 and a['eps_roundoff_term']>0 and a['eps_complex']>0,eps_complex=a['eps_complex'],prediction=a['prediction'].real.tolist())
r=same_uniform_decision((1,1),(2,2),0.,-.5)
check('known_rejection',r['action']=='REJECT' and r['logu']==-.5)
a=same_uniform_decision((2,2),(1,1),0.,-.5)
check('known_acceptance',a['action']=='ACCEPT')
r=same_uniform_decision((0,2),(0,2),0.,-.5)
check('uncertain_requires_refinement',r['action']=='REFINE' and r['retain_same_uniform'])
r=same_uniform_decision(None,(0,1),0.,-.5)
check('unavailable_interval_requires_refinement',r['action']=='REFINE')
a=dual_corrected_prediction(np.eye(2),np.ones((2,1)),np.ones((1,2)),np.eye(2),np.ones((2,1)),None)
check('missing_beta_no_projection_work',a['status']=='BOUND_UNAVAILABLE' and a['counts']['matrix_products']==0)
check('no_physical_counts',all(t['counts']['new_physical_rhs']==0 and t['counts']['new_maxwell_solves']==0 for t in trials))
report=dict(status='PASS_SMALL_MATRIX_CONDITIONAL_BOUND_CONTROLS_ONLY',seed=2026102251,checks=checks,numerical_assumptions='float64 gamma_k model, not directed interval arithmetic',new_maxwell=0,new_physical_rhs=0,new_physical_factorizations=0)
(OUT/'results.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
