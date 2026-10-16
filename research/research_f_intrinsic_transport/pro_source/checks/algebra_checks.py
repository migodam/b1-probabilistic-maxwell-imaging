"""Research F: small algebraic checks only.
No Maxwell mesh, physical data generation, training, or historical replay.
Run: python algebra_checks.py
Requires numpy, scipy. Outputs algebra_checks_results.json.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp

def rel(a, b):
    return float(np.linalg.norm(a-b)/max(np.linalg.norm(b),1e-15))

def main():
    rng=np.random.default_rng(20261016)
    n,p=6,3
    G=(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))*0.07
    np.fill_diagonal(G,0)
    chi=0.3+rng.random(n)*0.2+0.08j
    h=rng.normal(size=n)*0.1+1j*rng.normal(size=n)*0.01
    X=np.diag(chi); H=np.diag(h); I=np.eye(n)
    e=rng.normal(size=(n,p))+1j*rng.normal(size=(n,p))
    S=rng.normal(size=(4,n))+1j*rng.normal(size=(4,n))
    L=I-X@G; R=np.linalg.inv(L)
    j=R@X@e; E=e+G@j
    Ln=I-(X+H)@G
    jn=np.linalg.solve(Ln,(X+H)@e)
    b=R@H@E; K=R@H@G
    dj=np.linalg.solve(I-K,b)
    out={"scope":"ALGEBRA_ONLY_NOT_MAXWELL_SIMULATION_OR_HISTORICAL_REPLAY"}
    out["finite_resolvent_relative_error"]=rel(dj,jn-j)
    orders={}
    pred=np.zeros_like(b); term=b.copy()
    for k in range(1,7):
        pred+=term
        rem=np.linalg.solve(I-K,np.linalg.matrix_power(K,k)@b)
        orders[str(k)]={"prediction_relative_error":rel(pred,dj),
                       "remainder_identity_absolute_error":float(np.linalg.norm(dj-pred-rem))}
        term=K@term
    out["neumann_orders"]=orders
    # Exact derivatives and numerical central differences.
    D1=R@H@E
    D2=2*R@H@G@D1
    D3=3*R@H@G@D2
    def jof(t):
        Xt=X+t*H
        return np.linalg.solve(I-Xt@G,Xt@e)
    dt=1e-3
    fd1=(jof(dt)-jof(-dt))/(2*dt)
    fd2=(jof(dt)-2*jof(0)+jof(-dt))/dt**2
    fd3=(jof(2*dt)-2*jof(dt)+2*jof(-dt)-jof(-2*dt))/(2*dt**3)
    out["derivative_FD_relative_errors"]={"first":rel(fd1,D1),"second":rel(fd2,D2),"third":rel(fd3,D3)}
    # Field propagator along a non-straight material path.
    h2=rng.normal(size=n)*0.03
    def Xt(t): return np.diag(chi+t*h+np.sin(2*np.pi*t)*h2)
    def Xdot(t): return np.diag(h+2*np.pi*np.cos(2*np.pi*t)*h2)
    A0=I-G@Xt(0)
    Uexact=np.linalg.solve(I-G@Xt(1),A0)
    def ode(t, flat):
        U=flat.reshape(n,n)
        return (np.linalg.solve(I-G@Xt(t),G@Xdot(t))@U).ravel()
    ivp=solve_ivp(ode,(0,1),I.astype(complex).ravel(),rtol=2e-11,atol=2e-13)
    out["ordered_field_propagator_endpoint_identity_error"]=rel(ivp.y[:,-1].reshape(n,n),Uexact)
    # Multi-source exact material realizability from current changes.
    a=L@dj; bfield=E+G@dj
    hrecovered=np.sum(a*np.conj(bfield),axis=1)/np.sum(abs(bfield)**2,axis=1)
    out["current_compatibility_recovered_material_error"]=rel(hrecovered,h)
    out["multi_source_collinearity_residual"]=float(np.linalg.norm(a-hrecovered[:,None]*bfield))
    rho=rng.normal(size=dj.shape)+1j*rng.normal(size=dj.shape)
    khat=dj+0.002*rho
    residual=H@E-Ln@khat
    Z=np.linalg.solve(Ln.conj().T,S.conj().T)
    out["primal_dual_finite_output_identity_error"]=rel(Z.conj().T@residual,S@(dj-khat))
    # Nonnormal norm>1 but convergent Neumann series.
    Kn=np.array([[.2,5.],[0,.3]])
    bn=np.array([0.,1.])
    exact=np.linalg.solve(np.eye(2)-Kn,bn)
    summ=np.zeros(2); term=bn.copy()
    for _ in range(40):
        summ+=term; term=Kn@term
    out["nonnormal_neumann"]={"norm":float(np.linalg.norm(Kn,2)),
        "spectral_radius":float(max(abs(np.linalg.eigvals(Kn)))),
        "forty_term_relative_error":rel(summ,exact)}
    # Invertible resolvent despite divergent Neumann expansion.
    out["divergent_series_invertible_resolvent"]={"K":2.,"b":1.,"exact":-1.,
        "twenty_terms":float(sum(2.**k for k in range(20)))}
    # Maxwell-form two-cell passive reciprocal example, NOT a continuum experiment.
    coupling=.2+.05j; a0=.7+.03j; q=.05
    G2=np.array([[0,coupling],[coupling,0]],complex)
    def twocell(q):
        X2=np.diag([a0+q,a0-q])
        return np.linalg.solve(np.eye(2)-X2@G2,X2@np.array([1.,0.]))[1]
    fplus=twocell(q); fminus=twocell(-q)
    out["two_cell_disconnected_prior_fiber"]={
        "equal_data_absolute_error":float(abs(fplus-fminus)),
        "material_pair_L2_distance":float(2*np.sqrt(2)*q),
        "loss_imaginary_part":float(a0.imag),
        "midpoint_data_mismatch":float(abs(twocell(0)-fplus)),
        "analytic_formula_error":float(abs(fplus-coupling*(a0*a0-q*q)/(1-coupling*coupling*(a0*a0-q*q)))),
        "claim_scope":"Fixed-total-contrast 1D prior; q and -q isolated iso-data points away from q=0."}
    # Coordinate-artifact example: entire material plane is flat.
    # g(u,v)=(u+c*v^2,v); F(x,y)=x.
    c=1000.; t=.01
    raw=np.array([0+c*t*t,t])
    geodesic=np.array([-c*t*t+c*t*t,t])
    out["decoder_curvature_artifact"]={"raw_decoder_second_derivative":2*c,
        "intrinsic_prior_second_fundamental_form":0.,
        "raw_coordinate_step_data_drift":float(raw[0]),
        "geodesic_coordinate_step_data_drift":float(geodesic[0]),
        "explanation":"Raw decoder curvature may be arbitrarily large while the prior manifold is flat."}
    # Gap-free soft weak transport inequality.
    A=rng.normal(size=(4,7))
    lam=.04
    P=lam*np.linalg.inv(A.T@A+lam*np.eye(7))
    out["soft_weak_filter"]={"operator_norm_AP":float(np.linalg.norm(A@P,2)),
                            "universal_bound_sqrt_lambda_over_2":float(np.sqrt(lam)/2)}
    # Sharp discrete Bayes risk identity for indistinguishable endpoints.
    truths=np.array([-1.,1.])
    meanrisk=float(np.mean(truths**2))
    samplerisk=float(np.mean((truths[:,None]-truths[None,:])**2))
    out["bayes_risk_example"]={"posterior_mean_MSE":meanrisk,"independent_posterior_sample_MSE":samplerisk}
    assert out["finite_resolvent_relative_error"]<1e-11
    assert out["current_compatibility_recovered_material_error"]<1e-11
    assert out["ordered_field_propagator_endpoint_identity_error"]<1e-8
    assert out["primal_dual_finite_output_identity_error"]<1e-10
    assert out["soft_weak_filter"]["operator_norm_AP"]<=np.sqrt(lam)/2+1e-12
    assert out["two_cell_disconnected_prior_fiber"]["equal_data_absolute_error"]<1e-12
    output=Path(__file__).with_name("algebra_checks_results.json")
    output.write_text(json.dumps(out,indent=2),encoding="utf-8")
    print(json.dumps(out,indent=2))

if __name__=="__main__":
    main()
