"""Pure numerical controls: no Target/Book, Maxwell, data, GT, or network."""
from pathlib import Path
import ast
import json
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(ROOT/'src'))
from p3l.spine import PolynomialSpine, SpinePosterior
from p3l.hmc import PreconditionedHMC
from p3l.kernels import ChainState


def difference(fun,x,h=1e-6):
    E=np.eye(len(x))*h
    return np.asarray([(fun(x+e)-fun(x-e))/(2*h) for e in E]).T


class AnalyticFullPosterior:
    """Synthetic nonquadratic likelihood plus the actual x Gaussian prior."""
    probit=True
    def __init__(self,A,b):self.A,self.b=A,b
    def residual(self,x):return self.A@x+self.b*np.sin(x[:len(self.b)])
    def __call__(self,x):
        u=self.residual(x);return -.5*float(u@u+x@x)
    def gradient(self,x):
        u=self.residual(x);g=-self.A.T@u-x
        g[:len(u)]-=self.b*np.cos(x[:len(u)])*u
        return g


def leapfrog_flip(q,p,C,h,n,gradient):
    q=q.copy();p=p.copy();p+=.5*h*gradient(q)
    for i in range(n):
        q+=h*(C@p);p+=(.5 if i==n-1 else 1)*h*gradient(q)
    return q,-p


def main():
    rng=np.random.default_rng(2026102297)
    coeff=rng.normal(size=(64,6))*.06
    spine=PolynomialSpine(coeff,np.array([.3,-.2]))
    # Independent polynomial formula and derivative include physical scale.
    r=rng.normal(size=66)*.25;r[64:]=np.array([.8,-.7])
    t=(r[64:]-spine.center)/.5
    beta=coeff[:,0]+coeff[:,1]*t[0]+coeff[:,2]*t[1]+coeff[:,3]*t[0]**2+coeff[:,4]*t[0]*t[1]+coeff[:,5]*t[1]**2
    beta_error=float(np.max(np.abs(beta-spine.beta(r[64:]))))
    x=spine.transform(r)
    inverse_error=float(np.max(np.abs(spine.inverse(x)-r)))
    DT=spine.jacobian(r);DTfd=difference(spine.transform,r)
    jacobian_error=float(np.max(np.abs(DT-DTfd)))
    # Independent LU determinant; avoids this environment's NumPy slogdet
    # warnings on a finite unit-triangular matrix.
    from scipy.linalg import lu
    P,L,U=lu(DT)
    assert np.array_equal(P,np.eye(66))
    diagonal=np.diag(U);sign=np.prod(np.sign(diagonal));logdet=np.log(np.abs(diagonal)).sum()
    prior_error=float(np.max(np.abs(spine.prior_gradient(r)-difference(spine.log_prior,r))))
    prior_vs_wrong=float(abs(spine.log_prior(r)+.5*r@r))
    A=rng.normal(size=(11,66))*.03;b=rng.normal(size=11)*.1
    base=AnalyticFullPosterior(A,b);post=SpinePosterior(base,spine)
    full_error=float(np.max(np.abs(post.gradient(r)-difference(post,r))))
    density_error=float(abs(post(r)-base(spine.transform(r))))
    # Explicit original Gaussian normalization under this exact global shear.
    normalization_identity_error=float(abs(post.spine.log_prior(r)+.5*x@x))
    p=rng.normal(size=66);D=rng.normal(size=(66,66))*.01
    C=D@D.T+np.eye(66)*.5
    q1,p1=leapfrog_flip(r,p,C,.015,12,post.gradient)
    q2,p2=leapfrog_flip(q1,p1,C,.015,12,post.gradient)
    reversal_q=float(np.max(np.abs(q2-r)));reversal_p=float(np.max(np.abs(p2-p)))
    # Exact kernel endpoint against an independent integrator, preserving the
    # RNG momentum and MH draw sequence. No posterior convergence claim.
    kernel=PreconditionedHMC(post,post.gradient,.015,C,12)
    seed=2026102298;krng=np.random.default_rng(seed);oracle=np.random.default_rng(seed)
    from scipy.linalg import solve_triangular
    p0=solve_triangular(np.linalg.cholesky(C).T,oracle.normal(size=66),lower=False)
    expected_q,expected_p=leapfrog_flip(r,p0,C,.015,12,post.gradient)
    ratio=post(expected_q)-post(r)+.5*p0@C@p0-.5*expected_p@C@expected_p
    accepted=np.log(max(oracle.random(),np.finfo(float).tiny))<min(0.,ratio)
    out=kernel.step(ChainState(r,post(r),post.gradient(r)),krng)
    kernel_error=float(np.max(np.abs(out.state.q-(expected_q if accepted else r))))
    ratio_error=float(abs(out.metadata['log_mh_ratio']-ratio))
    # Every r direction maps to a distinct x; arbitrarily large eta remains
    # defined and no local-grid clipping or rank reduction enters T.
    remote=r.copy();remote[64:]=[20.,-17.]
    remote_error=float(np.max(np.abs(spine.transform(spine.inverse(spine.transform(remote)))-spine.transform(remote))))
    runner=ROOT/'tools/run_p3l_spine.py'
    ast.parse(runner.read_text())
    checks=dict(beta_formula_error=beta_error,inverse_error=inverse_error,
        transform_jacobian_finite_difference_error=jacobian_error,
        jacobian_sign=float(sign),jacobian_log_abs_det=float(logdet),
        original_normal_prior_gradient_error=prior_error,
        wrong_standard_r_prior_log_density_difference=prior_vs_wrong,
        original_prior_normalization_identity_error=normalization_identity_error,
        full_target_pullback_gradient_error=full_error,full_density_composition_error=density_error,
        leapfrog_reversal_q_error=reversal_q,leapfrog_reversal_p_error=reversal_p,
        existing_HMC_endpoint_error=kernel_error,existing_HMC_energy_ratio_error=ratio_error,
        outside_grid_polynomial_roundtrip_error=remote_error,
        analytic_only=True,Maxwell_executed=False,GT_used=False,training=False,
        physics_jobs_executed=False,scientific_mixing_claim=False)
    assert beta_error<1e-12 and inverse_error<1e-12
    assert jacobian_error<1e-8 and sign==1 and abs(logdet)<1e-12
    assert prior_error<1e-8 and full_error<1e-8 and prior_vs_wrong>1e-3
    assert density_error==0 and normalization_identity_error==0
    assert reversal_q<1e-12 and reversal_p<1e-12
    assert kernel_error<1e-12 and ratio_error<1e-12 and remote_error<1e-10
    checks['status']='PASS_PURE_COORDINATE_AND_KERNEL_CONTROLS_ONLY'
    path=Path(__file__).with_name('results.json')
    path.write_text(json.dumps(checks,indent=2,allow_nan=False)+'\n')
    print(json.dumps(checks,indent=2,allow_nan=False))


if __name__=='__main__':main()
