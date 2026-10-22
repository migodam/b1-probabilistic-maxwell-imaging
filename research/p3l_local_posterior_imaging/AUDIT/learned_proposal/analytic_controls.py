"""Synthetic Gaussian controls only; no campaign arrays or physics imports."""
import json
from pathlib import Path
import sys
import numpy as np
from scipy.special import logsumexp
from scipy.stats import multivariate_normal

ROOT=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(ROOT/'src'))
from p3l.learned_proposal import LearnedFullProposal, GaussianWeakDensity


def main():
    rng=np.random.default_rng(2026102201)
    d,k,p=4,2,2;s=d-k
    frame=np.linalg.qr(rng.normal(size=(d,d)))[0];vw,vs=frame[:,:k],frame[:,k:]
    q0=np.r_[rng.normal(size=d),[.3,.7]]
    J=rng.normal(size=(9,d+p));J[:,d:]=0
    r=rng.normal(size=9)
    weak_mean=rng.normal(size=k+p)
    base=rng.normal(size=(k+p,k+p));weak_cov=base@base.T+.5*np.eye(k+p)
    network=GaussianWeakDensity(weak_mean,weak_cov,rng)
    proposal=LearnedFullProposal(network,np.zeros(3),np.zeros((1,2,2)),vw,vs,q0,J,r)
    H=np.eye(s)+(J[:,:d]@vs).T@(J[:,:d]@vs)
    C=np.linalg.inv(H)
    constant=-C@(vs.T@q0[:d]+(J[:,:d]@vs).T@r)
    B=-C@(J[:,:d]@vs).T@(J[:,:d]@vw)
    # z=q0+Vs constant+(Vw+Vs B)a+Vs epsilon, eta unchanged.
    transform=np.zeros((d+p,k+p));transform[:d,:k]=vw+vs@B;transform[d:,k:]=np.eye(p)
    joint_mean=np.r_[q0[:d]+vs@constant,np.zeros(p)]+transform@weak_mean
    joint_cov=transform@weak_cov@transform.T
    joint_cov[:d,:d]+=vs@C@vs.T
    x=rng.normal(size=(150,d+p))
    reference=multivariate_normal(joint_mean,joint_cov)
    learned_error=float(np.max(np.abs(proposal.learned_log_prob(x)-reference.logpdf(x))))
    prior=multivariate_normal(np.zeros(d+p),np.eye(d+p))
    exact_mix=logsumexp(np.stack((np.log(.9)+reference.logpdf(x),np.log(.1)+prior.logpdf(x))),axis=0)
    mix_error=float(np.max(np.abs(proposal.log_prob(x)-exact_mix)))
    a,b,eta=proposal.to_chart(x)
    inverse_error=float(np.max(np.abs(proposal.from_chart(a,b,eta)-x)))
    assert learned_error<1e-10 and mix_error<1e-10 and inverse_error<1e-12
    draws=proposal.sample(60000,rng=rng)
    expected_mean=.9*joint_mean
    expected_cov=.9*(joint_cov+np.outer(joint_mean,joint_mean))+.1*np.eye(d+p)-np.outer(expected_mean,expected_mean)
    se=np.sqrt(np.diag(expected_cov)/len(draws))
    mean_max_z=float(np.max(np.abs(draws.mean(axis=0)-expected_mean)/se))
    covariance_relative_error=float(np.linalg.norm(np.cov(draws,rowvar=False)-expected_cov)/np.linalg.norm(expected_cov))
    assert mean_max_z<6 and covariance_relative_error<.05
    # Correlated strong noise propagated into pixel marginals, not diagonal C.
    A=rng.normal(size=(5,d));outer_a=np.array([.2,-.4]);u=np.array([.25,.8])
    fm,fs=proposal.conditional_field_moments(outer_a,u,A,offset=.2)
    zmean=q0[:d]+vw@outer_a+vs@proposal.strong_mean(outer_a,u)
    field_mean_error=float(np.max(np.abs(fm-(.2+A@zmean))))
    field_var_error=float(np.max(np.abs(fs*fs-np.diag(A@vs@C@vs.T@A.T))))
    assert field_mean_error<1e-12 and field_var_error<1e-12
    mm=proposal.conditional_material_moments(outer_a,u,A,2+.4j,offset=.2)
    strong=rng.multivariate_normal(proposal.strong_mean(outer_a,u),C,size=100000)
    field=.2+(q0[:d]+vw@outer_a+strong@vs.T)@A.T
    t=np.clip((field+.5),0,1);fraction=t**3*(10+t*(-15+6*t))
    chi=.1+.04j+(2+.4j)*fraction
    material_mean_error=float(np.max(np.abs(chi.mean(axis=0)-mm.mean)))
    material_second_error=float(np.max(np.abs(np.mean(np.abs(chi)**2,axis=0)-mm.second_abs_moment)))
    assert material_mean_error<.025 and material_second_error<.06
    # Actual torch protocol, including no-gradient execution and dtype bridge.
    import torch
    class Fixed(torch.nn.Module):
        def __init__(self):
            super().__init__();self.register_buffer('loc',torch.zeros(k+p,dtype=torch.float64))
        def sample(self,n,scalar,image):
            return torch.randn((1,n,k+p),dtype=self.loc.dtype)
        def log_prob(self,values,scalar,image):
            assert not torch.is_grad_enabled()
            return -.5*(values.square().sum(-1)+(k+p)*np.log(2*np.pi))
    fixed=Fixed().eval()
    torch_control=LearnedFullProposal(fixed,np.zeros(3),np.zeros((1,2,2)),vw,vs,q0,J,r)
    torch_result=torch_control.log_prob(torch_control.sample(20,rng=rng))
    assert torch_result.shape==(20,) and np.isfinite(torch_result).all()
    # Full 66D campaign layout is exercised with synthetic zero J, never data.
    full=LearnedFullProposal(GaussianWeakDensity(np.zeros(10),np.eye(10),rng),np.zeros(1154),np.zeros((4,24,24)),np.eye(64)[:,:8],np.eye(64)[:,8:],np.r_[np.zeros(64),[.5,.5]],np.zeros((512,66)),np.zeros(512))
    control_x=rng.normal(size=(8,66))
    full_prior_error=float(np.max(np.abs(full.log_prob(control_x)-(-.5*(np.sum(control_x**2,axis=-1)+66*np.log(2*np.pi))))))
    assert full_prior_error<1e-11
    result=dict(status='PASS_SYNTHETIC_NORMALIZATION_CHART_MOMENT_CONTROLS',new_physical_RHS=0,truth_access=False,training_runs=0,
        learned_logprob_max_error=learned_error,defensive_logprob_max_error=mix_error,inverse_chart_max_error=inverse_error,
        mixture_samples=60000,mixture_mean_max_standard_errors=mean_max_z,mixture_covariance_relative_error=covariance_relative_error,
        field_mean_max_error=field_mean_error,field_variance_max_error=field_var_error,material_mc_samples=100000,
        material_mean_max_mc_error=material_mean_error,material_second_abs_max_mc_error=material_second_error,
        torch_interface='PASS',full66_prior_max_error=full_prior_error,
        limitations='Synthetic proposal controls only; not Maxwell posterior, exact MH wrapper, learned quality or scientific validation.')
    output=Path(__file__).with_name('RESULTS.json');output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))


if __name__=='__main__':main()
