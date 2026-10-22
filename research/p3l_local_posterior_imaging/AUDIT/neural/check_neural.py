"""CPU random/analytic density controls only; no campaign inputs or training."""
from pathlib import Path
import sys,json
import numpy as np
import torch
from scipy.stats import multivariate_normal
from scipy.special import logsumexp
ROOT=Path(__file__).resolve().parents[4];sys.path.insert(0,str(ROOT/'src'))
from p3l.neural import build_model,MODEL_TYPES,parameter_count
OUT=Path(__file__).resolve().parent
SEED=2026102252
torch.manual_seed(SEED);torch.set_num_threads(1)
checks={};counts={}
def check(name,ok,**detail):
    checks[name]=dict(passed=bool(ok),**detail)
    if not ok:raise AssertionError(name)
scalar=torch.randn(3,1154,dtype=torch.float64)
image=torch.randn(3,4,24,24,dtype=torch.float64)
for kind in MODEL_TYPES:
    model=build_model(kind).double().cpu()
    sample=model.sample(7,scalar,image)
    lp=model.log_prob(sample,scalar,image)
    loss=-model.log_prob(sample[:,0],scalar,image).mean()
    loss.backward()
    gradients=[p.grad for p in model.parameters() if p.grad is not None]
    check(kind+'_shape_gradient',sample.shape==(3,7,10) and lp.shape==(3,7) and torch.isfinite(lp).all().item() and all(torch.isfinite(g).all().item() for g in gradients),parameters=parameter_count(model),sample_shape=list(sample.shape))
    counts[kind]=parameter_count(model)
check('cnn_mlp_full_parameter_match',abs(counts['cnn_full']-counts['mlp_full'])/counts['mlp_full']<.01,relative_difference=abs(counts['cnn_full']-counts['mlp_full'])/counts['mlp_full'])

# 2D fixed Gaussian/mixture density compared to SciPy and integrated directly.
a=np.linspace(-10,10,201);xx,yy=np.meshgrid(a,a,indexing='ij')
points=np.c_[xx.ravel(),yy.ravel()]
values=torch.tensor(points,dtype=torch.float64)[None,:,:]
s=torch.zeros(1,1154,dtype=torch.float64);im=torch.zeros(1,4,24,24,dtype=torch.float64)
def integrate(logp):
    density=np.exp(logp).reshape(len(a),len(a))
    return float(np.trapezoid(np.trapezoid(density,a,axis=1),a,axis=0))
for kind in ['mlp_diag','mlp_full','mlp_mdn']:
    model=build_model(kind,output_dim=2).double()
    with torch.no_grad():
        model.head.weight.zero_();model.head.bias.zero_()
        if kind=='mlp_mdn':
            # 4*(logit + mean2 + lowertriangle3)
            raw=model.head.bias.reshape(4,6)
            raw[:,0]=torch.tensor([-.5,0,.3,.1],dtype=torch.float64)
            raw[:,1:3]=torch.tensor([[-1.,0],[0,1],[1,0],[0,-1]],dtype=torch.float64)
        logp=model.log_prob(values,s,im).numpy()[0]
        if kind=='mlp_mdn':
            logits,means,chol=model.distribution_parameters(s,im)
            logits,means,chol=logits.numpy()[0],means.numpy()[0],chol.numpy()[0]
            ref=logsumexp(np.array([multivariate_normal.logpdf(points,mean=means[k],cov=chol[k]@chol[k].T)+logits[k]-logsumexp(logits) for k in range(4)]),axis=0)
        else:
            mean,chol=model.distribution_parameters(s,im)
            mean,chol=mean.numpy()[0],chol.numpy()[0]
            ref=multivariate_normal.logpdf(points,mean=mean,cov=chol@chol.T)
    integral=integrate(logp)
    error=float(np.max(abs(logp-ref)))
    check(kind+'_density_normalized',abs(integral-1)<1e-9 and error<1e-10,grid_integral=integral,max_scipy_logdensity_error=error)

# Nontrivial random10D flow: inverse+logJacobian checked independently with
# a full automatic-differentiation Jacobian determinant.
model=build_model('mlp_flow').double()
z=torch.randn(3,7,10,dtype=torch.float64)
y,ld=model.forward_from_base(z,scalar,image)
back,invld=model.inverse_to_base(y,scalar,image)
check('flow_random_inverse',torch.max(abs(back-z)).item()<1e-12 and torch.max(abs(ld+invld)).item()<1e-12,inverse_error=torch.max(abs(back-z)).item(),logdet_cancellation=torch.max(abs(ld+invld)).item())
x0=z[0,0].detach().requires_grad_(True)
def mapping(x):return model.forward_from_base(x[None,:],scalar[:1],image[:1])[0][0]
J=torch.autograd.functional.jacobian(mapping,x0)
_,direct_ld=torch.linalg.slogdet(J)
_,predicted_ld=model.forward_from_base(x0[None,:],scalar[:1],image[:1])
check('flow_full_autograd_jacobian',abs(direct_ld.item()-predicted_ld.item())<1e-11,logdet_error=abs(direct_ld.item()-predicted_ld.item()))
logp=model.log_prob(y,scalar,image)
logbase=-.5*(z.square()+np.log(2*np.pi)).sum(-1)
check('flow_change_of_variables',torch.max(abs(logp+ld-logbase)).item()<1e-11,maximum_error=torch.max(abs(logp+ld-logbase)).item())

# A2D nonidentity affine flow has an independent analytic normal density and
# direct finite-grid normalization, not merely an inverse code self-check.
model=build_model('mlp_flow',output_dim=2).double()
with torch.no_grad():
    for k,layer in enumerate(model.layers):
        layer.net[-1].weight.zero_()
        layer.net[-1].bias.copy_(torch.tensor([.1*(k+1),-.07*(k+1),.06,-.04],dtype=torch.float64))
    center,_=model.forward_from_base(torch.zeros(1,2,dtype=torch.float64),s,im)
    e1,_=model.forward_from_base(torch.tensor([[1.,0.]],dtype=torch.float64),s,im)
    e2,_=model.forward_from_base(torch.tensor([[0.,1.]],dtype=torch.float64),s,im)
    affine=torch.stack([e1[0]-center[0],e2[0]-center[0]],dim=-1).numpy()
    logp=model.log_prob(values,s,im).numpy()[0]
    ref=multivariate_normal.logpdf(points,mean=center.numpy()[0],cov=affine@affine.T)
integral=integrate(logp);err=float(np.max(abs(logp-ref)))
check('flow_affine_normalization',abs(integral-1)<1e-9 and err<1e-10,grid_integral=integral,max_scipy_logdensity_error=err)
# Both lower-triangle head families preserve the .02 floor under large
# negative diagonal logits, with no clipping of labels or parameters.
model=build_model('mlp_full',output_dim=2).double()
with torch.no_grad():
    model.head.weight.zero_();model.head.bias.zero_();model.head.bias[2]=-100;model.head.bias[4]=-100
    _,chol=model.distribution_parameters(s,im)
check('cholesky_floor',torch.min(torch.diagonal(chol,dim1=-2,dim2=-1)).item()>=.02)
report=dict(status='PASS_CPU_RANDOM_ANALYTIC_DENSITY_CONTROLS_ONLY',seed=SEED,parameter_counts=counts,checks=checks,training_runs=0,campaign_data_reads=0,physical_calls=0,device='cpu',torch_version=torch.__version__)
(OUT/'results.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
