"""Pure full66 coordinate/Jacobian replay; no optimizer or Maxwell calls."""
import json
from pathlib import Path
import resource
import sys
import time

import numpy as np
import torch

from p3r4.conditional_coordinates import ConditionalCoordinateMap, fixed_random_rotation, native_to_gaussianized_jacobian, covariance_cholesky_from_factor
from p3r4.anchor_density import AnchorDensityPredictor, context_features, FORMAT
from p3r4.train_density import fit_normalization
from p3r4.density import make_model

wall,cpu=time.monotonic(),time.process_time()
torch.set_num_threads(1)
rng=np.random.default_rng(59)
n=8
cx=rng.normal(size=(n,66));J=rng.normal(size=(n,512,66))*.05
y=rng.normal(size=(n,512));x=cx+rng.normal(size=(n,66))*.15
R=fixed_random_rotation(20261025113)
results={}
for mode in ('gn','random'):
 coordinate=ConditionalCoordinateMap(mode,cx,J,R if mode=='random' else None)
 Jx=native_to_gaussianized_jacobian(cx,J)
 H=np.eye(66)+Jx.transpose(0,2,1)@Jx
 error=np.max(abs(coordinate.upper.transpose(0,2,1)@coordinate.upper-H))
 assert error<1e-12,error
 delta=x-cx;r=coordinate.forward(delta);back=coordinate.inverse(r)
 err64=float(np.max(abs(delta-back)));assert err64<1e-12
 # Quantize actual stored upper and transformed residual to float32; inverse
 # computation in this audit is genuinely float32 (production usesfloat64 CPU).
 upper32=torch.tensor(coordinate.upper,dtype=torch.float32)
 rot32=torch.tensor(R,dtype=torch.float32)
 d32=torch.tensor(delta,dtype=torch.float32)
 raw32=(d32@rot32 if mode=='random' else d32)
 r32=(upper32@raw32[:,:,None]).squeeze(-1)
 back32=torch.linalg.solve_triangular(upper32,r32[:,:,None],upper=True).squeeze(-1)
 if mode=='random':back32=back32@rot32.T
 err32=float(torch.max(abs(back32-d32)));assert err32<2e-6
 mat=torch.tensor(coordinate.matrix()[0],dtype=torch.float64)
 z0=torch.tensor(delta[0],dtype=torch.float64,requires_grad=True)
 jac=torch.autograd.functional.jacobian(lambda z:mat@z,z0)
 ld=float(torch.linalg.slogdet(jac).logabsdet)
 ld_err=abs(ld-coordinate.logabsdet[0]);assert ld_err<1e-12
 # Random rawbase covariance is a conjugate of the GN covariance; all66
 # eigenvalues remain represented. Compute by triangular solves on basis.
 basis=np.broadcast_to(np.eye(66),(n,66,66))
 raw_factor_rows=coordinate.inverse(basis)
 cov=raw_factor_rows.transpose(0,2,1)@raw_factor_rows
 if mode=='gn':gn_cov=cov.copy()
 else:
  assert np.max(abs(cov-R@gn_cov@R.T))<1e-12
  assert np.max(abs(np.linalg.eigvalsh(cov)-np.linalg.eigvalsh(gn_cov)))<1e-12
 results[mode]=dict(gram_factor_max=float(error),inverse_float64_max=err64,
                    inverse_float32_max=err32,autodiff_full66_logdet_error=ld_err,
                    min_covariance_eigenvalue=float(np.linalg.eigvalsh(cov).min()))
 feature=context_features(y,cx)
 stats=fit_normalization(r,feature)
 for kind in ('gaussian','diag_gaussian','spline','triangular','cfm'):
  config=dict(dimension=66,context_dim=578,width=16)
  if kind!='cfm':config.update(layers=4,output_rank=8,diag_hidden=32)
  model=make_model(kind,**config)
  for base_space in (('coordinate','raw') if kind=='cfm' else ('coordinate',)):
   ac=dict(j_summary='none',residual_center=True,anchor_scale=.3,coordinate_preconditioner=mode,
           alignment_rotation=R.tolist() if mode=='random' else None,preconditioner_base=base_space)
   payload=dict(format=FORMAT,kind=kind,model_config=config,state_dict=model.state_dict(),normalization=stats,
                contract=dict(cfm_base='anchor' if kind=='cfm' else 'standard'),anchor_contract=ac)
   predictor=AnchorDensityPredictor(payload)
   kwargs=dict(context_J=J)
   if kind=='cfm':
    lp=predictor.log_prob(x,y,cx,allow_approximate=True,nfe=4,**kwargs)
    v=r/.3 if base_space=='coordinate' else (x-cx)/.3
    expected=-.5*(v*v+np.log(2*np.pi)).sum(1)-66*np.log(.3)
    if base_space=='coordinate':expected+=coordinate.logabsdet
    samples=predictor.sample(y,cx,n_samples=3,nfe=4,batch_size=4,**kwargs)
   else:
    lp=predictor.log_prob(x,y,cx,**kwargs)
    z=(r-stats['x_mean'])/stats['x_scale']
    expected=-.5*(z*z+np.log(2*np.pi)).sum(1)-np.log(stats['x_scale']).sum()+coordinate.logabsdet
    samples=predictor.sample(y,cx,n_samples=3,batch_size=4,**kwargs)
   err=float(np.max(abs(lp-expected)));assert err<3e-4,(mode,kind,base_space,err)
   assert samples.shape==(n,3,66) and np.isfinite(samples).all()
   results[f'{mode}_{kind}_{base_space}']=dict(raw_logdensity_reference_max=err,all66_samples_finite=True)

# Rank-deficient data Jacobian is legitimate; identity prior retains all66 axes.
Jrank=np.zeros((1,512,66));Jrank[0,:,0]=1/np.sqrt(512);Jrank[0,:,1]=1e8/np.sqrt(512)
rankmap=ConditionalCoordinateMap('gn',cx[:1],Jrank)
a=rng.normal(size=(1,66));r=rankmap.forward(a)
back=rankmap.inverse(r)
# Explicitly preserve float32 quantization failure evidence at extreme scaling.
u32=torch.tensor(rankmap.upper,dtype=torch.float32);a32=torch.tensor(a,dtype=torch.float32)
r32=(u32@a32[:,:,None]).squeeze(-1)
b32=torch.linalg.solve_triangular(u32,r32[:,:,None],upper=True).squeeze(-1)
results['ill_conditioned_rank1_stress']=dict(J_scale=1e8,
    float64_roundtrip_max=float(np.max(abs(back-a))),
    float32_roundtrip_max=float(torch.max(abs(b32-a32))),
    minimum_upper_diagonal=float(np.diagonal(rankmap.upper,axis1=1,axis2=2).min()),
    interpretation='finite fullrank coordinate map; float32 accuracy can degrade under extreme mixed scaling; no H floor/rank truncation introduced')
assert np.isfinite(back).all() and results['ill_conditioned_rank1_stress']['minimum_upper_diagonal']>0

# Stable covariance factor QR retains tiny positive singular values without
# forming A A.T, including a strongly mixed factor.
a=np.array([[[1.,1e8],[0.,1.]]])
chol=covariance_cholesky_from_factor(a)
assert np.diagonal(chol,axis1=1,axis2=2).min()>0
results['stable_covariance_factor_positive']=True
results.update(status='PASS_WITH_RECORDED_FLOAT32_STRESS_LIMIT',physics_calls=0,optimizer_updates=0,
               wall_seconds=time.monotonic()-wall,cpu_seconds=time.process_time()-cpu,
               peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024))
(Path(__file__).parent/'NUMERICAL_COORDINATES.json').write_text(json.dumps(results,indent=2)+'\n')
# Self-contained zero-update data interface fixtures, no physical simulator.
for split,count in [('train',8),('validation',4)]:
 np.savez_compressed(Path(__file__).parent/f'synthetic_{split}.npz',x=x[:count],y=y[:count],context_x=cx[:count],context_J=J[:count],base_mean=cx[:count],parent=np.array([f'{split}_{i}' for i in range(count)]))
print(json.dumps(results,indent=2))
