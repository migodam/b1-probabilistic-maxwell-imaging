"""Zero-optimizer/zero-physics algebra replay of explicit anchor densities."""
import json
from pathlib import Path
import resource
import sys
import time

import numpy as np
import torch

from p3r4.anchor_density import FORMAT, AnchorDensityPredictor, context_features, base_parameters
from p3r4.density import make_model
from p3r4.train_density import fit_normalization
from p3r4.train_anchor_density import validate_metric, cfm_error_loss

wall,cpu=time.monotonic(),time.process_time()
torch.set_num_threads(1)
rng=np.random.default_rng(12)
n=8
y=rng.normal(size=(n,512));cx=rng.normal(size=(n,66));bmean=cx+.2
x=bmean+rng.normal(size=(n,66))*.1
J=rng.normal(size=(n,512,66))
frame=np.linalg.qr(rng.normal(size=(66,66)))[0]
results={}
for js in ('none','log_column_norm'):
 feature=context_features(y,cx,J,js)
 assert feature.shape==(n,578 if js=='none' else 644)
 if js!='none':
  deriv=np.ones_like(cx);deriv[:,64:]=np.exp(-.5*cx[:,64:]**2)/np.sqrt(2*np.pi)
  reference=np.log1p(np.sum((J*deriv[:,None])**2,axis=1))
  assert np.max(abs(feature[:,578:]-reference))<1e-12
 for centered in (False,True):
  center=bmean if centered else np.zeros_like(bmean)
  stats=fit_normalization(x-center,feature,frame)
  for kind in ('gaussian','diag_gaussian','spline','triangular','cfm'):
   config=dict(dimension=66,context_dim=feature.shape[1],width=16)
   if kind!='cfm':config.update(layers=4,output_rank=8,diag_hidden=32)
   model=make_model(kind,**config)
   p=dict(format=FORMAT,kind=kind,model_config=config,state_dict=model.state_dict(),normalization=stats,
          contract=dict(cfm_base='anchor' if kind=='cfm' else 'standard'),
          anchor_contract=dict(j_summary=js,residual_center=centered,anchor_scale=.3))
   predictor=AnchorDensityPredictor(p)
   kwargs=dict(context_J=J,base_mean=bmean)
   if kind=='cfm':
    lp=predictor.log_prob(x,y,cx,allow_approximate=True,nfe=4,**kwargs)
    expected=(-.5*(((x-bmean)/.3)**2+np.log(2*np.pi)).sum(1)-66*np.log(.3))
    samples=predictor.sample(y,cx,n_samples=3,nfe=4,batch_size=4,**kwargs)
   else:
    lp=predictor.log_prob(x,y,cx,**kwargs)
    z=((x-center)@frame-stats['x_mean'])/stats['x_scale']
    expected=-.5*(z*z+np.log(2*np.pi)).sum(1)-np.log(stats['x_scale']).sum()+stats['frame_logabsdet']
    samples=predictor.sample(y,cx,n_samples=3,batch_size=4,**kwargs)
   err=float(np.max(abs(lp-expected)))
   assert err<3e-4,(js,centered,kind,err)
   assert samples.shape==(n,3,66) and np.isfinite(samples).all()
   results[f'{js}_{centered}_{kind}']=dict(raw66_logprob_reference_max=err,sample_shape=list(samples.shape))
metric=np.broadcast_to(np.eye(66),(n,66,66)).copy()
_,check=validate_metric(metric,n)
model=make_model('cfm',width=16,context_dim=578)
a=torch.randn(n,66);b=torch.randn(n,578);mu=torch.zeros(n,66);chol=torch.eye(66)[None].repeat(n,1,1)
g1=torch.Generator().manual_seed(2);g2=torch.Generator().manual_seed(2)
unweighted=cfm_error_loss(model,a,b,g1,mu,chol)
weighted=cfm_error_loss(model,a,b,g2,mu,chol,torch.tensor(metric,dtype=torch.float32))
assert torch.allclose(unweighted,weighted,atol=1e-6)
weighted.backward()
assert all(torch.isfinite(p.grad).all() for p in model.parameters() if p.grad is not None)
try:validate_metric(metric*.1,n);raise AssertionError('out-of-bound metric accepted')
except ValueError:pass
for split,nn in [('train',8),('validation',4)]:
 np.savez_compressed(Path(__file__).parent/f'synthetic_{split}.npz',x=x[:nn],y=y[:nn],context_x=cx[:nn],base_mean=bmean[:nn],context_J=J[:nn],velocity_metric=metric[:nn],parent=np.array([f'{split}_{i}' for i in range(nn)]))
(Path(__file__).parent/'METRIC_CONTRACT.json').write_text(json.dumps(dict(kind='random',coordinate='network_standardized_full66_velocity',endpoint_label_dependent=False,TRAIN_only_tuning=True,bounds=[.5,2.],note='identity fixture only, no scientific geometry metric'),indent=2)+'\n')
results.update(status='PASS',metric_identity_equals_plain_loss=True,metric_out_of_bounds_rejected=True,
               weighted_error_gradients_finite=True,physics_calls=0,optimizer_updates=0,
               wall_seconds=time.monotonic()-wall,cpu_seconds=time.process_time()-cpu,
               peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024))
(Path(__file__).parent/'ALGEBRA_REPLAY.json').write_text(json.dumps(results,indent=2)+'\n')
print(json.dumps({k:v for k,v in results.items() if not isinstance(v,dict)},indent=2))
