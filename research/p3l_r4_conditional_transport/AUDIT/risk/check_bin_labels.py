"""Substantive mathematical/serialization checks; zero optimization or physics."""
import json,time,resource,sys,copy
from pathlib import Path
import numpy as np,torch
from p1_material.prior import saturation
from p3r4.risk_bin import *
from p3r4.risk import RiskActionPredictor,feature_image,full_material
from p3r2.residual_risk import read_anchor
from p3r2.risk_models import observation_features
wall,cpu=time.monotonic(),time.process_time();torch.set_num_threads(1)
rng=np.random.default_rng(20261025041)
s=np.broadcast_to(np.linspace(-2,2,1001),(64,1001)).copy()
u=rng.uniform(size=(64,2));delta=np.c_[.8+u[:,0],.03+.15*u[:,1]]
base=Path('research/p3l_r4_conditional_transport/AUDIT/risk')
payload=torch.load(base/'zero_update_interface/BEST.pt',map_location='cpu',weights_only=False);g=payload['geometry']
context=read_anchor(base/'check_validation_context.npz');y=observation_features(context['observation'])
image=feature_image(y,g,context,payload['feature_normalization'])
rows={}
for k in (8,16):
 labels=joint_bin_labels_from_s(s,delta,k);exact=delta[:,:,None]*saturation(s)[:,None,:]
 error=float(np.max(np.abs(labels.mean(2)-exact)));assert error<3e-15
 assert np.all(np.diff(labels,axis=2)<=1e-14) and np.all(labels>=0)
 assert np.all(labels<=delta[:,:,None,None]+1e-14)
 # Recover the old finite point-threshold rule on the same grid.
 t,w=np.polynomial.legendre.leggauss(k);w*=15*(1-t*t)**2/16
 pointlabels=delta[:,:,None,None]*(s[:,None,None,:]>t[None,None,:,None])
 old=np.sum(pointlabels*w[None,None,:,None],2)
 model=BinRisk('phase'+str(k),g);count=sum(p.numel() for p in model.parameters());assert count==(6114 if k==8 else 6386)
 with torch.inference_mode():mean,mom=model(torch.tensor(image),return_moments=True)
 assert torch.all(torch.diff(mom,dim=2)<=0)
 model(torch.tensor(image)).sum().backward();assert all(p.grad is not None and torch.isfinite(p.grad).all() for p in model.parameters())
 p=copy.deepcopy(payload);p.update(kind='phase'+str(k),state_dict=model.state_dict(),moment_rule=MOMENT_RULE)
 path=base/f'BIN{k}_INTERFACE_ONLY.pt';torch.save(p,path)
 predicted=RiskActionPredictor(path).predict(y,context)
 expected=full_material(mean.numpy(),g);replay=float(np.max(np.abs(predicted-expected)));assert replay==0
 labels32=labels.astype('float32');error32=float(np.max(np.abs(labels32.mean(2)-exact)))
 rows[str(k)]=dict(parameter_count=count,float64_pointwise_identity_max_error=error,float32_label_integration_max_error=error32,
                  old_point_quadrature_max_component_error=float(np.max(np.abs(old-exact))),
                  old_point_quadrature_component_mse=float(np.mean((old-exact)**2)),
                  new_bin_label_component_mse=float(np.mean((labels.mean(2)-exact)**2)),
                  serialization_replay_max_error=replay,weights_sum=float(model.weights.sum()),
                  monotonicity='PASS',bounds='PASS',finite_gradients='PASS')
result=dict(status='PASS',grid_points=1001,standardized_s_range=[-2,2],joint_attributes=64,seed=20261025041,
            checks=rows,physics_calls=0,training_updates=0,actual_validation_training='NOT_RUN_PARENT_SCHEDULES',
            wall_seconds=time.monotonic()-wall,cpu_seconds=time.process_time()-cpu,
            peak_rss_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)*(1 if sys.platform=='darwin' else 1024))
(base/'BIN_LABEL_CORRECTNESS.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
