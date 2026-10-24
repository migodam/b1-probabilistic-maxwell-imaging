"""Numerical tests only: no trainer call, optimizer, or physics execution."""
import json
from pathlib import Path
import resource
import sys
import time

import numpy as np
import torch

from p3r4.field_moments import *
from p3r4.train_field_moments import fit_statistics
from p3r4.probability_eval import decoder_from_geometry,decode_gaussianized

wall,cpu=time.monotonic(),time.process_time()
torch.set_num_threads(1)
rng=np.random.default_rng(17)
root=Path(__file__).resolve().parents[4]
stage=root/'research/p3l_r4_conditional_transport'
cfg=json.loads((stage/'configs/TARGET.json').read_text())
points=json.loads((root/'research/p1_material_target_20261018/inputs/geometry.json').read_text())['points']
g=decoder_from_geometry(cfg,points);roi=np.array(g['roi'])
material=decode_gaussianized(rng.normal(size=(16,66)),g)[0]
y=rng.normal(size=(16,512)).astype(np.float32)
obs=unpack_observation(y);assert np.array_equal(packed_observation(obs),y)
data=dict(y=y,material=material)
results={}
for moments in (False,True):
 stats=fit_statistics(data,roi,moments,batch_size=3)
 labels=moment_labels(material,roi,moments)
 assert np.max(abs(stats['target_mean']-labels.mean(0)))<1e-14
 assert np.max(abs(stats['y_mean']-y.mean(0,dtype=np.float64)))<1e-14
 assert stats['target_scale'][0,0]==stats['target_scale'][1,0]
 for kind in ('mlp_small','mlp_medium','cnn_small','cnn_medium'):
  model=FieldMomentRegressor(kind,256,moments)
  payload=dict(format=FORMAT,architecture=kind,second_moments=moments,geometry=g,**stats,state_dict=model.state_dict())
  predictor=FieldMomentPredictor(payload)
  out=predictor.predict_moments(y,batch_size=3)
  mean=out['material_mean']
  assert np.allclose(mean[:,roi],material[:,roi].mean(0)[None],atol=1e-14,rtol=0)
  outside=np.setdiff1d(np.arange(576),roi)
  assert np.array_equal(mean[:,outside],np.full((16,len(outside)),complex(*cfg['host_chi'])))
  true_mse=np.mean(abs(mean-material)**2,1)
  reference=np.sum(abs(mean[:,roi]-material[:,roi])**2,1)/576
  assert np.max(abs(true_mse-reference))<1e-15
  inp=torch.from_numpy(y[:2]);label=torch.from_numpy(((labels[:2]-stats['target_mean'])/stats['target_scale']).astype(np.float32))
  loss=(model(inp)-label).square().mean();loss.backward()
  assert all(torch.isfinite(p.grad).all() for p in model.parameters() if p.grad is not None)
  entry=dict(parameters=sum(p.numel() for p in model.parameters()),update0_train_mean=True,
             full576_MSE_reference_error=float(np.max(abs(true_mse-reference))),finite_gradients=True)
  if moments:
   covariance=out['raw_pointwise_covariance'][:,roi]
   centered=np.stack((material[:,roi].real-labels[:,0].mean(0),material[:,roi].imag-labels[:,1].mean(0)),axis=-1)
   reference_cov=np.einsum('nri,nrj->rij',centered,centered)/16
   cov_error=float(np.max(abs(covariance-reference_cov[None])))
   assert cov_error<1e-14
   assert np.array_equal(out['raw_pointwise_covariance'][:,outside],np.zeros((16,len(outside),2,2)))
   assert np.linalg.eigvalsh(out['psd_pointwise_covariance']).min()>-1e-12
   entry.update(train_moment_covariance_reference_error=cov_error,
                psd_projection_warning=out['psd_projection_warning'])
  else:assert 'psd_pointwise_covariance' not in out
  results[kind+('_moments' if moments else '')]=entry
bad=np.array([[[[-1.,2.],[2.,-3.]]]])
fixed,warning=project_pointwise_psd(bad)
assert np.linalg.eigvalsh(fixed).min()>-1e-12 and warning['negative_eigenvalue_count']>0
results['indefinite_covariance_projection']=dict(status='PASS',**warning)
# Exact old R2 checkpoint adapter when locally available; never rewrite it.
reference=root.parent/'b1-p3l-r2-risk-posterior/research/p3l_r2_risk_posterior_exploration/A/training/n4096/mlp_small_s11/BEST.pt'
if reference.exists():
 payload=torch.load(reference,map_location='cpu',weights_only=False)
 adapted=load_field_moments(reference)
 try:from p3r2.risk_models import RiskPredictor
 except ModuleNotFoundError:from src.p3r2.risk_models import RiskPredictor
 historical=RiskPredictor(payload)
 expected=historical.predict_full(obs[:2],batch_size=1)
 got=adapted.predict(y[:2],batch_size=1)
 assert np.array_equal(expected,got)
 assert adapted.predict_moments(y[:2])['second_moment_status']=='NOT_ESTIMATED_HISTORICAL_MEAN_ONLY'
 results['historical_R2_reference']=dict(status='PASS',path=str(reference),predictions_bitwise_equal=True,weights_retrained=False)
else:results['historical_R2_reference']=dict(status='LOCAL_CHECKPOINT_UNAVAILABLE',interface_implemented=True)
results.update(status='PASS',trainer_invoked=False,physics_calls=0,optimizer_updates=0,
               wall_seconds=time.monotonic()-wall,cpu_seconds=time.process_time()-cpu,
               peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024))
(Path(__file__).parent/'NUMERICAL_AUDIT.json').write_text(json.dumps(results,indent=2)+'\n')
print(json.dumps(results,indent=2))
