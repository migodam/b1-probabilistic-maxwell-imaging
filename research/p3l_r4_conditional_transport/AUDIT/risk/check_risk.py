"""Zero-physics structural and numerical checks; never trains."""
import json,time,resource,sys
from pathlib import Path
import numpy as np,torch
from scipy.special import ndtri
from p3r4.risk import *
from p3r2.residual_risk import read_anchor
wall,cpu=time.monotonic(),time.process_time();torch.set_num_threads(1)
root=Path('../b1-p3l-r2-risk-posterior/research/p3l_r2_risk_posterior_exploration')
cp=root/'A/residual_training/n512/local_geometry_cnn_s33/BEST.pt';payload=torch.load(cp,map_location='cpu',weights_only=False);g=payload['geometry']
a=read_anchor(root/'A/anchor_data/validation.npz');a={k:v[:2] for k,v in a.items()}
y=observation_features(a['observation']);image=feature_image(y,g,a,payload['normalization'])
base=roi_components(a['chi0'],g['roi']);r2=load_r2_baseline(cp)
# Exact historical weights map one-to-one to the common local trunk and head.
m=LocalRisk('local_cnn',g);m.local.load_state_dict(r2.model.local.state_dict());m.head.load_state_dict(r2.model.head.state_dict())
with torch.inference_mode():out=m(torch.tensor(image),torch.tensor(base)).numpy()
# R2 residual_scale is checkpoint-specific and must be retained in exact replay.
ref=roi_components(r2.predict_full(a),g['roi'])
scaled=base+(out-base)*payload['residual_scale']
replay_error=float(np.max(np.abs(scaled-ref)));assert replay_error<2e-7
with np.load(root/'A/data/validation.npz') as z:q=z['q'][:3];truth=z['chi'][:3]
x=np.c_[q[:,:64],ndtri(q[:,64:])]
checks={};roi=g['roi'];cfg=g['config']
for kind in ('local_cnn','no_geometry_cnn','phase8','phase16','direct_moment'):
 model=LocalRisk(kind,g);count=sum(p.numel() for p in model.parameters());assert count<=30000
 row={'parameters':count}
 if model.k:
  mean,mom=model(torch.tensor(image),return_moments=True);assert torch.all(torch.diff(mom,dim=2)<=0)
  assert torch.all(mom>=0) and torch.all(mom<=model.bounds[None,:,None,None])
  lab=phase_labels(x,g,model.thresholds.numpy());assert np.all(np.diff(lab,axis=2)<=0)
  approx=np.array(cfg['host_chi'])[None,:,None]+np.sum(lab*model.weights.numpy()[None,None,:,None],2)
  exact=roi_components(truth,roi);row.update(quadrature_label_max_error=float(np.max(np.abs(approx-exact))),weights_sum=float(model.weights.sum()))
  mean.sum().backward();assert all(p.grad is not None and torch.isfinite(p.grad).all() for p in model.parameters())
 checks[kind]=row
result={'status':'PASS','physics_calls':0,'training_updates':0,'exact_r2_scaled_architecture_error':replay_error,'checks':checks,'wall_seconds':time.monotonic()-wall,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)*(1 if sys.platform=='darwin' else 1024)}
Path('research/p3l_r4_conditional_transport/AUDIT/risk/CORRECTNESS.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
