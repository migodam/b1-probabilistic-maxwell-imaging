"""No-fit fusion architecture/permutation/serialization checks."""
from pathlib import Path
import json,time,resource,sys
import numpy as np,torch
from p3r4.risk_fusion import *
from p3r4.evaluation import checkpoint_identity
from p3r2.residual_risk import read_anchor
from p3r2.risk_models import observation_features
wall,cpu=time.monotonic(),time.process_time();torch.set_num_threads(1);rng=np.random.default_rng(20261025123)
base=Path('research/p3l_r4_conditional_transport/AUDIT/fusion');momentdir=base/'SCHEMA_MOMENTS';momentdir.mkdir(exist_ok=True)
oldbase=Path('research/p3l_r4_conditional_transport/AUDIT/risk');old=torch.load(oldbase/'zero_update_interface/BEST.pt',map_location='cpu',weights_only=False);g=old['geometry'];online=read_anchor(oldbase/'check_validation_context.npz');parents=np.load(oldbase/'check_validation_context.npz')['parent'].astype(str);n=len(parents)
raw=rng.normal(size=(n,4,24,24)).astype('float32');raw[:,2:]=abs(raw[:,2:]);identity=dict(path='SYNTHETIC_MODEL_NOT_A_CHECKPOINT',bytes=0,mtime_ns=0);artifacts=[]
for i,parent in enumerate(parents):
 path=momentdir/f'{parent}.npz';np.savez(path,parent=np.array(parent),moment_image=raw[i],K=np.array(32),checkpoint_identity=np.array(json.dumps(identity)))
 artifacts.append(checkpoint_identity(path))
(momentdir/'MOMENTS_COMPLETE.json').write_text(json.dumps(dict(status='MOMENTS_COMPLETE',truth_read=False,K=32,checkpoint_identity=identity,artifacts=artifacts)))
read=read_moments(momentdir,parents,identity);assert np.array_equal(raw,read)
mean=raw.mean((0,2,3),keepdims=True);scale=raw.std((0,2,3),keepdims=True);a=transform_moments(raw,mean,scale,parents,'moments');b=transform_moments(raw,mean,scale,parents,'permuted')
assert np.array_equal(np.sort(a.reshape(n,4,576),axis=2),np.sort(b.reshape(n,4,576),axis=2)) and not np.array_equal(a,b)
assert np.array_equal(b,transform_moments(raw,mean,scale,parents,'permuted'))
model=FusionRisk(g);assert sum(v.numel() for v in model.parameters())==6418
image=torch.tensor(rng.normal(size=(n,8,24,24)),dtype=torch.float32);anchor=torch.tensor(roi_components(online['chi0'],g['roi']))
assert torch.equal(model(image,torch.tensor(a),anchor),anchor)
model(image,torch.tensor(a),anchor).square().mean().backward();assert all(p.grad is not None and torch.isfinite(p.grad).all() for p in model.parameters())
payload=dict(format=FORMAT,geometry=g,control='moments',state_dict=model.state_dict(),feature_normalization=old['feature_normalization'],moment_mean=mean,moment_scale=scale,moment_checkpoint_identity=identity)
path=base/'FUSION_INTERFACE_ONLY.pt';torch.save(payload,path);predictor=FusionPredictor(path);y=observation_features(online['observation']);prediction=predictor.predict(y,online,raw,parents);reference=full_material(anchor.numpy(),g);assert np.array_equal(prediction,reference)
try:predictor.predict(y,online)
except ValueError:pass
else:raise AssertionError('online momentarrays not required')
result=dict(status='PASS',parameter_count=6418,zerohead_anchor_exact=True,perparent_permutation_same_amplitudes=True,permutation_deterministic=True,predictor_replay_exact=True,finite_gradients=True,online_moments_required=True,
 physics_RHS=0,training_updates=0,featuregeneration_launched=False,wall_seconds=time.monotonic()-wall,cpu_seconds=time.process_time()-cpu,peak_rss_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)*(1 if sys.platform=='darwin' else 1024))
(base/'FUSION_MATH_CHECK.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
