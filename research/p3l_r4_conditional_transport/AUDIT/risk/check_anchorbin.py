"""Decoder-only revision2 math/replay checks; no fitting or physics."""
from pathlib import Path
import json,time,resource,sys
import numpy as np,torch
from p3r4.risk_anchorbin import *
from p3r4.risk_bin import joint_bin_labels_from_s
from p1_material.prior import MaterialPrior,saturation
wall,cpu=time.monotonic(),time.process_time();torch.set_num_threads(1);rng=np.random.default_rng(20261025045)
base=Path('research/p3l_r4_conditional_transport/AUDIT/risk');payload=torch.load(base/'zero_update_interface/BEST.pt',map_location='cpu',weights_only=False);g=payload['geometry'];cfg=g['config']
points=json.loads(Path('research/p1_material_target_20261018/inputs/geometry.json').read_text())['points'];prior=MaterialPrior(cfg,points)
q=prior.draw(rng,32);chi=np.array([prior.decode(a)[0] for a in q]);roi=np.array(g['roi']);s=(cfg['field_mean']+q[:,:64]@np.asarray(g['A'])[roi].T)/cfg['phase_transition_halfwidth']
delta=np.stack((.8+q[:,64],.03+.15*q[:,65]),1);image=torch.tensor(rng.normal(size=(32,8,24,24)),dtype=torch.float32);qt=torch.tensor(q,dtype=torch.float32);ct=torch.tensor(chi,dtype=torch.complex64)
rows={}
for k in (8,16):
 model=AnchorBinRisk('phase'+str(k),g);count=sum(p.numel() for p in model.parameters());assert count==(6114 if k==8 else 6386)
 masses,baseline,phase=model.baseline(qt,ct)
 assert torch.all(masses>0) and torch.allclose(masses.sum(2),torch.ones_like(masses[:,:,0]),atol=2e-7)
 target=torch.tensor(joint_bin_labels_from_s(s,delta,k),dtype=torch.float32)
 label_error=float(torch.max(torch.abs(baseline-target)));assert label_error<2e-5
 with torch.inference_mode():mean,mom=model(image,qt,ct,True)
 assert torch.all(torch.diff(mom,dim=2)<=0) and torch.all(mom>=0) and torch.all(mom<=model.bounds[None,:,None,None]+1e-7)
 anchor=torch.stack((ct[:,roi].real,ct[:,roi].imag),1)
 expected=(1-EPSILON)*anchor+EPSILON*(model.host[None,:,None]+model.bounds[None,:,None]/2)
 init_error=float(torch.max(torch.abs(mean-expected)));assert init_error<8*torch.finfo(torch.float32).eps*float(model.bounds.max())
 full=full_material(mean.numpy(),g);rms=float(np.sqrt(np.mean(abs(full-chi)**2)));assert rms<1e-4
 # Full tilt head gradients are finite; label formula is unchanged.
 output,mom=model(image,qt,ct,True);loss=((mom-target)**2).mean();loss.backward();assert all(p.grad is not None and torch.isfinite(p.grad).all() for p in model.parameters())
 zero=torch.zeros((2,2));host=model.host;bounds=model.bounds
 zm,zt,zphase=anchor_mass_from_components(host[None,:,None].expand(2,2,7).clone(),zero,host,bounds,k)
 assert torch.all(zt==0) and torch.all(torch.isfinite(torch.log(zm))) and torch.allclose(zm.sum(2),torch.ones((2,2,7)))
 for scale in (1e-2,1e-4,1e-8):
  small=torch.full((2,2),scale,dtype=torch.float64);small_anchor=host.double()[None,:,None]+small[:,:,None]*torch.full((2,2,7),.37,dtype=torch.float64)
  sm,st,_=anchor_mass_from_components(small_anchor,small,host.double(),bounds.double(),k)
  assert torch.isfinite(torch.log(sm)).all() and float(st.max())<=scale*(1+1e-10)
 ztail=torch.flip(torch.cumsum(torch.flip(zm[:,:,1:],dims=[2]),dim=2),dims=[2])*bounds[None,:,None,None]
 zmean=ztail.mean(2);assert torch.allclose(zmean,EPSILON*bounds[None,:,None]/2,atol=2e-8)
 # Format-specific predictor must preserve old loader routes and raw anchors.
 p=dict(payload);p.update(format=FORMAT,kind='phase'+str(k),anchor_epsilon=EPSILON,state_dict=model.state_dict())
 path=base/f'ANCHORBIN{k}_INTERFACE_ONLY.pt';torch.save(p,path);predictor=AnchorBinPredictor(path)
 from p3r2.residual_risk import read_anchor
 from p3r2.risk_models import observation_features
 online=read_anchor(base/'check_validation_context.npz');yy=observation_features(online['observation']);im=feature_image(yy,g,online,p['feature_normalization'])
 with torch.inference_mode():reference=full_material(model(torch.tensor(im),torch.tensor(online['q0'],dtype=torch.float32),torch.tensor(online['chi0'],dtype=torch.complex64)).numpy(),g)
 replay=predictor.predict(yy,online);assert np.array_equal(reference,replay)
 try:model(image)
 except ValueError:pass
 else:raise AssertionError('missing raw anchor accepted')
 broken=ct.clone();broken[:,roi]+=0.01j
 try:model.baseline(qt,broken)
 except ValueError:pass
 else:raise AssertionError('inconsistent Re/Im phase accepted')
 rows[str(k)]=dict(parameter_count=count,baseline_label_max_error=label_error,update0_smoothed_anchor_max_error=init_error,update0_full_RMS_from_anchor=rms,
                  baseline_masses_sum_max_error=float(torch.max(torch.abs(masses.sum(2)-1))),monotonic_bounds='PASS',finite_gradients='PASS',
                  zero_delta_joint_moments_zero='PASS',near_zero_delta_finite_bounded_limit='PASS',zero_delta_smoothed_epsBhalf_limit='PASS',raw_anchor_required='PASS',ReIm_disagreement_rejected='PASS',predictor_replay_max_error=float(np.max(abs(reference-replay))))
result=dict(status='PASS',checks=rows,epsilon=EPSILON,raw_decoder_anchors=32,physics_calls=0,training_updates=0,actual_validation='NOT_RUN_PARENT_SCHEDULES',
 wall_seconds=time.monotonic()-wall,cpu_seconds=time.process_time()-cpu,peak_rss_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)*(1 if sys.platform=='darwin' else 1024))
(base/'ANCHORBIN_REVISION2_CHECK.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
