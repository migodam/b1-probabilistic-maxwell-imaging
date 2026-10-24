"""Provider schema and model compatibility only; no featuregeneration orfits."""
from pathlib import Path
import json,tempfile,time,resource,sys,importlib.util
import numpy as np,torch
from p3r4.risk_fusion import *
from p3r4.evaluation import checkpoint_identity
from p3r2.residual_risk import read_anchor
from p3r2.risk_models import observation_features
wall,cpu=time.monotonic(),time.process_time();torch.set_num_threads(1);rng=np.random.default_rng(20261025151)
root=Path.cwd();audit=root/'research/p3l_r4_conditional_transport/AUDIT/fusion/providers';work=Path(tempfile.mkdtemp(prefix='SCHEMA_ONLY_',dir=audit))
oldbase=root/'research/p3l_r4_conditional_transport/AUDIT/risk';old=torch.load(oldbase/'zero_update_interface/BEST.pt',map_location='cpu',weights_only=False);g=old['geometry'];online=read_anchor(oldbase/'check_validation_context.npz');parents=np.load(oldbase/'check_validation_context.npz')['parent'].astype(str);y=observation_features(online['observation']);n=len(parents)
raw=rng.normal(size=(n,4,24,24)).astype('float32');raw[:,2:]=abs(raw[:,2:]);model=FusionRisk(g);assert sum(p.numel() for p in model.parameters())==6418
identities={};gates={};checks=[]
for provider,count in PROVIDERS.items():
 directory=work/provider;directory.mkdir();identity=dict(path='SYNTHETIC_PROVIDER_'+provider,bytes=0,mtime_ns=0);identities[provider]=identity;artifacts=[]
 for i,parent in enumerate(parents):
  arrays=dict(parent=np.array(parent),moment_image=raw[i],K=np.array(count),checkpoint_identity=np.array(json.dumps(identity)))
  if provider!='spline':arrays['provider']=np.array(provider)
  if count:arrays['x_samples']=np.zeros((32,66))
  path=directory/f'{parent}.npz';np.savez(path,**arrays);artifacts.append(checkpoint_identity(path))
 gate=dict(status='MOMENTS_COMPLETE',truth_read=False,K=count,checkpoint_identity=identity,artifacts=artifacts)
 if provider!='spline':gate['provider']=provider
 (directory/'MOMENTS_COMPLETE.json').write_text(json.dumps(gate));gates[provider]=gate
 values=read_moments(directory,parents,identity,provider,count);assert np.array_equal(values,raw)
 payload=dict(format=FORMAT,geometry=g,control='moments',state_dict=model.state_dict(),feature_normalization=old['feature_normalization'],moment_mean=np.zeros((1,4,1,1),dtype='float32'),moment_scale=np.ones((1,4,1,1),dtype='float32'),moment_checkpoint_identity=identity)
 if provider!='spline':payload.update(moment_provider=provider,moment_K=count,moment_contract_version=2)
 path=work/(provider+'.pt');torch.save(payload,path);predictor=FusionPredictor(path)
 prediction=predictor.predict(y,online,values,parents,moment_provider=provider,moment_K=count)
 assert np.array_equal(prediction,full_material(roi_components(online['chi0'],g['roi']),g))
 if provider=='spline':assert np.array_equal(predictor.predict(y,online,values,parents),prediction)
 else:
  try:predictor.predict(y,online,values,parents)
  except ValueError:pass
  else:raise AssertionError('newprovider metadata omitted')
 for wrong_provider,wrong_count in [('spline' if provider!='spline' else 'field_moments',count),(provider,0 if count else 32)]:
  try:read_moments(directory,parents,identity,wrong_provider,wrong_count)
  except ValueError:pass
  else:raise AssertionError('wrongprovider/K accepted')
  try:predictor.predict(y,online,values,parents,moment_provider=wrong_provider,moment_K=wrong_count)
  except ValueError:pass
  else:raise AssertionError('predictor wrongprovider/K accepted')
 matching_moment_gates(gate,gate)
 checks.append(dict(provider=provider,K=count,oldSplineLegacyCompatible=provider=='spline',wrongprovider_K_rejected=True,predictor_replay_exact=True))
try:matching_moment_gates(gates['diag_gaussian'],gates['field_moments'])
except ValueError:pass
else:raise AssertionError('TRAIN/validationprovider mismatch accepted')
bad=dict(gates['field_moments'],K=32)
try:moment_contract(bad)
except ValueError:pass
else:raise AssertionError('fakefieldK32 accepted')
# Fieldprovider must be deterministic: no latentdraw arrays even ifK0 isdeclared.
directory=work/'field_moments';path=directory/f'{parents[0]}.npz'
with np.load(path) as z:arrays={k:z[k] for k in z.files}
arrays['x_samples']=np.zeros((32,66));np.savez(path,**arrays);gate=gates['field_moments'];gate['artifacts'][0]=checkpoint_identity(path);(directory/'MOMENTS_COMPLETE.json').write_text(json.dumps(gate))
try:read_moments(directory,parents,identities['field_moments'],'field_moments',0)
except ValueError:pass
else:raise AssertionError('inventedfield samples accepted')
# Mock the deterministic provider API: it must never invoke latent sampling.
spec=importlib.util.spec_from_file_location('momentgenerator',root/'tools/r4_fusion_moments.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
class FieldStub:
 def sample(self,*a,**kw):raise AssertionError('field sample must not be called')
 def predict_moments(self,*a,**kw):
  cov=np.zeros((1,576,2,2));cov[:,:,0,0]=.2;cov[:,:,1,1]=.3;cov[:,:,0,1]=cov[:,:,1,0]=.1
  return dict(material_mean=np.full((1,576),.4+.2j),psd_pointwise_covariance=cov,raw_pointwise_covariance=cov,psd_projection_warning={},second_moment_status='LEARNED_POINTWISE_JOINT_RE_IM_MOMENTS_ONLY')
image,extras=module.provider_image(FieldStub(),'field_moments',y[:1],123,None);assert image.shape==(4,24,24) and 'x_samples' not in extras and np.allclose(image[2],.2) and np.allclose(image[3],.3)
result=dict(status='PASS',checks=checks,train_validation_mismatch_rejected=True,field_fakeK32_rejected=True,invented_field_samples_rejected=True,
 deterministic_field_no_sampling=True,parameter_count=6418,legacySpline_payload_unchanged=True,physics_RHS=0,training_updates=0,actual_featuregeneration_calls=0,
 wall_seconds=time.monotonic()-wall,cpu_seconds=time.process_time()-cpu,peak_rss_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)*(1 if sys.platform=='darwin' else 1024))
(audit/'PROVIDER_COMPATIBILITY_CHECK.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
