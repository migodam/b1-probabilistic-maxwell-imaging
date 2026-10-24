import sys,json,tempfile,zipfile,pathlib
from unittest.mock import patch
import numpy as np
ROOT=pathlib.Path("/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r4-conditional-transport");sys.path[:0]=[str(ROOT/'src'),str(ROOT/'tools')]
from p3r4.evaluation import write_new_json,checkpoint_identity
from p3r4.fresh_risk import canonical_infer
from r4_fresh_risk import component_check,bound
from types import SimpleNamespace
A=ROOT/'research/p3l_r4_conditional_transport/AUDIT/risk/fresh_freeze_v2';c=json.loads((A/'RISK_COMPONENT_FROZEN.json').read_text());component_check(c)
assert len(c['methods'])==45 and [v['K'] for v in c['providers'].values()]==[32,32,0]
assert c['cohorts']['locked']['n']==512 and c['cohorts']['confirmation']['n']==256
for cohort in ['locked','confirmation']:
 try:bound(SimpleNamespace(component=A/'RISK_COMPONENT_FROZEN.json',top_freeze=A/('TOP_LEVEL_'+cohort.upper()+'_TEMPLATE.json'),cohort=cohort))
 except ValueError:pass
 else:raise AssertionError('inactive template must reject')
parents=[f'SYNTHETIC_{i:05d}' for i in range(64)];f=dict(methods=[dict(id='good',loader='r4',checkpoint='unused'),dict(id='broken',loader='r4',checkpoint='unused_bad'),dict(id='anchor',loader='anchor')],parents=parents,providers={},canonical_batch_size=64,array_codec='stored')
class Model:
 geometry=dict(config={},roi=np.arange(576),A=np.zeros((576,64)))
 def __init__(self,bad=False):self.bad=bad
 def predict(self,y,anchor,batch_size):
  assert batch_size==64 and len(y)==64
  if self.bad:raise ValueError('synthetic failure retained')
  return np.arange(64*576).reshape(64,576)*1e-6+1j*np.ones((64,576))*.01
on=dict(observation=np.zeros((1,32,8),complex),residual=np.zeros((1,32,8),complex),q0=np.zeros((1,66)),chi0=np.full((1,576),.2+.1j),J=np.zeros((1,512,66)))
with tempfile.TemporaryDirectory() as d:
 d=pathlib.Path(d);ctx=d/'contexts';ctx.mkdir();features=d/'features';features.mkdir()
 obs=d/'observations';obs.mkdir()
 for p in parents:
  np.savez(ctx/(p+'.npz'),parent=np.array(p))
  np.savez(obs/(p+'.npz'),observation=np.zeros((32,8),complex))
 write_new_json(features/'FEATURES_COMPLETE.json',dict(freeze=f,gates=[],observation_sources=[checkpoint_identity(obs/(p+'.npz')) for p in parents]))
 with patch('p3r4.fresh_risk.read_online_context',lambda p,parent:on),patch('p3r4.risk_anchorbin.load_risk_compatible',lambda path,device:Model(path.endswith('bad'))),patch('p3r2.fusion.passive_box',lambda pred,cfg,roi:pred.copy()):
  canonical_infer(f,ctx,d/'native',features,d/'out',lambda f:None)
 with np.load(d/'out/good.npz') as z, np.load(d/'out/batches/good/batch_000.npz') as b:
  assert np.array_equal(z['prediction_raw'],b['prediction_raw']);assert z['parent'].astype(str).tolist()==parents
 with np.load(d/'out/broken.npz') as z:assert z['fallback_anchor'].sum()==64 and np.array_equal(z['prediction_raw'],np.repeat(on['chi0'],64,0))
 with zipfile.ZipFile(d/'out/good.npz') as z:assert all(k.compress_type==zipfile.ZIP_STORED for k in z.infolist())
 gate=json.loads((d/'out/ONLINE_COMPLETE.json').read_text());assert len(gate['failures']['broken'])==64 and gate['truth_read'] is False and gate['physics_calls']==0
write_new_json(A/'CONTRACT_CHECKS.json',dict(status='SYNTHETIC_CONTRACT_CHECKS_ONLY',actual_fresh_observations=0,GT_reads=0,physics_RHS=0,training_updates=0,component_sources_verified=True,methods=45,providers_K=[32,32,0],inactive_top_templates_rejected=2,synthetic_batch64_stored_roundtrip_bitwise=True,synthetic_failure64_preserved=True,source_context_parent_guard=True,sealed_original_observation_contract=True,representative_total_RHS_cap=64*5*256,raw_diagnostic_RHS_512=512*6*8,raw_diagnostic_RHS_256=256*6*8,scope='serialization/batching/denominator/release contract only; no scientific validation, no change to originalFAIL'))
print('Synthetic contract checks completed; no fresh observation or physical work')
