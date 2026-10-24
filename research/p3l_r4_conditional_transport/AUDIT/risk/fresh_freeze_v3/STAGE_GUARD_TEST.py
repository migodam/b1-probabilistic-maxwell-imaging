from pathlib import Path
import sys,json,tempfile
from unittest.mock import patch
from types import SimpleNamespace
ROOT=Path("/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r4-conditional-transport");sys.path[:0]=[str(ROOT/'src'),str(ROOT/'tools')]
import r4_fresh_representatives as rep,r4_raw_material_forward as raw
from p3r4.evaluation import checkpoint_identity,write_new_json
from r4_fresh_risk import component_check
S=ROOT/'research/p3l_r4_conditional_transport';A=S/'AUDIT/risk/fresh_freeze_v3';c=json.loads((A/'RISK_COMPONENT_FROZEN.json').read_text());old=json.loads((S/'AUDIT/risk/fresh_freeze_v2/RISK_COMPONENT_FROZEN.json').read_text());component_check(c)
assert c['risk_rule']==old['risk_rule'] and c['providers']==old['providers']
assert [(m['id'],m.get('checkpoint_identity')) for m in c['methods']]==[(m['id'],m.get('checkpoint_identity')) for m in old['methods']]
checks=[]
class Intercept(Exception):pass
called=[]
def no_book(stage,*args,**kwargs):called.append(stage);raise Intercept('Book intercepted before construction; no ledger/physics')
def reject(fn):
 try:fn()
 except ValueError:return
 raise AssertionError('stage mismatch was not rejected')
with tempfile.TemporaryDirectory() as d:
 d=Path(d);manifest=d/'m.json';manifest.write_text('{}');req=d/'p.json';rel=d/'release.json';spec=d/'raw.json'
 request=dict(stage='IV',manifest=dict(stage='IV',parents=['R4L1_00000']),manifest_identity=checkpoint_identity(manifest),source_identities=[],qualified=[],maximum_total_RHS=81920)
 req.write_text(json.dumps(request));release=dict(stage='IV',released=True,experiment='fresh_locked64_representatives_v1',prepared_request_identity=checkpoint_identity(req));rel.write_text(json.dumps(release))
 a=SimpleNamespace(stage='II',prepared_request=req,release=rel,job='test_no_run',max_wall=60,out=d/'repout')
 reject(lambda:rep.run(a));checks.append('representative CLI II vs IV request rejected before Book')
 a.stage='IV';release['stage']='II';rel.write_text(json.dumps(release));reject(lambda:rep.run(a));checks.append('representative release II vs IV request rejected before Book')
 release['stage']='IV';rel.write_text(json.dumps(release))
 with patch.object(rep,'online_check',lambda p:({},{})),patch('p3r4.core.Book',no_book):
  try:rep.run(a)
  except Intercept:pass
  else:raise AssertionError('Book not intercepted')
 assert called[-1]=='IV';checks.append('matching representative stage forwarded IV; Book construction intercepted')
 with patch.object(rep,'online_check',lambda p:({'stage':'IV'},{})):
  reject(lambda:rep.prepare(SimpleNamespace(stage='II',manifest=manifest)));checks.append('fresh representative prepare requires IV')
 for stage,parent,argstage,relstage,expected in [('IV','R4L1_00000','II','IV','reject'),('IV','R4L1_00000','IV','II','reject'),('II','R4L1_00000','II','II','reject'),('IV','R4L1_00000','IV','IV','IV'),(None,'R4DEV_00000','II',None,'II')]:
  payload=dict(status='FROZEN_ARBITRARY_MATERIAL_FORWARD',rows=[],declared_parents=[parent],planned_RHS=0)
  if stage:payload['stage']=stage
  spec.write_text(json.dumps(payload));r=dict(released=True,experiment=raw.EXPERIMENT)
  if relstage:r['stage']=relstage
  rel.write_text(json.dumps(r));a=SimpleNamespace(stage=argstage,freeze=spec,release=rel,job='test_no_run',max_rhs=0,max_wall=60,out=d/'rawout')
  if expected=='reject':reject(lambda:raw.run(a))
  else:
   with patch('p3r4.core.Book',no_book):
    try:raw.run(a)
    except Intercept:pass
    else:raise AssertionError('Book not intercepted')
   assert called[-1]==expected
  checks.append('raw stage contract '+str((stage,parent,argstage,relstage,expected)))
 physical=d/'physical';physical.mkdir();(physical/'RESULTS.json').write_text(json.dumps(dict(stage='II')));receipt=d/'CLOSED.json';receipt.write_text('{}');(d/'START.json').write_text('{}')
 a=SimpleNamespace(stage='IV',prepared_request=req,physical=physical,receipt=receipt,out=d/'score',truth=d/'offline_only')
 with patch.object(rep,'online_check',lambda p:({},{})):reject(lambda:rep.score(a))
 checks.append('score rejects physical II for requested IV before truth access')
write_new_json(A/'STAGE_GUARD_CHECKS.json',dict(status='ZERO_PHYSICS_CONTRACT_CHECKS',checks=checks,Book_constructions=0,Book_arguments_intercepted=called,Maxwell_RHS=0,GT_reads=0,observations_generated=0,same45models_providers_validation_rule=True,hardcaps=dict(per_arm_states=16,per_arm_RHS=256,representatives=81920,rawL1=24576,rawConfirmation=12288)))
print('Stage guards checked; all Book calls intercepted, 0 physical calls')
