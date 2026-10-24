#!/usr/bin/env python3
"""Wait at60s cadence, then run authorized frozen risk online steps only.

Separate from all frozen inference/model sources. No GT/scoring/physical launch,
no automatic retries, and no inference before complete global context marker.
"""
from pathlib import Path
import json,sys,time,subprocess,os,traceback
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3r4.evaluation import checkpoint_identity,write_new_json
S=ROOT/'research/p3l_r4_conditional_transport';F=S/'results/fresh';COMP=S/'AUDIT/risk/fresh_freeze_v3/RISK_COMPONENT_FROZEN.json';MARKER=S/'results/FRESH_PHYSICS_AND_CONTEXTS_COMPLETE.json';OUT=F/'risk_online_queue_v1'
def main():
 if OUT.exists():raise FileExistsError('one background risk queue only; no automatic retry')
 OUT.mkdir(parents=True);commands=[];releases={};start=time.time()
 for cohort in ['locked','confirmation']:
  top=S/'configs'/('FROZEN_'+cohort.upper()+'.json');release=S/'configs'/('FRESH_RISK_'+cohort.upper()+'_ONLINE_RELEASE_V1.json')
  spec=dict(released=True,experiment='fresh_risk_panels_v1',cohort=cohort,component_identity=checkpoint_identity(COMP),top_freeze_identity=checkpoint_identity(top),allow_features_infer=True,allow_offline_score=False,authorization='parent20261025: background zero-physics features/infer/export only after globalCOMPLETE marker; no heldoutGT/no newfits/no physics',prerequisite_marker=str(MARKER.resolve()))
  if release.exists():
   if json.loads(release.read_text())!=spec:raise ValueError('existing release differs; do not overwrite')
  else:write_new_json(release,spec)
  releases[cohort]=checkpoint_identity(release)
  cmd=[sys.executable,str(ROOT/'tools/queue_r4_fresh_risk.py'),'run','--cohort',cohort,'--component',str(COMP),'--top-freeze',str(top),'--freeze',str(F/cohort/'BOUND_RISK.json'),'--release',str(release),'--out',str(F/cohort/'risk45'),'--stage','IV'];commands.append(cmd)
 write_new_json(OUT/'WAIT_FROZEN.json',dict(status='WAIT_FOR_GLOBAL_COMPLETE',pid=os.getpid(),start_unix=start,marker=str(MARKER),cadence_seconds=60,component_identity=checkpoint_identity(COMP),release_identities=releases,commands=commands,heldout_GT_allowed=False,physical_calls_allowed=False,training_allowed=False))
 while not MARKER.exists():
  print(json.dumps(dict(status='WAITING_FOR_GLOBAL_PHYSICS_CONTEXT_COMPLETE',utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()))),flush=True);time.sleep(60)
 gate=json.loads(MARKER.read_text())
 if gate.get('status')!='COMPLETE':raise ValueError('global marker exists but not COMPLETE; preserve failure, no retry')
 write_new_json(OUT/'PREREQUISITE_COMPLETE.json',dict(status='GLOBAL_COMPLETE_OBSERVED',marker_identity=checkpoint_identity(MARKER),observed_unix=time.time()))
 for cmd in commands:subprocess.run(cmd,check=True,cwd=ROOT)
 # Actual observed risk outputs are sealed by ONLINE_COMPLETE. Preparing a
 # representative request reads only those outputs and existing online caches.
 manifest=F/'locked/risk45/exports/REPRESENTATIVE_MANIFEST.json';prepared=F/'locked/representatives_prepared'
 cmd=[sys.executable,str(ROOT/'tools/run_r4_job.py'),'--stage','IV','--job','fresh_locked_representatives_prepare_v1','--max-wall','120','--rss-limit-mib','2048','--',sys.executable,str(ROOT/'tools/r4_fresh_representatives.py'),'prepare','--manifest',str(manifest),'--stage','IV','--out',str(prepared)]
 subprocess.run(cmd,check=True,cwd=ROOT)
 # Inactive release templates only. Parent owns any physical launch/caps/wall.
 for cohort in ['locked','confirmation']:
  exports=F/cohort/'risk45/exports';raw=exports/'RAW_FORWARD_FROZEN.json';r=json.loads(raw.read_text())
  write_new_json(exports/'RAW_FORWARD_RELEASE_TEMPLATE.json',dict(released=False,experiment='raw_material_forward_diagnostic_v1',stage='IV',freeze_identity=checkpoint_identity(raw),planned_RHS=r['planned_RHS'],max_wall='PARENT_MUST_SET',scope='arbitraryimage response only; decoder-realizability NOT_PROVEN; no GT'))
 write_new_json(prepared/'PHYSICAL_RELEASE_TEMPLATE.json',dict(released=False,experiment='fresh_locked64_representatives_v1',stage='IV',prepared_request_identity=checkpoint_identity(prepared/'FROZEN_REQUEST.json'),maximum_total_RHS=81920,max_wall='PARENT_MUST_SET',scope='fixedSpline_s1+R2raw/clipped, sharedordinary chargedonce; noGT'))
 write_new_json(OUT/'COMPLETE.json',dict(status='BOTH_COHORT_RISK_ONLINE_COMPLETE',completed_unix=time.time(),prerequisite_marker_identity=checkpoint_identity(MARKER),parents=dict(locked=512,confirmation=256),methods=45,outputs={cohort:str((F/cohort/'risk45').resolve()) for cohort in ['locked','confirmation']},representative_prepared=str((prepared/'FROZEN_REQUEST.json').resolve()),heldout_GT_opened=False,scoring_run=False,new_fits=0,Maxwell_RHS=0,all_risk_and_density_inference_required_before_first_score='parent global gate; this worker grants no scoring permission'))
 print(json.dumps(dict(status='BOTH_COHORT_RISK_ONLINE_COMPLETE',path=str(OUT/'COMPLETE.json'))),flush=True)
if __name__=='__main__':
 try:main()
 except BaseException as exc:
  if OUT.exists() and not (OUT/'ERROR.json').exists():write_new_json(OUT/'ERROR.json',dict(status='FAILED_NO_AUTORETRY',error=repr(exc),traceback=traceback.format_exc(),heldout_GT_opened=False))
  raise
