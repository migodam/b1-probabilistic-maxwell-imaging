"""Prepare a fixed exposed24 x eleven-method R15 manifest; never launches runs."""
from __future__ import annotations
import argparse,json,time
from copy import deepcopy
from datetime import datetime,timezone
from pathlib import Path
STATEFUL=('g19','g19_s19','twofold','twofold_hop','twofold_dr','twofold_ap','twofold_full_dr','twofold_reset')
BASELINES=('cie','dual_cie','dr')
METHODS=STATEFUL+BASELINES


def write_new(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('x') as f:json.dump(value,f,indent=2,allow_nan=False);f.write('\n')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--parent-root',type=Path);p.add_argument('--campaign',default='twofold_exposed24_v1');p.add_argument('--git-identity',default=None);a=p.parse_args()
    c=time.process_time();w=time.perf_counter();root=a.root.resolve();parent=(a.parent_root or root.parent/'research_max_struct_r14').resolve()
    if not a.campaign or Path(a.campaign).name!=a.campaign or a.campaign in ('.','..') or any(x in a.campaign for x in '/\\'):raise ValueError('Safe campaign name required')
    source=parent/'configs/independent24_v1.json';original=json.loads(source.read_text());protocol=deepcopy(original['protocol'])
    if protocol['CPU_cap']!=120 or protocol['proposal_fraction']!=.6:raise ValueError('Original120/72 contract required')
    # Preserve every original protocol field; new fixed knobs are namespaced.
    specific=dict(tauS=.1,tauG=.1,numerical_rank_rtol=1e-12,period=16,omega=1.,scale=6.,DRscale=6.,AP_preparation=0,onehop_resolvent_corrections=1,
        cycle_physical_audit='Every stateful arm, same matched opportunity; implemented by parent solver',
        cycle_fullstate_snapshots='First4 cycle before/after, solver handles; fixed operator/background cache once, not per checkpoint',
        periodic_physical_audit_CPU_seconds=2.,checkpoint_iterations=[1,32,128,512,1500],checkpoint_stage_end=True,checkpoint_final=True)
    protocol['r15']=specific
    by={(t['case_id'],t['method']):t for t in original['tasks']};cases=[f'r14_new_{i:04d}' for i in range(24)];tasks=[];public=[]
    for cid in cases:
        path=root/'data_public/independent24'/(cid+'.npz')
        if not path.is_file():raise FileNotFoundError('Parent must copy unchanged public input and record provenance before freeze: '+str(path))
        public.append(dict(case_id=cid,path=str(path),bytes=path.stat().st_size,source_reference=str(parent/'data_public/independent24'/(cid+'.npz'))))
        for method in METHODS:
            template=by[cid,method if method in BASELINES else 'g19'];t=deepcopy(template)
            t.update(task_id=cid+'__'+method,campaign=a.campaign,method=method,algorithm=method,cohort='EXPOSED24_FOLLOWUP',
                input_file='data_public/independent24/'+cid+'.npz',runner_family='FROZEN_NATIVE_C2' if method in BASELINES else 'PARENT_STATEFUL_G19_TWOFOLD',CPU_cap=120)
            t['observation_contract']['input_file']=t['input_file'];tasks.append(t)
    if len(tasks)!=264 or len({t['task_id'] for t in tasks})!=264:raise ValueError('Full264 frozen matrix required')
    value=dict(campaign=a.campaign,status='FROZEN_BEFORE_R15_PREDICTIONS',exposure_status='EXPOSED24_FOLLOWUP',new_materials=False,new_independent_or_blind_claim=False,
        source_manifest=str(source),source_campaign=original['campaign'],protocol=protocol,original_protocol_exact_except_namespaced_r15={k:v for k,v in protocol.items() if k!='r15'}==original['protocol'],
        methods=list(METHODS),stateful_methods=list(STATEFUL),native_baselines=list(BASELINES),cases=cases,expected_tasks=264,git_identity=a.git_identity,
        git_identity_null_reason=None if a.git_identity else 'Not supplied; no identity fabricated',material_real_DOF=576,
        truth_access='NONE',input_policy='Parent-copied unchanged public R14 independent24 observations at R15-local data_public/independent24; separate copy provenance required',
        selector='Original paid legal/eligible minPhi; original minR-CN99 secondary, fallback primary only in offline deployment reporting',
        stateful_scope='Explicit R15 G19/current-state transplant; native frozen CIE/Dual/DR separate controls, not claimed identical finite trajectories',tasks=tasks)
    target=root/'configs'/(a.campaign+'.json');receipt=root/'audit'/('MANIFEST_FREEZE_'+a.campaign+'.json')
    if target.exists() or receipt.exists():raise FileExistsError('Never overwrite a frozen manifest/receipt')
    write_new(target,value);write_new(receipt,dict(status='PREPARED_NOT_EXECUTED',created_utc=datetime.now(timezone.utc).isoformat(),manifest=str(target),source=str(source),source_bytes=source.stat().st_size,
        public_inputs=public,cases=24,tasks=264,methods=list(METHODS),protocol_exact_except_namespaced_r15=value['original_protocol_exact_except_namespaced_r15'],
        git_identity=a.git_identity,truth_access='NONE',physics_calls=0,launches=0,CPU_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-w))
    print(json.dumps(dict(manifest=str(target),tasks=264,truth_access='NONE',launches=0)))

if __name__=='__main__':main()
