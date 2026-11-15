"""Finite fixed56 full-state interventions; invoke only after main queue closure.

No retries, daemon, GT, scoring or changes to the intervention mathematics.
Four task lanes0,2,4,6, one BLAS thread, NORMAL; full-child external CPU.
"""
from __future__ import annotations
import os
THREAD_KEYS=('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS')
for _k in THREAD_KEYS:os.environ[_k]='1'
import argparse,json,importlib.util,re,subprocess,sys,time,traceback
from pathlib import Path
from collections import deque
from datetime import datetime,timezone
MAIN='twofold_exposed24_v1'
CAMPAIGN='fullstate_interventions_v1'
CORES=(0,2,4,6)
METHODS=('twofold_hop','twofold_dr','twofold_ap','twofold_full_dr','twofold_reset')


def read(path):return json.loads(Path(path).read_text())
def utc():return datetime.now(timezone.utc).isoformat()
def write_new(path,value):
    with Path(path).open('x') as f:json.dump(value,f,indent=2,allow_nan=False);f.write('\n')


def complete_checkpoint(path):
    """Inspect typed JSON/index metadata only, never load numerical state arrays."""
    try:
        import numpy as np
        with np.load(path,allow_pickle=False) as z:
            tree=json.loads(str(z['__state_tree_json__'].item()));members=set(z.files)
            if tree.get('type')!='dict':raise ValueError('Missing typed state dictionary')
            state=dict(tree['items'])
            required={'m','A','J','Q','old_grad','old_direction','stage','iteration','iterations','r15'}
            if not required<=set(state):raise ValueError('Incomplete PR/current/material state fields')
            r15=dict(state['r15']['items'])
            if r15.get('pending_cycle') is not None:raise ValueError('Pending cycle is not a completed checkpoint')
            for key in ('method','memory_X','memory_Q','cycle_index','fixed_cache_relative'):
                if key not in r15:raise ValueError('Missing R15 complete-state field '+key)
            if r15['method']!='twofold_hop':raise ValueError('Source checkpoint is not shared twofold_hop')
            def arrays(node):
                if isinstance(node,dict):
                    if node.get('type')=='ndarray' and node.get('key') not in members:raise ValueError('Missing typed array entry')
                    for v in node.values():arrays(v)
                elif isinstance(node,list):
                    for v in node:arrays(v)
            arrays(tree)
            fixed=path.parent/r15['fixed_cache_relative']
            if not fixed.is_file():raise FileNotFoundError('Required fixed cache missing: '+str(fixed))
            background=r15.get('background_cache_relative')
            if background is not None and not (path.parent/background).is_file():raise FileNotFoundError('Background cache missing')
        return dict(complete=True,typed_schema='R14 recursive state plus R15 PR/reflection metadata',metadata_only=True,bytes=path.stat().st_size,fixed_cache=str(fixed),background_cache=str(path.parent/background) if background else None)
    except Exception as exc:return dict(complete=False,error=repr(exc),metadata_only=True)


def source_checkpoint(root,cid):
    folder=root/'RESULTS'/MAIN/(cid+'__twofold_hop');preferred=folder/'CHECKPOINT_00032.npz';checked=[]
    if preferred.is_file():
        info=complete_checkpoint(preferred);checked.append(dict(path=str(preferred),**info))
        if info['complete']:return dict(path=str(preferred),requested_iteration=32,chosen_iteration=32,fallback=False,checked=checked)
    candidates=[]
    for p in folder.glob('CHECKPOINT_*.npz'):
        match=re.fullmatch(r'CHECKPOINT_(\d+)\.npz',p.name)
        if match and p!=preferred:candidates.append((int(match.group(1)),p))
    for iteration,p in sorted(candidates):
        info=complete_checkpoint(p);checked.append(dict(path=str(p),**info))
        if info['complete']:return dict(path=str(p),requested_iteration=32,chosen_iteration=iteration,fallback=True,fallback_reason='Requested32 missing/incomplete; earliest complete numeric checkpoint; CYCLE_BEFORE/STAGE/FINAL excluded',checked=checked)
    return dict(path=None,requested_iteration=32,chosen_iteration=None,fallback=True,fallback_reason='No complete numeric source checkpoint',checked=checked)


def geometry_map(root):
    folder=root/'audit/rank_all24/real';mapping={}
    for p in sorted(folder.glob('EXECUTION_*.json')):
        row=read(p)
        if row.get('status')!='COMPLETE':continue
        cid=str(row.get('public_path','')).replace('\\','/').rsplit('/',1)[-1].removesuffix('.npz')
        if not re.fullmatch('r14_new_\\d{4}',cid):continue
        group=row.get('operator_group','')
        if not re.fullmatch('operator_group_\\d+',group):raise ValueError('Unsafe geometry group')
        geometry=folder/group/'SPACES.npz';metadata=folder/group/'METADATA.json'
        if not geometry.is_file() or not metadata.is_file():raise FileNotFoundError('Missing real geometry cache '+group)
        if read(metadata).get('scalar_domain')!='real_material':raise ValueError('Material perturbations require real-material geometry')
        if cid in mapping:raise ValueError('Duplicate geometry execution for '+cid)
        mapping[cid]=dict(path=str(geometry),metadata=str(metadata),execution=str(p),operator_group=group,source='Existing real-material rank characterization exact-G/S equality grouping; no new rank selection or GT')
    return mapping


def fixed_tasks(root):
    geometry=geometry_map(root);cases=[f'r14_new_{i:04d}' for i in (9,20,21,23)];sources={c:source_checkpoint(root,c) for c in cases};tasks=[]
    specifications=[(c,m,'base',1) for c in cases for m in METHODS]
    specifications += [(f'r14_new_{i:04d}',m,band,sign) for i in (21,23) for m in ('twofold_hop','twofold_dr','twofold_reset') for band in ('U11','U01','U00') for sign in (-1,1)]
    for cid,method,band,sign in specifications:
        ident=cid+'__'+method+'__'+band+(str(sign) if band!='base' else '')
        source=sources[cid];geo=geometry.get(cid);public=root/'data_public/independent24'/(cid+'.npz')
        reasons=[]
        if source['path'] is None:reasons.append('NO_COMPLETE_SHARED_SOURCE_CHECKPOINT')
        if geo is None:reasons.append('NO_FROZEN_REAL_GEOMETRY_MAPPING')
        if not public.is_file():reasons.append('PUBLIC_INPUT_MISSING')
        tasks.append(dict(task_id=ident,case_id=cid,method=method,band=band,sign=sign,source_checkpoint=source,geometry=geo,public=str(public),preflight_failures=reasons,CPU_cap=60.,requested_updates=64 if band=='base' else 16))
    if len(tasks)!=56 or len({x['task_id'] for x in tasks})!=56:raise ValueError('Fixed20+36 task matrix mismatch')
    return tasks


def require_main_terminal(root):
    path=root/'queues'/MAIN/'RECEIPT.json'
    if not path.is_file():raise RuntimeError('Main terminal receipt absent; run this finite queue after main completes, no waiting service')
    q=read(path);manifest=read(root/'configs'/(MAIN+'.json'));expected={t['task_id'] for t in manifest['tasks']}
    if len(expected)!=264 or q.get('status') not in ('COMPLETED','COMPLETED_WITH_ERRORS') or q.get('running') or q.get('pending_tasks',0) or q.get('pending_cases',0):raise RuntimeError('Main264 queue is not terminal')
    if {j['task_id'] for j in q.get('jobs',[])}!=expected:raise RuntimeError('Main terminal receipt must account for all264 IDs')
    return path,q


def fresh_counts(cumulative,inherited):
    if not isinstance(cumulative,dict) or not isinstance(inherited,dict):return None
    return {k:cumulative.get(k,0)-inherited.get(k,0) for k in set(cumulative)|set(inherited) if isinstance(cumulative.get(k,0),(int,float)) and isinstance(inherited.get(k,0),(int,float))}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--legacy-repo',type=Path,required=True);a=p.parse_args();root=a.root.resolve();legacy=a.legacy_repo.resolve()
    start=time.perf_counter();cpu=0. # Absolute controller process CPU includes imports and checkpoint/geometry preparation.
    main_path,main_receipt=require_main_terminal(root)
    runner=root/'src/intervene_fullstate.py';cpu_source=legacy/'research/research_max_struct_r11/experiments/twod_queue.py'
    if not runner.is_file() or not cpu_source.is_file():raise FileNotFoundError('Actual intervention child/R11 external CPU helper required')
    q=root/'queues'/CAMPAIGN
    if q.exists():raise FileExistsError('Immutable intervention queue already exists: '+str(q))
    tasks=fixed_tasks(root)
    for task in tasks:
        destination=root/'RESULTS'/CAMPAIGN/task['task_id']
        if destination.exists():raise FileExistsError('Existing intervention never overwritten: '+str(destination))
    import psutil
    controller=psutil.Process()
    if not hasattr(controller,'cpu_affinity') or not set(CORES)<=set(controller.cpu_affinity()):raise RuntimeError('Pinned cores0,2,4,6 unavailable')
    if os.name=='nt':controller.nice(psutil.NORMAL_PRIORITY_CLASS)
    spec=importlib.util.spec_from_file_location('r15_intervention_external_cpu',cpu_source);old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
    q.mkdir(parents=True);(q/'logs').mkdir();(q/'TERMINAL_JOBS').mkdir()
    write_new(q/'FROZEN_TASKS.json',dict(campaign=CAMPAIGN,tasks=tasks,requested_tasks=56,main_receipt=str(main_path),main_status=main_receipt['status'],source_policy='Shared main twofold_hop checkpoint32, earliest complete numeric fallback recorded',truth_access='NONE',automatic_retry=False))
    write_new(q/'MAIN_TERMINAL_SNAPSHOT.json',main_receipt)
    env=os.environ.copy();env.update({k:'1' for k in THREAD_KEYS});env['PYTHONUNBUFFERED']='1';env['PYTHONPATH']=os.pathsep.join([str(legacy),str(legacy/'vendor')])
    pending=deque(tasks);running={};jobs=[];seq=0;last=None;events=(q/'EVENTS.jsonl').open('x')
    def emit(row):events.write(json.dumps(dict(utc=utc(),**row))+'\n');events.flush();print(json.dumps(row),flush=True)
    def save(status):
        nonlocal seq,last
        signature=(status,len(pending),len(jobs),tuple((p.pid,j['task']['task_id']) for p,j in running.items()))
        if signature==last:return
        last=signature;seq+=1;values=[j['external'].get('child_process_CPU_seconds') for j in jobs]
        write_new(q/('RECEIPT.json' if status!='RUNNING' else f'STATUS_{seq:06d}.json'),dict(status=status,campaign=CAMPAIGN,requested_tasks=56,pending_tasks=len(pending),running=[dict(pid=p.pid,task_id=j['task']['task_id'],cpu=j['cpu']) for p,j in running.items()],jobs=jobs,
            queue_wall_seconds=time.perf_counter()-start,controller_CPU_seconds=time.process_time()-cpu,selected_cores=list(CORES),threads=1,priority='NORMAL',automatic_retry=False,
            finished_child_CPU_seconds=sum(values) if all(x is not None for x in values) else None,measured_finished_child_CPU_seconds=sum(x for x in values if x is not None),
            missing_external_CPU_tasks=[j['task_id'] for j in jobs if j['external'].get('child_process_CPU_seconds') is None],external_cpu_helper=str(cpu_source),
            source_main_receipt=str(main_path),truth_access='NONE',cost_scope='Whole child OS CPU plus controller are disjoint; checkpoint inherited counters are not newly paid; internal child clock is nested, not added',wall_scope='This queue interval only; do not add concurrent child wall times'))
    def terminal(row):write_new(q/'TERMINAL_JOBS'/(row['task_id']+'.json'),row);jobs.append(row);emit(dict(event='END',**row))
    while pending or running:
        occupied={j['cpu'] for j in running.values()}
        for core in CORES:
            if core in occupied or not pending:continue
            task=pending.popleft();base=dict(task_id=task['task_id'],case_id=task['case_id'],method=task['method'],band=task['band'],sign=task['sign'],cpu=core,source_checkpoint=task['source_checkpoint'],geometry=task['geometry'])
            logpath=q/'logs'/(task['task_id']+'.log');log=logpath.open('x')
            if task['preflight_failures']:
                log.write(json.dumps(task['preflight_failures']));log.close()
                terminal(dict(base,status='PREFLIGHT_FAILED',returncode=None,external=dict(child_process_CPU_seconds=0.,CPU_scope='No child created; preparation charged to controller'),wall_seconds=0.,child_status=None,error=task['preflight_failures'],log=str(logpath)));continue
            argv=[str(runner),'--root',str(root),'--legacy-repo',str(legacy),'--source-checkpoint',task['source_checkpoint']['path'],'--public',task['public'],'--geometry',task['geometry']['path'],'--case-id',task['case_id'],'--method',task['method'],'--band',task['band'],'--sign',str(task['sign']),'--cpu',str(core)]
            boot=('import os,sys,runpy,psutil;'+f'psutil.Process().cpu_affinity([{core}]);'+'psutil.Process().nice(psutil.NORMAL_PRIORITY_CLASS) if os.name=="nt" else None;'+f'sys.argv={argv!r};runpy.run_path({str(runner)!r},run_name="__main__")')
            try:process=subprocess.Popen([sys.executable,'-c',boot],env=env,cwd=root,stdout=log,stderr=subprocess.STDOUT)
            except Exception as exc:
                log.write(traceback.format_exc());log.close();terminal(dict(base,status='LAUNCH_FAILED',returncode=None,external=dict(child_process_CPU_seconds=0.,CPU_scope='No child created'),wall_seconds=0.,child_status=None,error=repr(exc),log=str(logpath)));continue
            running[process]=dict(task=task,cpu=core,log=log,start=time.perf_counter(),logpath=logpath,base=base)
            emit(dict(event='START',task_id=task['task_id'],cpu=core,pid=process.pid))
        for process,j in list(running.items()):
            code=process.poll()
            if code is None:continue
            try:external=old.external_cpu(process)
            except Exception as exc:external=dict(child_process_CPU_seconds=None,CPU_scope='EXTERNAL_CPU_CAPTURE_FAILED',error=repr(exc))
            j['log'].close();rp=root/'RESULTS'/CAMPAIGN/j['task']['task_id']/'RECEIPT.json';child={};error=None
            try:child=read(rp) if rp.is_file() else {}
            except Exception as exc:error=repr(exc)
            good=code==0 and child.get('status') in ('COMPLETE','SCHEDULE_COMPLETE','WHOLE_PROCESS_CPU') and (os.name!='nt' or external.get('child_process_CPU_seconds') is not None)
            terminal(dict(j['base'],pid=process.pid,status='FINISHED' if good else 'FAILED',returncode=code,external=external,wall_seconds=time.perf_counter()-j['start'],child_status=child.get('status'),receipt=str(rp),receipt_read_error=error,
                child_internal_CPU_seconds=child.get('whole_process_CPU_seconds'),physical_cost=child.get('physics_cost'),auxiliary_inherited=child.get('auxiliary_inherited'),auxiliary_cumulative=child.get('auxiliary_cumulative'),
                auxiliary_new=fresh_counts(child.get('auxiliary_cumulative'),child.get('auxiliary_inherited')),fixed_projection_inherited=child.get('fixed_projection_inherited'),fixed_projection_cumulative=child.get('fixed_projection_cumulative'),
                fixed_projection_new=fresh_counts(child.get('fixed_projection_cumulative'),child.get('fixed_projection_inherited')),cost_scope='Only cumulative-minus-inherited is new auxiliary/background work; whole child CPU is actual fresh cost; physical counters are fresh child model counts',log=str(j['logpath'])))
            del running[process]
        save('RUNNING')
        if running:time.sleep(.5)
    final='COMPLETED' if all(j['status']=='FINISHED' for j in jobs) else 'COMPLETED_WITH_ERRORS';save(final);events.close()
    return 0 if final=='COMPLETED' else 1

if __name__=='__main__':sys.exit(main())
