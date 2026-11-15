"""Finite two-child Apple CPU queue for the frozen shared-readout control."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
import argparse
from concurrent.futures import ThreadPoolExecutor,as_completed
from datetime import datetime,timezone
import json
from pathlib import Path
import platform
import resource
import shutil
import subprocess
import sys
import time


def dump(path,value):Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True)
    p.add_argument('--legacy-repo',type=Path,required=True);p.add_argument('--manifest',type=Path,required=True)
    p.add_argument('--workers',type=int,choices=(1,2),default=2);a=p.parse_args()
    cpu,wall=time.process_time(),time.perf_counter();root=a.root.resolve();repo=a.legacy_repo.resolve();manifest=json.loads(a.manifest.read_text())
    campaign=manifest['campaign'];tasks=manifest['tasks'];queue=root/'queues'/campaign/'attempt_02';queue.mkdir(parents=True,exist_ok=False)
    (queue/'logs').mkdir();frozen=queue/'FROZEN_MANIFEST.json';dump(frozen,manifest)
    for task in tasks:
        if (root/'RESULTS'/campaign/task['task_id']).exists():raise FileExistsError('Existing results retained: '+task['task_id'])
    audit=root/'audit'/campaign;source=audit/'source_snapshot_attempt_02';source.mkdir(exist_ok=False)
    paths=[root/'src/shared_readout_control.py',root/'src/run_shared_material.py',root/'src/run_shared_material_queue.py',
        root/'src/run_experiment.py',root/'src/dualmetric_core.py',root/'RESEARCH_LOOP.md',a.manifest.resolve(),
        root/'audit/shared_readout4_v1/CONTRACT.md',
        repo/'research/research_max_struct_r13/src/full_current.py',
        repo/'research/research_max_struct_r4/src/cie_onehop.py',repo/'research/research_max_struct_r2/src/cie_global.py',
        repo/'research/research_fbg_r3_experiments/src/cie_tsom_baseline.py',
        repo/'research/research_max_struct_r8/src/benchmark_solver.py']
    provenance=[]
    for path in paths:
        destination=source/path.relative_to(repo);destination.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(path,destination)
        provenance.append(dict(source=str(path),snapshot=str(destination),bytes=destination.stat().st_size,mtime_ns=path.stat().st_mtime_ns))
    dump(source/'PROVENANCE.json',dict(status='FROZEN_BEFORE_FIRST_CHILD',utc=datetime.now(timezone.utc).isoformat(),files=provenance,
        commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip(),
        contract='same core step/current/PR/stage; private two-material-function replacement only',truth_access='NONE'))
    before=resource.getrusage(resource.RUSAGE_CHILDREN);jobs=[]
    runtime=audit/'runtime_existing_psutil'
    if not (runtime/'psutil/__init__.py').exists():raise FileNotFoundError('Recorded existing ABI3psutil privatepath is required')
    previous=json.loads((root/'queues'/campaign/'RECEIPT.json').read_text())
    dump(queue/'STARTED.json',dict(tasks=len(tasks),workers=a.workers,interpreter=sys.executable,platform=platform.platform(),
        machine=platform.machine(),BLAS_threads=1,affinity='OS scheduling; no pinning/exclusivity claim',source_snapshot=str(source),truth_access='NONE'))
    def child(task):
        log=queue/'logs'/(task['task_id']+'.log');started=time.perf_counter()
        argv=[sys.executable,str(root/'src/run_shared_material.py'),'--root',str(root),'--legacy-repo',str(repo),
            '--manifest',str(frozen),'--task-id',task['task_id']]
        env=os.environ.copy();env['PYTHONPATH']=str(runtime)+os.pathsep+str(repo)
        with log.open('x') as stream:completed=subprocess.run(argv,stdout=stream,stderr=subprocess.STDOUT,check=False,env=env)
        folder=root/'RESULTS'/campaign/task['task_id'];receipt=folder/'RECEIPT.json'
        record=json.loads(receipt.read_text()) if receipt.exists() else None
        native_path=folder/'NATIVE_HISTORY.json';native=json.loads(native_path.read_text()) if native_path.exists() else []
        shared_count=0 if record is None else record.get('auxiliary_cost',{}).get('shared_material_LBFGS_readouts',0)
        arc_count=0 if record is None else record.get('auxiliary_cost',{}).get('aplus_readouts',0)
        verified=bool(shared_count>0 and arc_count==0 and native and all(row.get('material_solver')=='shared_LBFGS60' for row in native))
        return dict(task_id=task['task_id'],case_id=task['case_id'],pair=task['pair'],returncode=completed.returncode,
            child_status=None if record is None else record['status'],
            child_reported_CPU_seconds=None if record is None else record['CPU_seconds'],wall_seconds=time.perf_counter()-started,
            shared_solver_tag_verified=verified,shared_material_readouts=shared_count,Arc_readouts=arc_count,
            receipt=str(receipt),log=str(log))
    with ThreadPoolExecutor(max_workers=a.workers) as executor:
        futures={executor.submit(child,task):task for task in tasks}
        for future in as_completed(futures):
            task=futures[future]
            try:row=future.result()
            except Exception as exc:row=dict(task_id=task['task_id'],returncode=None,child_status='CONTROLLER_FAILURE',error=repr(exc))
            jobs.append(row);dump(queue/'PROGRESS.json',dict(expected=len(tasks),terminal=len(jobs),jobs=jobs))
            print(json.dumps(row),flush=True)
    after=resource.getrusage(resource.RUSAGE_CHILDREN)
    good=lambda row:row.get('returncode')==0 and row.get('child_status')=='COMPLETE' and row.get('shared_solver_tag_verified')
    receipt=dict(status='ALL_EXPECTED_TASKS_TERMINAL',expected=len(tasks),terminal=len(jobs),successful=sum(good(row) for row in jobs),
        failed_or_missing=sum(not good(row) for row in jobs),jobs=jobs,workers=a.workers,BLAS_threads=1,
        solver_children_whole_OS_CPU_seconds=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,
        controller_CPU_seconds=time.process_time()-cpu,controller_wall_seconds=time.perf_counter()-wall,
        source_snapshot=str(source),truth_access='NONE',GT_scoring='NOT_AUTHORIZED_OR_RUN',
        timing_scope='childOS CPU includes imports/material/current/physicalcorrector/serialization/failures; controllerCPU disjoint; wall notadded',
        no_absolute_Mac_Windows_speed_comparison=True)
    receipt.update(explicit_retry_after_setup_failure=True,prior_failed_attempt_receipt=str(root/'queues'/campaign/'RECEIPT.json'),
        prior_failed_attempt_OS_CPU_seconds=previous['solver_children_whole_OS_CPU_seconds'],
        cumulative_solver_children_OS_CPU_seconds=receipt['solver_children_whole_OS_CPU_seconds']+previous['solver_children_whole_OS_CPU_seconds'],
        existing_psutil_runtime_receipt=str(audit/'PSUTIL_RUNTIME_RECEIPT.json'),runtime_PYTHONPATH=str(runtime),
        prior_attempt_physics='all16failed at psutil import before model/solver setup; zero physical LU/RHS')
    dump(queue/'RECEIPT.json',receipt)
    inventory=[]
    for task in tasks:
        folder=root/'RESULTS'/campaign/task['task_id'];files=[]
        if folder.exists():
            files=[dict(path=str(path.relative_to(root)),bytes=path.stat().st_size,mtime_ns=path.stat().st_mtime_ns)
                   for path in sorted(folder.rglob('*')) if path.is_file()]
        inventory.append(dict(task_id=task['task_id'],case_id=task['case_id'],pair=task['pair'],files=files,missing=not bool(files)))
    seal=root/'predictions'/campaign/'attempt_02';seal.mkdir(parents=True,exist_ok=False)
    dump(seal/'PREDICTION_SEAL.json',dict(status='ALL_TERMINAL_PREDICTIONS_SEALED_WITH_FAILURES_RETAINED',
        expected_tasks=len(tasks),terminal_tasks=len(jobs),successful=receipt['successful'],failed_or_missing=receipt['failed_or_missing'],
        queue_receipt=str(queue/'RECEIPT.json'),inventory=inventory,truth_access_before_seal='NONE',GT_scoring_authorized=False,
        mechanism_scope='four exposed cases only; commonLBFGS60 material solver; not a newmainalgorithm or blind recovery confirmation',
        original_main_and_independent_algorithms_modified=False,source_snapshot=str(source)))
    print(json.dumps(dict(queue_receipt=str(queue/'RECEIPT.json'),seal=str(seal/'PREDICTION_SEAL.json'),successful=receipt['successful'],failed=receipt['failed_or_missing'])))
    return 0 if receipt['failed_or_missing']==0 else 1


if __name__=='__main__':sys.exit(main())
