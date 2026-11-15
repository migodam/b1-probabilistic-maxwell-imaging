"""Same-process private G19 factory injection into the complete frozen C2 main."""
from __future__ import annotations
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
import argparse,importlib.util,json,runpy,sys,time,types
from datetime import datetime,timezone
from pathlib import Path


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--root',type=Path,required=True);parser.add_argument('--legacy-repo',type=Path,required=True);parser.add_argument('--manifest',type=Path,required=True);parser.add_argument('--task-id',required=True);parser.add_argument('--cpu',type=int)
    a=parser.parse_args();s=a.root.resolve();W=a.legacy_repo.resolve();R=s.parent;wall=time.perf_counter()
    cfg=json.loads(a.manifest.read_text());matches=[t for t in cfg['tasks'] if t['task_id']==a.task_id]
    if len(matches)!=1:raise ValueError('Exactly one fixed task is required')
    task=matches[0];arm=task['method']
    if task.get('algorithm')!='g19' or task.get('basis_rank')!=19 or cfg.get('exposure_status')!='EXPOSED24_FOLLOWUP':raise ValueError('Original g19 dispatch and fixed19 exposed followup contract required')
    if float(cfg['protocol']['CPU_cap'])!=120 or float(cfg['protocol']['proposal_fraction'])!=.6:raise ValueError('Original total120/proposal72 budget required')
    out=s/'RESULTS'/task['campaign']/task['task_id']
    if out.exists():raise FileExistsError(out)
    runner=R/'runtime/equal_budget_C2_runner.py';frozen=R/'runtime/frozen_C2_runner.py'
    if not runner.is_file() or not frozen.is_file():raise FileNotFoundError('Original R14 runtime C2 sources required')
    expected=frozen.read_text().replace('deadline=min(cap,time.process_time()+cap*(1-proposal_fraction))','deadline=cap')
    if runner.read_text()!=expected:raise ValueError('C2 source differs from its original one-line absolute-budget extension')
    paths=[s/'src',s/'runtime',R/'runtime',W/'research/research_max_struct_r13/runtime',W/'research/research_max_struct_r13/src',W/'research/research_max_struct_r12/src',W/'research/research_max_struct_r11/src',W/'research/research_max_struct_r10/src',W,W/'vendor']
    sys.path[:0]=[str(p) for p in paths if p.is_dir()]
    # Fail on real dependency absence; never substitute virtual path/module stubs.
    dependencies={}
    for module in ('field_current','material_events'):
        spec=importlib.util.find_spec(module)
        if spec is None:raise ModuleNotFoundError('Required inherited dependency '+module+' absent from declared paths')
        dependencies[module]=spec.origin
    import numpy as np
    from observable19_factory import ARMS,bind_runner
    from research.research_max_struct_r13.src import full_current as original
    if arm not in ARMS:raise ValueError('Unregistered fixed arm '+arm)
    input_path=(s/task['input_file']).resolve()
    if input_path.parent!=(R/'data_public/independent24').resolve():raise ValueError('Only the original public independent24 NPZ is permitted')
    with np.load(input_path,allow_pickle=False) as z:sigma=float(z['sigma']);anchor_m0=z['m_init'].copy()
    run,context=bind_runner(arm,sigma,anchor_m0);private=types.ModuleType('full_current');private.run_full_current=run
    original_prepare=original.prepare_full;original_material=original.full_material_readout
    previous_module=sys.modules.get('full_current');previous_argv=sys.argv
    sys.modules['full_current']=private
    command=[str(runner),'--root',str(s),'--legacy-repo',str(W),'--manifest',str(a.manifest.resolve()),'--task-id',a.task_id]
    if a.cpu is not None:command.extend(['--cpu',str(a.cpu)])
    record=dict(created_utc=datetime.now(timezone.utc).isoformat(),task_id=a.task_id,arm=arm,dispatch='Complete original C2 via runpy same process; private import slot restored in finally; no CPU clock restart',runner=str(runner),runner_bytes=runner.stat().st_size,C2_declared_budget_only_source_match=True,declared_paths=[dict(path=str(p),exists=p.is_dir()) for p in paths],resolved_dependencies=dependencies,original_full_current=original.__file__,input=str(input_path),input_keys_read_by_wrapper=['sigma','m_init'],sigma=sigma,wrapper_process_CPU_before_C2=time.process_time(),truth_access='NONE',protocol_deepcopy_unchanged=True)
    audit=s/'audit/dispatch_provenance'/task['campaign'];audit.mkdir(parents=True,exist_ok=True)
    with (audit/(a.task_id+'.json')).open('x') as f:json.dump(record,f,indent=2)
    try:
        sys.argv=command;runpy.run_path(str(runner),run_name='__main__')
    finally:
        sys.argv=previous_argv
        if previous_module is None:sys.modules.pop('full_current',None)
        else:sys.modules['full_current']=previous_module
        record.update(whole_process_CPU_after_C2=time.process_time(),wrapper_wall_seconds=time.perf_counter()-wall,original_package_globals_unchanged=original.prepare_full is original_prepare and original.full_material_readout is original_material,actual_setup=context['setup'],actual_readout_calls=context['readout_calls'])
        if out.is_dir():
            with (out/'WRAPPER_PROVENANCE.json').open('x') as f:json.dump(record,f,indent=2,allow_nan=False)

if __name__=='__main__':main()
