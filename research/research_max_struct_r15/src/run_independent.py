"""R15 fixed dispatcher: stateful parent solver or exact inherited C2, same process."""
from __future__ import annotations
import os
for _k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[_k]='1'
import argparse,importlib.util,json,runpy,sys,time
from datetime import datetime,timezone
from pathlib import Path
STATEFUL=('g19','g19_s19','twofold','twofold_hop','twofold_dr','twofold_ap','twofold_full_dr','twofold_reset')
BASELINES=('cie','dual_cie','dr')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('root','legacy-repo','manifest'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--task-id',required=True);p.add_argument('--cpu',type=int);a=p.parse_args();root=a.root.resolve();legacy=a.legacy_repo.resolve()
    cfg=json.loads(a.manifest.read_text());matches=[t for t in cfg['tasks'] if t['task_id']==a.task_id]
    if len(matches)!=1:raise ValueError('Exactly one registered task required')
    t=matches[0];method=t['method']
    if method not in STATEFUL+BASELINES or t.get('algorithm')!=method:raise ValueError('Unregistered method/algorithm')
    for key in ('campaign','task_id'):
        x=t[key]
        if not x or Path(x).name!=x or x in ('.','..') or any(c in x for c in '/\\'):raise ValueError('Unsafe task identity')
    if t['input_file']!='data_public/independent24/'+t['case_id']+'.npz':raise ValueError('Public input must use parent-copied R15 path')
    if float(cfg['protocol']['CPU_cap'])!=120 or float(cfg['protocol']['proposal_fraction'])!=.6:raise ValueError('Frozen120/72 budget required')
    output=root/'RESULTS'/t['campaign']/t['task_id']
    if output.exists():raise FileExistsError('Frozen predictions never overwritten: '+str(output))
    baseline=method in BASELINES;runner=root/('runtime/equal_budget_C2_runner.py' if baseline else 'src/run_stateful.py')
    if not runner.is_file():raise FileNotFoundError(runner)
    reference=root/'runtime/frozen_C2_runner.py';source_match=None
    if baseline:
        expected=reference.read_text().replace('deadline=min(cap,time.process_time()+cap*(1-proposal_fraction))','deadline=cap')
        source_match=runner.read_text()==expected
        if not source_match:raise ValueError('Inherited C2 differs from original single budget-line extension')
    record=dict(created_utc=datetime.now(timezone.utc).isoformat(),task_id=t['task_id'],case_id=t['case_id'],campaign=t['campaign'],method=method,
        runner_family='FROZEN_NATIVE_C2' if baseline else 'PARENT_STATEFUL_G19_TWOFOLD',runner=str(runner),runner_bytes=runner.stat().st_size,
        manifest=str(a.manifest.resolve()),input_file=t['input_file'],git_identity=cfg.get('git_identity'),
        dispatch='runpy same process, no CPU clock restart',declared_budget_only_source_match=source_match,
        whole_process_CPU_at_dispatch=time.process_time(),truth_access='NONE',portable_modules={})
    directory=root/'audit/dispatch_provenance'/t['campaign'];directory.mkdir(parents=True,exist_ok=True)
    with (directory/(t['task_id']+'.json')).open('x') as f:json.dump(record,f,indent=2)
    previous_argv=sys.argv;old_modules={}
    try:
        if baseline:
            sys.path[:0]=[str(root/'runtime'),str(root/'src'),str(legacy),str(legacy/'vendor')]
            # Actual copied sources, never empty substitute packages or a remote
            # dependency on an absent R13 package. All imports remain paid.
            modules=[('field_current','field_current.py'),('material_events','material_events.py')]
            if method in ('cie','dual_cie'):modules.append(('full_current','frozen_full_current.py'))
            for name,filename in modules:
                source=root/'runtime'/filename
                if not source.is_file():raise FileNotFoundError(source)
                old_modules[name]=sys.modules.get(name)
                spec=importlib.util.spec_from_file_location(name,source);module=importlib.util.module_from_spec(spec)
                sys.modules[name]=module;spec.loader.exec_module(module)
                record['portable_modules'][name]=dict(source=str(source),bytes=source.stat().st_size)
        sys.argv=[str(runner),*previous_argv[1:]]
        runpy.run_path(str(runner),run_name='__main__')
    finally:
        sys.argv=previous_argv
        for name,old in old_modules.items():
            if old is None:sys.modules.pop(name,None)
            else:sys.modules[name]=old
        record['whole_process_CPU_after_runner']=time.process_time()
        # Independent of prediction-directory availability: failed imports/setup
        # still receive a final dispatch receipt and external queue CPU charge.
        with (directory/(t['task_id']+'__FINAL.json')).open('x') as f:json.dump(record,f,indent=2)
        if output.is_dir():
            with (output/'WRAPPER_PROVENANCE.json').open('x') as f:json.dump(record,f,indent=2)

if __name__=='__main__':main()
