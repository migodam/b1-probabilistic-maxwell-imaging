#!/usr/bin/env python3
"""POST-GT low-contrast physics-regime sensitivity; solver reads anonymous Y only.

A separate generator applies the disclosed uniform rescaling to both original
materials. This solver never reads those materials, support, labels or peaks.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
import argparse,json,shutil,time,traceback
from pathlib import Path
from datetime import datetime,timezone
from threadpoolctl import threadpool_limits
from src.solver import Run,load_public,json_safe
ROOT=Path(__file__).resolve().parent
CONFIG=ROOT/'audit/FULLJ_DIAGNOSTIC_CONTRAST_CONFIG.json'
SNAPSHOT=ROOT/'audit/FULLJ_DIAGNOSTIC_CONTRAST_SOURCE_SNAPSHOT'


def write(path,data):
    path=Path(path)
    if path.exists():raise FileExistsError('Preserving '+str(path))
    path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(json_safe(data),indent=2,allow_nan=False)+'\n')


def freeze():
    if SNAPSHOT.exists():raise FileExistsError('Contrast control snapshot already exists')
    SNAPSHOT.mkdir(parents=True)
    rows=[]
    for name in ('contrast_control.py','audit/FULLJ_DIAGNOSTIC_CONTRAST_CONFIG.json','src/solver.py','src/physics.py','src/multiscale.py','src/__init__.py'):
        out=SNAPSHOT/name;out.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/name,out);rows.append(dict(path=name,bytes=out.stat().st_size))
    write(SNAPSHOT/'FREEZE.json',dict(frozen_utc=datetime.now(timezone.utc).isoformat(),files=rows,
          status='POST_GT_CONTROLLED_PHYSICS_REGIME_SENSITIVITY_NOT_MAIN_RECOVERY_OR_BLIND',
          private_material_access_by_solver=False,comparison_limit='new low-contrast observations cannot replace original high-contrast main cases',hashes='NOT_COMPUTED'))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--freeze',action='store_true');p.add_argument('--public');p.add_argument('--operators');p.add_argument('--output')
    a=p.parse_args()
    if a.freeze:freeze();return
    if not all((a.public,a.operators,a.output)):p.error('Supply public/operators/output')
    manifest=json.loads((SNAPSHOT/'FREEZE.json').read_text())
    for row in manifest['files']:
        if (ROOT/row['path']).read_bytes()!=(SNAPSHOT/row['path']).read_bytes():raise RuntimeError('Source changed aftercontrastfreeze')
    contract=json.loads(CONFIG.read_text());live=None;begun,cpu=time.perf_counter(),time.process_time()
    try:
        with threadpool_limits(limits=1):
            live=Run(load_public(a.public,a.operators),'fullj_native',contract['config']);live.started,live.cpu_started=begun,cpu
            try:result=live.run()
            except BaseException as exc:
                live.status='FAILED_PAID_UNEXPECTED';live.stop_reason=type(exc).__name__;live.failures[type(exc).__name__]+=1
                live.events.append(dict(kind='UNEXPECTED_FAILURE_PAID',error=repr(exc),seconds=live.elapsed()));result=live.finish()
    except BaseException as exc:
        result=dict(method='contrast_native',status='FAILED_INITIALIZATION',error=repr(exc),traceback=traceback.format_exc(),
              wall_seconds=time.perf_counter()-begun,cpu_seconds=time.process_time()-cpu,accounting={'counts':None},top1=None)
    result['method']='contrast_native';result['development_status']='POST_GT_CONTROLLED_PHYSICS_REGIME_SENSITIVITY_NOT_MAIN_RECOVERY_OR_BLIND'
    result['truth_access']='SOLVER_NONE; anonymous rescaled observations plus original public operators only'
    result['control_contract']=contract;result['source_freeze_utc']=manifest['frozen_utc'];result['public_input_path']=str(Path(a.public).resolve())
    result['comparison_limit']='Low contrast changes forward data and inverse regime; not a replacement score for original8-18 material study'
    write(a.output,result)
    print(json.dumps(dict(output=a.output,status=result['status'],wall_seconds=result['wall_seconds'],
         residual=(result.get('top1') or {}).get('relative_residual'),accounting=result.get('accounting'))),flush=True)


if __name__=='__main__':main()
