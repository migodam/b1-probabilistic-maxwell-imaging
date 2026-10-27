"""NEW development strong baselines: prior x0=zeros66 -> Original Probit64/128.

There is no Native optimizer. The observation-independent initial native point
[zeros64,.5,.5] pays8RHS; Original probit_optimize is then called unchanged and
reuses that Point at x0. Caller-owned Target/central Book pays all physical work.
Neither is Original Native96+Probit64 or a locked-D arm. Root alone executes.
"""
from __future__ import annotations
import os
for _name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_name]='1'
import argparse,gc,hashlib,json,sys,time
from dataclasses import asdict
from pathlib import Path
import numpy as np
OUT=Path(__file__).resolve().parents[1];TRACK=OUT/'track_c';sys.path.insert(0,str(OUT))
from shared.physics import Book,Target,probit_optimize,write
from track_c.anchors import AnchorConfig,_ObservedTarget,_Trace,_accepted_iteration_hook,build_online_context

CAPS=(64,128)


def run(target,y,parent,destination,cap=64):
    """Return standard observation/anchor context and explicit NEW-baseline report."""
    y=np.asarray(y)
    if y.shape!=(32,8) or not np.iscomplexobj(y) or not np.isfinite(y).all():raise ValueError('Original complex32x8 observation')
    destination=Path(destination)
    if destination.exists():raise FileExistsError('prior Probit path never silently overwrites/replays')
    destination.mkdir(parents=True)
    if cap not in CAPS:raise ValueError('frozen development Probit cap64 or128 required')
    name=f'prior_x0_probit{cap}'
    config=AnchorConfig(name,1,cap)
    observed=_ObservedTarget(target,destination)
    trace=_Trace(observed,y,parent,destination/'path',config)
    initial_ledger=target.book.totals();started=time.perf_counter()
    try:
        q0=np.r_[np.zeros(64),[.5,.5]]
        initial=observed.evaluate(q0,purpose=f'{parent}_prior_x0_no_Native_optimizer')
        initial_chi=initial.state.chi.copy()
        with _accepted_iteration_hook(probit_optimize,trace,'probit'):
            point,x,H,statistics=probit_optimize(observed,y,str(parent),destination/'original_optimizer/probit',initial,cap=cap)
        context=build_online_context(observed,y,parent,point,x,H)
        context.update(prior_initial_q=q0,prior_initial_chi=initial_chi)
        np.savez_compressed(destination/'CONTEXT.npz',**context)
        ledger=target.book.totals()
        report=dict(parent=str(parent),config=dict(name=name,Probit_cap=cap,initial_point_fresh_state_cap=1),
            definition=f'DIRECT_PRIOR_X0_PLUS_ORIGINAL_PROBIT{cap}_NEW_ANCHOR',Native_optimizer_ran=False,
            not_original_full_baseline=True,not_locked_D_arm=True,initial_native_q=q0.tolist(),initial_Probit_x=np.zeros(66).tolist(),
            GT_read=False,original_Probit_math_unchanged=True,original_probit_report=statistics,
            total_fresh_forward_attempts=observed.attempts,total_fresh_forward_attempt_cap=1+cap,worst_RHS_per_parent=40*(1+cap),
            path_resource_cost={key:ledger[key]-initial_ledger[key] for key in ('RHS','forward_RHS','derivative_RHS','states','LU','physics_wall_seconds')},
            wall_seconds=time.perf_counter()-started,path_complete=True,new_fallback=False,
            online_context='observation,residual,q0,chi0,J plus paid prediction/x/H; no truth',output_context=str(destination/'CONTEXT.npz'))
        if observed.attempts>1+cap or report['path_resource_cost']['RHS']>40*(1+cap):raise RuntimeError('prior Probit declared path cap exceeded')
        write(destination/'ANCHOR_REPORT.json',report)
        return context,report
    except BaseException as exc:
        write(destination/'FAILURE.json',dict(error=repr(exc),parent=str(parent),GT_read=False,no_implicit_retry=True));raise
    finally:trace.close()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest',type=Path,default=OUT/'configs/parents_CDEV.json')
    parser.add_argument('--online-root',type=Path,default=OUT/'raw/online')
    parser.add_argument('--output-root',type=Path,default=TRACK/'contexts/CDEV')
    parser.add_argument('--caps',nargs='+',type=int,choices=CAPS,default=list(CAPS))
    args=parser.parse_args()
    rows=json.loads(args.manifest.read_text())['parents'];parents=[str(r['parent']) for r in rows]
    if len(parents)!=64 or any(not p.startswith('R6_CDEV_') for p in parents):raise ValueError('explicit exposed CDEV64 development only')
    if not args.output_root.resolve().is_relative_to(TRACK.resolve()):raise ValueError('Track C outputs only')
    freeze=json.loads((TRACK/'probit_prior_baseline/DEVELOPMENT_SOURCE_FREEZE.json').read_text())
    for row in freeze['files']:
        content=Path(row['path']).read_bytes()
        if len(content)!=row['bytes'] or hashlib.blake2b(content).hexdigest()!=row['blake2b']:raise RuntimeError('development baseline frozen source/protocol changed')
    pending=[(p,cap) for p in parents for cap in args.caps if not (args.output_root/p/f'prior_x0_probit{cap}'/'REPORT.json').exists()]
    auth=json.loads((OUT/'configs/AUTHORIZATION.json').read_text());ledger=Book.totals()
    available=min(auth['stage_caps']['C']-ledger['stages'].get('C',{}).get('RHS',0),auth['total_complex_RHS']-ledger['RHS'])
    reserve=sum(40*(1+cap) for _,cap in pending)
    budget=dict(paths_pending=len(pending),caps=args.caps,worst_RHS_per_cap={'64':2600,'128':5160},worst_remaining_RHS=reserve,available_RHS=available,
                status='RUN' if reserve<=available else 'HOLD',GT_read=False,not_locked_D_arm=True)
    write(TRACK/'probit_prior_baseline/PRE_RUN_BUDGET.json',budget)
    if budget['status']!='RUN':raise RuntimeError('HOLD: complete dev64 worst-case reserve unavailable')
    for index,(parent,cap) in enumerate(pending):
        name=f'prior_x0_probit{cap}'
        destination=args.output_root/parent/name
        if (destination/'REPORT.json').exists():continue
        with np.load(args.online_root/f'{parent}.npz',allow_pickle=False) as z:y,sigma=z['observation'].copy(),float(z['sigma'])
        book=Book('C',f'anchor_{parent}_{name}',parent,'DEVELOPMENT_PRIOR_PROBIT_BASELINE',max_wall=120)
        status,error='COMPLETE',None;tick=time.perf_counter()
        try:
            target=Target(book,sigma=sigma);context,report=run(target,y,parent,destination,cap)
            write(destination/'REPORT.json',dict(report,end_to_end_anchor_wall_seconds=time.perf_counter()-tick,sigma=sigma))
            target.backend()._transfers.clear();del target,context,report;gc.collect()
        except BaseException as exc:status,error='FAILED',repr(exc);raise
        finally:book.close(status,error)
        if (index+1)%16==0:print(json.dumps(dict(event='PRIOR_PROBIT_DEVELOPMENT',completed=index+1,RHS=Book.totals()['RHS'])),flush=True)
    write(TRACK/'probit_prior_baseline/ANCHORS_COMPLETE.json',dict(status='COMPLETE',n=64,parents=parents,configs=[f'prior_x0_probit{cap}' for cap in args.caps],caps=args.caps,
        GT_read=False,Native_optimizer_ran=False,not_locked_D_arm=True))


if __name__=='__main__':main()
