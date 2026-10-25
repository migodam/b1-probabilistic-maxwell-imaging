"""Known-base E1 intrinsic compensation negative control; no work on import."""
from pathlib import Path
from dataclasses import replace
import argparse
import gc
import json
import time

import numpy as np
from scipy.optimize import least_squares

from wslam.accounting import BudgetBook, BudgetExceeded, json_default
from wslam.forward import Evaluator
from wslam.physics import acquisition, pack, packj
from wslam.world import World

COLS=np.array([2,3,4,6])
INITIAL=np.array([.2,2.,3.,.12])
LOWER=np.array([-1.,1.6,2.,.08])
UPPER=np.array([1.,2.6,5.,.18])
SCALE=np.array([2.,1.,3.,.1])


def write_json(path,value):
    temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value,indent=2,default=json_default,allow_nan=False)+'\n')
    temporary.replace(path)


def save_npz(path,arrays):
    temporary=path.with_suffix('.tmp.npz')
    np.savez_compressed(temporary,**arrays);temporary.replace(path)


def changed_world(base,direction,step,x):
    epsilon=np.asarray(base.epsinf).copy()
    epsilon[:,0]+=step*direction[0];epsilon[:,1]+=step*direction[1]
    return replace(base,epsinf=epsilon,pose=np.asarray(x[:2]).copy(),contrast=float(x[2]),radius=float(x[3]),target_loss=.02)


def profile(base,direction,step,target,tx,rx,book,directory):
    directory.mkdir()
    solver=Evaluator(book,'E0E1',tx=tx,rx=rx,h=.17,n=256,grid_refinement=2)
    cache={};trials=[];paid=[];before=book.totals('E0E1');began=time.perf_counter()
    def persist():
        arrays={}
        for i,row in enumerate(trials):
            for key,value in row.items():arrays[f'trial_{i:03d}_{key}']=value
        for i,row in enumerate(paid):
            for key,value in row.items():arrays[f'paid_{i:03d}_{key}']=value
        save_npz(directory/'ALL_TRIALS.npz',arrays)
    def get(x,want_j=False):
        key=tuple(float(v) for v in x)
        item=cache.get(key)
        if item is None or want_j and item['J'] is None:
            world=changed_world(base,direction,step,x)
            ev=solver.evaluate(world,1.,jacobian=want_j,kind='left_right')
            data=np.asarray(ev.data).copy()
            jacobian=None if not want_j else np.asarray(ev.jacobian)[:,COLS].copy()
            if data.shape!=(2,1) or want_j and jacobian.shape!=(2,4):
                raise ValueError('Unexpected minimal two-Tx / Rx8 data or Jacobian shape')
            paid_row=dict(x=np.asarray(x).copy(),data=data.copy(),wall_epsinf=np.asarray(world.epsinf),
                metadata_json=np.array(json.dumps(ev.metadata,default=json_default)))
            if jacobian is not None:paid_row['complex_J']=jacobian.copy()
            paid.append(paid_row)
            if item is None:
                item=dict(data=data,J=None);cache[key]=item
            if jacobian is not None:item['J']=jacobian
            persist()
        return item
    def residual(x):
        row=dict(kind=np.array('residual'),x=np.asarray(x).copy(),status=np.array('ATTEMPTED'))
        trials.append(row);persist()
        try:
            item=get(x)
            value=pack(item['data']-target)
            row.update(data=item['data'].copy(),residual=value.copy(),status=np.array('RETURNED'))
            return value
        except Exception as exc:
            row.update(status=np.array('FAILED'),error=np.array(type(exc).__name__+': '+str(exc)));raise
        finally:persist()
    def jacobian(x):
        row=dict(kind=np.array('jacobian'),x=np.asarray(x).copy(),status=np.array('ATTEMPTED'))
        trials.append(row);persist()
        try:
            item=get(x,True);value=packj(item['J'])
            row.update(data=item['data'].copy(),complex_J=item['J'].copy(),white_J=value.copy(),status=np.array('RETURNED'))
            return value
        except Exception as exc:
            row.update(status=np.array('FAILED'),error=np.array(type(exc).__name__+': '+str(exc)));raise
        finally:persist()
    result=dict(wall_step=step,status='NOT_RUN',success=False,final=None,
        scope='Known local noiseless finite-compensation control; not online reconstruction')
    try:
        solution=least_squares(residual,INITIAL.copy(),jac=jacobian,method='trf',
            bounds=(LOWER,UPPER),x_scale=SCALE,max_nfev=8)
        # Existing cached final forward data: no unbudgeted endpoint reevaluation.
        x=solution.x.copy();item=cache[tuple(float(v) for v in x)]
        white=pack(item['data']-target)
        result.update(status='CONVERGED' if solution.success else 'MAX_NFEV_OR_NONCONVERGENCE',
            success=bool(solution.success),optimizer_status=int(solution.status),message=solution.message,
            nfev=int(solution.nfev),njev=int(solution.njev or 0),primary_residual_white=float(np.linalg.norm(white)),
            primary_residual_energy=float(white@white),final=x,
            final_wall_epsinf=changed_world(base,direction,step,x).epsinf,
            finite_nonconvergence_retained=not bool(solution.success))
        arrays=dict(x=x,primary_data=item['data'],reference=target,primary_residual=white,
            wall_epsinf=result['final_wall_epsinf'],target_loss_fixed=np.array(.02))
        if item['J'] is not None:arrays['final_complex_J']=item['J']
        save_npz(directory/'FINAL_PRIMARY.npz',arrays)
    except BudgetExceeded as exc:
        result.update(status='STOPPED_BUDGET',error=str(exc));write_json(directory/'RESULT.json',result);raise
    except Exception as exc:
        result.update(status='FAILED',error_type=type(exc).__name__,error=str(exc),
            completed_residuals=sum(r['kind'].item()=='residual' and r['status'].item()=='RETURNED' for r in trials))
    finally:
        persist();solver.clear();gc.collect()
        result.update(stage_cost_diff={k:book.totals('E0E1')[k]-before[k] for k in before},wall_seconds=time.perf_counter()-began)
        write_json(directory/'RESULT.json',result)
    return result


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1]);args=parser.parse_args();root=args.root.resolve()
    supplement=json.loads((root/'E1_INTRINSIC_NEGATIVE_SUPPLEMENT.json').read_text())
    if json.loads((root/'results/e0/gate.json').read_text()).get('decision')!='PASS':raise RuntimeError('E0 PASS required')
    dynamic=root/'results/e1/dynamic';rank=root/'results/e1/rank_audit'
    saved=json.loads((dynamic/'RESULTS.json').read_text())['rows']
    directions=[r['wall_direction'] for r in saved if r.get('case')=='profile_direction_frozen']
    if len(directions)!=1:raise ValueError('One existing frozen wall direction required')
    direction=np.asarray(directions[0],float)
    if not np.array_equal(direction,np.asarray(supplement['frozen_wall_direction'])):
        raise ValueError('Preregistered direction differs from saved E1 mechanism direction')
    with np.load(dynamic/'pose0_f1.npz',allow_pickle=False) as raw:
        primary_reference=raw['data'][:2,[8]].copy()
        j0=raw['J'].reshape(4,16,7)[:2,8,:][:,COLS].copy()
    with np.load(rank/'JACOBIANS.npz',allow_pickle=False) as raw:
        independent_reference=raw['independent_data'][:2,[8]].copy()
        j1=raw['independent'].reshape(4,16,7)[:2,8,:][:,COLS].copy()
    b0,b1=packj(j0),packj(j1)
    diagnostic=dict(B_primary_singular_values=np.linalg.svd(b0,compute_uv=False),
        B_independent_singular_values=np.linalg.svd(b1,compute_uv=False),
        B_operator_grid_discrepancy=float(np.linalg.norm(b1-b0,2)),
        base_mean_grid_difference_white=float(np.linalg.norm(pack(independent_reference-primary_reference))),
        source='Saved rank_audit independent Jacobian and dynamic primary, no new solves',
        interpretation='Empirical discrete uncertainty only; no rank gate or continuum certificate',
        realification='All two-port real values then all imaginary values; sigma .01',shape=list(b0.shape))
    out=root/'results/e1/intrinsic_negative';out.mkdir(parents=True,exist_ok=False)
    book=BudgetBook(root/'results/budget_ledger.json')
    report=dict(status='LOCAL_RANK_REPORTED_BEFORE_PHYSICAL_ENDPOINTS',supplement=supplement,
        local_rank_diagnostic=diagnostic,profiles=[dict(wall_step=s,status='NOT_RUN') for s in (-.01,.01)],
        general_unknown_wall_reconstruction=False,exact_PDE_equality_claim=False,posterior_claim=False,
        overall_gate='NOT_ASSIGNED_PARENT_REVIEW',budget_start=book.totals('E0E1'))
    write_json(out/'RESULTS.json',report)
    print(json.dumps(dict(case='local_4x4_B_before_profiles',**diagnostic),default=json_default),flush=True)
    save_npz(out/'LOCAL_B_AUDIT.npz',dict(B_primary=b0,B_independent=b1,J_primary=j0,J_independent=j1,
        primary_reference=primary_reference,independent_reference=independent_reference,wall_direction=direction))
    remaining=book.snapshot('E0E1')['E0E1']['remaining']
    if remaining['rhs']<52 or remaining['factorizations']<18 or remaining['process_cpu_seconds']<=0:
        report.update(status='NOT_RUN',reason='INSUFFICIENT_EXISTING_E0E1_HEADROOM_NO_BORROWING')
        write_json(out/'RESULTS.json',report);return
    tx,rx=acquisition();tx=tx[:2];rx=rx[[8]]
    base=World(np.ones((2,2),bool),np.array([[3.,4.],[3.,4.]]),.08,(.2,2.),radius=.12,contrast=3.,target_loss=.02)
    stopped=False
    for index,step in enumerate((-.01,.01)):
        directory=out/f'wall_step_{step:+.2f}'
        if stopped:report['profiles'][index]['reason']='SHARED_E0E1_BUDGET_STOP';continue
        try:report['profiles'][index]=profile(base,direction,step,primary_reference,tx,rx,book,directory)
        except BudgetExceeded as exc:
            stopped=True;report['profiles'][index].update(status='STOPPED_BUDGET',error=str(exc))
        write_json(out/'RESULTS.json',report)
    # Forward-only independent endpoints after both primary profiles are durable.
    for row in report['profiles']:
        row['independent_status']='NOT_RUN'
        if stopped:
            row['independent_reason']='SHARED_E0E1_BUDGET_STOP';continue
        if row.get('final') is None:
            row['independent_reason']='NO_FINAL_PRIMARY_WORLD';continue
        before=book.totals('E0E1');solver=None
        try:
            step=row['wall_step'];x=np.asarray(row['final']);world=changed_world(base,direction,step,x)
            solver=Evaluator(book,'E0E1',tx=tx,rx=rx,h=.17,n=384,grid_refinement=3)
            ev=solver.evaluate(world,1.,jacobian=False)
            fine=np.asarray(ev.data).copy();metadata=dict(ev.metadata)
            directory=out/f'wall_step_{step:+.2f}'
            with np.load(directory/'FINAL_PRIMARY.npz',allow_pickle=False) as raw:primary=raw['primary_data'].copy()
            fine_residual=pack(fine-independent_reference)
            endpoint_difference=pack(fine-primary)
            delta_endpoint=float(np.linalg.norm(endpoint_difference))
            delta_base=diagnostic['base_mean_grid_difference_white']
            row.update(independent_status='SCORED',independent_residual_white=float(np.linalg.norm(fine_residual)),
                independent_residual_energy=float(fine_residual@fine_residual),
                endpoint_mean_grid_difference_white=delta_endpoint,base_mean_grid_difference_white=delta_base,
                observed_endpoint_plus_base_grid_discrepancy=delta_endpoint+delta_base,
                uncertainty_interpretation='Observed two-grid mean discrepancy sum, not a rigorous continuous/PDE error bound',
                independent_metadata=metadata)
            save_npz(directory/'FINAL_INDEPENDENT.npz',dict(x=x,primary=primary,independent=fine,
                independent_base=independent_reference,independent_residual=fine_residual,
                endpoint_grid_difference_white=endpoint_difference,wall_epsinf=world.epsinf,
                metadata_json=np.array(json.dumps(metadata,default=json_default))))
        except BudgetExceeded as exc:stopped=True;row.update(independent_status='STOPPED_BUDGET',error=str(exc))
        except Exception as exc:row.update(independent_status='FAILED',error_type=type(exc).__name__,error=str(exc))
        finally:
            if solver is not None:solver.clear()
            gc.collect();row['independent_stage_cost']={k:book.totals('E0E1')[k]-before[k] for k in before}
            write_json(out/'RESULTS.json',report)
    report.update(status='STOPPED_E0E1_BUDGET' if stopped else 'DIAGNOSTIC_COMPLETE_PARENT_REVIEW_REQUIRED',budget_end=book.totals('E0E1'))
    write_json(out/'RESULTS.json',report)
    write_json(out/'gate.json',dict(decision='INCONCLUSIVE',execution_status=report['status'],
        reason='Finite compensation residuals and empirical independent-grid uncertainty only; parent scientific review required',exact_PDE_equality_claim=False))


if __name__=='__main__':main()
