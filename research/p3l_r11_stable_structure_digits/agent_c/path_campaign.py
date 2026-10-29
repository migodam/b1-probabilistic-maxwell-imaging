"""Actual global/local residual-gated fixed Green-space R10 path replay.

Reads only frozen setup, occupancies, observation-only y, and cached direct F.
All new approximate and fallback solves are charged against 20,000 RHS.
"""
from __future__ import annotations
import os
for _n in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'):
    os.environ[_n] = '1'
import csv, json, sys, time, traceback
from pathlib import Path
import numpy as np

D=Path(__file__).resolve().parent
sys.path.insert(0,str(D.parent/'shared'))
sys.path.insert(0,str(D.parent/'agent_b'))
sys.path.insert(0,str(D))
import core
from controlled_solver import GlobalControlledSolver
from fixed_total_field import CostLedger
from dissipative_bound import GlobalDissipativeBound


class CaptureBound(GlobalDissipativeBound):
    def __init__(self,*a,**kw):
        super().__init__(*a,**kw);self.capture=[]
    def __call__(self,a,delta,residual,approxF,y=None):
        result=super().__call__(a,delta,residual,approxF,y)
        self.capture.append((approxF.copy(), residual.copy(),result.copy()))
        return result


def csv_write(path, rows):
    if not rows:return
    keys=list(dict.fromkeys(k for row in rows for k in row))
    with Path(path).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)


def main():
    started=time.perf_counter();out=D/'path_results';out.mkdir(exist_ok=True)
    if (out/'SUMMARY.json').exists():raise RuntimeError('preserve completed result; do not overwrite')
    setup=core.load_setup(); C,O,E,Fh=(setup[k] for k in ('C','O','E','Fh'))
    ledger=CostLedger(out/'COST_LEDGER.jsonl',rhs_cap=20000)
    rows=[];attempt_rows=[];shape_rows=[];solvers={};status='FAILED';error=None
    try:
        for mode in ('global','local'):
            s=GlobalControlledSolver(C,O,E,Fh,ledger=ledger,ranks=(32,64,96,128),bound_mode=mode)
            # Capture already calculated outputs for an independent cached-reference audit.
            s.bound_object=CaptureBound(C,O,core.SIGMA,mode=mode)
            s.bound=s.bound_object;solvers[mode]=s
        for pi in range(8):
            parent=f'R10BI_FRESH_{pi:03d}'
            base=core.R10/'followup_shape_reuse'
            anchor=np.load(base/f'{parent}_ANCHOR.npz')
            a0=anchor['fraction'][setup['ROI']]
            z=np.load(core.R10/'results/candidates_fresh10/inputs/online'/f'{parent}.npz')
            if set(z.files)!={'observation','sigma','parent'}:raise ValueError('observation-only input contract')
            y=z['observation']
            for ti in range(3):
                z=np.load(base/f'{parent}_t{ti}.npz');a=z['fraction'][setup['ROI']]
                new=(a0==0)&(a>0)
                shape_rows.append(dict(parent=parent,t=float(z['t']),new_support_cells=int(new.sum()),
                    new_input_energy_fraction=float(np.linalg.norm(np.sqrt(a[new])[:,None]*E[new])**2/
                        np.linalg.norm(np.sqrt(a)[:,None]*E)**2)))
                for ni,delta in enumerate(z['delta']):
                    ref=z['DIRECT_LU'][ni];ref_phi=core.phi(ref,y)
                    for mode in (('global','local') if (pi+ti+ni)%2==0 else ('local','global')):
                        s=solvers[mode];s.bound_object.capture=[]
                        before=dict(ledger.summary());result=s.evaluate(a,delta,y=y);diag=result['diagnostics']
                        after=ledger.summary()
                        recorded=s.bound_object.capture
                        for k,(fa,res,b) in enumerate(recorded):
                            actual=float(np.sqrt(2)/core.SIGMA*np.linalg.norm(fa-ref))
                            actual_pot=abs(core.phi(fa,y)-ref_phi)
                            attempt_rows.append(dict(parent=parent,t_index=ti,t=float(z['t']),node=ni,
                                mode=mode,rank=diag['attempts'][k]['rank'],white_bound=b['white_bound'],
                                potential_bound=b['potential_bound'],actual_white_error=actual,
                                actual_potential_error=actual_pot,accepted=b['accepted'],
                                bound_dominates_sampled_error=bool(b['white_bound'] is not None and
                                    actual<=b['white_bound']*(1+1e-8)+1e-10),
                                potential_bound_dominates_sampled_error=bool(b['potential_bound'] is not None and
                                    actual_pot<=b['potential_bound']*(1+1e-8)+1e-10),
                                numerical_coercivity=b['numerical_coercivity'],
                                residual_weighted_norm=b['residual_weighted_norm'],
                                solve_seconds=diag['attempts'][k]['solve_seconds'],
                                bound_seconds=diag['attempts'][k]['bound_seconds'],
                                new_support_cells=int(new.sum())))
                        rows.append(dict(parent=parent,t_index=ti,t=float(z['t']),node=ni,mode=mode,
                            delta_real=float(delta.real),delta_imag=float(delta.imag),
                            new_support_cells=int(new.sum()),chosen_rank=diag['chosen_rank'],
                            fallback=diag['fallback'],returned_white_error=float(np.sqrt(2)/core.SIGMA*np.linalg.norm(result['F']-ref)),
                            returned_potential_error=abs(core.phi(result['F'],y)-ref_phi),
                            wall_seconds=diag['wall_seconds'],
                            conservative_RHS=after['RHS_conservative_total']-before['RHS_conservative_total'],
                            exact_seconds=diag.get('exact_seconds',0.),
                            projection_seconds=diag.get('current_shape_projection_seconds',0.),
                            last_white_bound=recorded[-1][2]['white_bound'] if recorded else None,
                            last_potential_bound=recorded[-1][2]['potential_bound'] if recorded else None))
            core.write(out/'PARTIAL.json',dict(completed_parents=pi+1,node_arm_rows=len(rows)))
        totals={}
        for mode,s in solvers.items():
            rr=[r for r in rows if r['mode']==mode];aa=[r for r in attempt_rows if r['mode']==mode]
            totals[mode]=dict(nodes=len(rr),qualified=sum(not r['fallback'] for r in rr),
                fallback=sum(r['fallback'] for r in rr),
                chosen_ranks={str(k):sum(r['chosen_rank']==k for r in rr) for k in (32,64,96,128,256)},
                new_conservative_RHS=sum(r['conservative_RHS'] for r in rr),
                node_wall_seconds=sum(r['wall_seconds'] for r in rr),
                full_cold_task_seconds=s.setup['basis_seconds']+s.bound_object.setup_seconds+sum(r['wall_seconds'] for r in rr),
                basis_setup_seconds=s.setup['basis_seconds'],bound=s.bound_object.summary(),
                model_setup_seconds=s.model_setup_seconds,
                max_returned_white_error=max(r['returned_white_error'] for r in rr),
                max_returned_potential_error=max(r['returned_potential_error'] for r in rr),
                all_sampled_white_bounds_hold=all(r['bound_dominates_sampled_error'] for r in aa),
                all_sampled_potential_bounds_hold=all(r['potential_bound_dominates_sampled_error'] for r in aa))
        status='COMPLETE'
        summary=dict(status=status,scope='EXPOSED_R10_DEV_PATHS_8_PARENTS_24_SHAPES_192_NODES',
            bases='fixed Green block Krylov from C and E; no anchor-shape dependence',
            ranks=[32,64,96,128],modes=totals,shapes=shape_rows,
            actual_y_used=True,GT_access=False,RHS576=0,new_physics_RHS_cap=20000,
            wall_seconds=time.perf_counter()-started,
            qualification='NUMERICAL_FLOAT64_DISCRETE_CONDITIONAL_BOUND_NOT_PDE_OR_INTERVAL_CERTIFICATE',
            cached_exact_reference='historical R10 followup DIRECT_LU; no new reference solves',
            timing='one matched rotated two-arm controlled window; no statistical speedup claim')
        core.write(out/'SUMMARY.json',summary)
    except BaseException as exc:
        error=repr(exc);core.write(out/'ERROR.json',dict(error=error,traceback=traceback.format_exc()));raise
    finally:
        csv_write(out/'NODE_RESULTS.csv',rows);csv_write(out/'ATTEMPT_RESULTS.csv',attempt_rows)
        csv_write(out/'SHAPE_SUPPORT.csv',shape_rows)
        core.write(out/'RECEIPT.json',ledger.close(status,error=error))


if __name__=='__main__':main()
