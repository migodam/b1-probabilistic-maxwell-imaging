"""Optional rank128 and inversion-interface qualification, within remaining kernel reserve."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'):os.environ[k]='1'
import sys,json,time,csv,traceback
from pathlib import Path
import numpy as np
D=Path(__file__).resolve().parent;P=D.parent;sys.path.insert(0,str(P/'shared'))
from core import load_setup,R10,SIGMA,write,ExactSolver
from fixed_total_field import TotalFieldSolver,CostLedger

def main():
    previous=json.loads((D/'KERNEL_RECEIPT.json').read_text());left=30000-previous['RHS_conservative_total']
    book=CostLedger(D/'SUPPLEMENT_COST_LEDGER.jsonl',left);rows=[];ss=[];status='FAILED';error=None
    try:
        S=load_setup();C,O,E,Fh,roi=S['C'],S['O'],S['E'],S['Fh'],S['ROI']
        for pi in range(8):
            parent=f'R10BI_FRESH_{pi:03d}';root=R10/'followup_shape_reuse';a0=np.load(root/f'{parent}_ANCHOR.npz')['fraction'][roi]
            solver=TotalFieldSolver(C,O,E,Fh,a_anchor=a0,ranks=(128,),basis_kind='anchor',ledger=book)
            y=np.load(R10/'results/candidates_fresh10/inputs/online'/f'{parent}.npz')['observation']
            for ti in range(3):
                saved=np.load(root/f'{parent}_t{ti}.npz');a=saved['fraction'][roi];ds=saved['delta'];ref=np.load(D/f'{parent}_t{ti}_RESULT.npz')['F_exact']
                ss.append(dict(parent=parent,t_index=ti,basis_seconds=solver.setup['basis_seconds']))
                approx=[];returned=[]
                for ni,delta in enumerate(ds):
                    out=solver.evaluate(a,delta,y=y);di=out['diagnostics'];at=di['attempts'][0]
                    # Recover candidate exactly from already projected model; no duplicate solve or free solve.
                    # During evaluate callback, capture all computed arrays rather than recompute them.
                    fa=solver.last_approx_F
                    ae=float(np.sqrt(2)/SIGMA*np.linalg.norm(fa-ref[ni]))
                    rows.append(dict(parent=parent,t_index=ti,node=ni,qualified=not di['fallback'],fallback=di['fallback'],
                                     white_bound=at['white_bound'],potential_bound=at['potential_bound'],actual_approx_white_error=ae,
                                     returned_white_error=float(np.sqrt(2)/SIGMA*np.linalg.norm(out['F']-ref[ni])),
                                     bound_holds=ae<=at['white_bound']*(1+1e-8)+1e-9,
                                     solve_seconds=at['solve_seconds'],bound_seconds=at['bound_seconds'],
                                     shape_projection_seconds=di['current_shape_projection_seconds'],wall_seconds=di['wall_seconds']))
                    approx.append(fa);returned.append(out['F'])
                np.savez_compressed(D/f'{parent}_t{ti}_R128.npz',F_approx=np.stack(approx),F_returned=np.stack(returned))
        # Inversion API: independent shared current-state implementation is the comparison.
        # Both algorithms' actual work is reserved before evaluation; no finite-difference tuning.
        a=np.load(R10/'followup_shape_reuse/R10BI_FRESH_000_t1.npz')['fraction'][roi];delta=1.31+.11j
        s=TotalFieldSolver(C,O,E,Fh,a_anchor=a0,ranks=(32,),basis_kind='anchor',ledger=book)
        out=s.evaluate(a,delta,need_jac=True)
        book.reserve(256,40,'INDEPENDENT_SHARED_EXACT_JACOBIAN_REFERENCE')
        ref=ExactSolver(C,O,E,Fh).evaluate(a,delta,need_jac=True)
        api=dict(status='PASS',observations='source-major8sourcesx32receivers',J_occ_shape=list(out['J_occ'].shape),J_delta_shape=list(out['J_delta'].shape),
                 F_white_difference=float(np.sqrt(2)/SIGMA*np.linalg.norm(out['F']-ref['F'])),
                 J_occ_relative_difference=float(np.linalg.norm(out['J_occ']-ref['J_occ'])/np.linalg.norm(ref['J_occ'])),
                 J_delta_relative_difference=float(np.linalg.norm(out['J_delta']-ref['J_delta'])/np.linalg.norm(ref['J_delta'])),
                 J_exact=True,approximate_J='NOT_CONTROLLED_NOT_USED',RHS_per_J_call=40,LU_per_J_call=1)
        if api['F_white_difference']>1e-8 or max(api['J_occ_relative_difference'],api['J_delta_relative_difference'])>1e-10:raise RuntimeError('API_COMPARISON_FAIL')
        write(D/'API_QUALIFICATION.json',api)
        basis=sum(s['basis_seconds'] for s in ss[::3]);online=sum(r['wall_seconds'] for r in rows)
        write(D/'R128_SUMMARY.json',dict(status='COMPLETE',nodes=192,qualified=sum(r['qualified'] for r in rows),fallback_nodes=sum(r['fallback'] for r in rows),
               max_actual_approx_white_error=max(r['actual_approx_white_error'] for r in rows),median_actual_approx_white_error=float(np.median([r['actual_approx_white_error'] for r in rows])),
               actual_error_below_point01=sum(r['actual_approx_white_error']<=.01 for r in rows),sampled_bounds_all_hold=all(r['bound_holds'] for r in rows),
               max_returned_white_error=max(r['returned_white_error'] for r in rows),fixed_basis_setup_seconds=basis,online_seconds=online,
               cold_task_including_fallback_seconds=basis+online,scope='Optional rank128 same24DEVshapes192nodes; original30k kernel reserve unchanged'))
        status='COMPLETE'
    except BaseException as e:error=repr(e);write(D/'SUPPLEMENT_ERROR.json',dict(error=error,traceback=traceback.format_exc()));raise
    finally:
        if rows:
            with (D/'R128_NODE_RESULTS.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
        rec=book.close(status,error=error);rec['previous_kernel_RHS']=previous['RHS_conservative_total'];rec['combined_RHS_conservative_total']=previous['RHS_conservative_total']+rec['RHS_conservative_total'];write(D/'SUPPLEMENT_RECEIPT.json',rec)
if __name__=='__main__':main()
