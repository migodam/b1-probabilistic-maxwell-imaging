"""Exactly thirty frozen DCTMC-style iterations; public observations only."""
from pathlib import Path
from collections import Counter
from datetime import datetime, timezone
import json
import resource
import sys
import time
import traceback
import numpy as np
from scipy import linalg as la
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parents[1]
sys.path.insert(0, str(REPO))
from research.research_fbg_r3_experiments.src.physics import MaxwellModel
from research.research_max_struct_r2.experiments.run_blind import dump, metric
from research.research_max_struct_r4.experiments.cooperative_paths import refine
from research.research_max_struct_r4.src.tmatrix_completion_control import TMatrixCompletionControl

# No shape/amplitude/host metadata is loaded; grid_shape defines discretization.
PUBLIC_FIELDS = ('G', 'E0', 'S', 'Y', 'grid_shape')


def utc():
    return datetime.now(timezone.utc).isoformat()


def run():
    started, cpu = time.perf_counter(), time.process_time()
    out = ROOT/'results/tmatrix_control'
    out.mkdir(parents=True, exist_ok=True)
    receipt_path = out/'RECEIPT.json'
    if receipt_path.exists():
        raise FileExistsError('Preserve existing T-matrix control disposition; no rerun or tuning')
    cfg = json.loads((ROOT/'configs/tmatrix_control.json').read_text())
    if cfg['iterations'] != 30 or cfg['corrector_calls'] != 1:
        raise ValueError('Frozen first control requires30 iterations and exactly one corrector')
    model, rows, materials = None, [], []
    receipt = dict(status='RUNNING', run_id='tmatrix_control', stage='DEV_KNOWN_METHOD_CONTROL',
                   started_utc=utc(), truth_access='NONE', public_fields=list(PUBLIC_FIELDS),
                   python_executable=sys.executable, startup_failure_receipt='STARTUP_FAILURE.json',
                   configuration=cfg, completed_iterations=0, corrector_started=False,
                   label='DCTMC-style variation, not a new invention or faithful paper reproduction',
                   literature_audit='Primary §III/IV/V.1 inspected: diagonal extraction is allowed; intermediate T/V generally nonsymmetric',
                   failures=[], limits='retained data compatibility and decreasing locality defect are not original material recovery or convergence certificates')
    dump(receipt_path, receipt)
    dump(out/'FROZEN_CONFIG.json', cfg)
    (out/'candidates').mkdir(exist_ok=True)
    try:
        data_path = REPO/'research/research_max_struct_r2/data_public'/f"{cfg['case_id']}.npz"
        with np.load(data_path, allow_pickle=False) as data:
            public = {key: data[key] for key in PUBLIC_FIELDS}
        model = MaxwellModel(public)
        Y = public['Y']
        if model.n != 576:
            raise ValueError('Frozen control expects576 independent material cells')
        control = TMatrixCompletionControl(model,Y,cfg['product_singular_relative_cutoff'],cfg['material_box'])
        receipt['data_projection'] = control.data_info
        receipt['historical_regularizer'] = dict(smooth=model.reg.smooth,tv=model.reg.tv,
                                                epsilon=model.reg.epsilon,edges=model.reg.edges)
        m = np.full(model.n,cfg['initial_material'])
        Tphysical,state,host = control.physical(m)
        T,completed = control.overwrite(Tphysical)
        host.update(iteration=0, **completed, candidate='candidates/iteration_000.npz')
        rows.append(host);materials.append(m.copy())
        np.savez_compressed(out/host['candidate'],m=m,Y_prediction=state['Y'])
        dump(out/'ITERATIONS.json',rows)
        receipt.update(last_completed_candidate=host,physical_and_dense_cost=dict(model.cost))
        dump(receipt_path,receipt)
        best_index = 0
        for iteration in range(1,cfg['iterations']+1):
            iteration_before = Counter(model.cost)
            m,V,locality = control.localize(T)
            Tphysical,state,physical = control.physical(m)
            T,completed = control.overwrite(Tphysical)
            row = dict(iteration=iteration,**locality,**physical,**completed,
                       candidate=f'candidates/iteration_{iteration:03d}.npz',
                       iteration_cost=dict(Counter(model.cost)-iteration_before))
            rows.append(row);materials.append(m.copy())
            np.savez_compressed(out/row['candidate'],m=m,Y_prediction=state['Y'])
            if row['physical_Phi'] < rows[best_index]['physical_Phi']:
                best_index = iteration
            dump(out/'ITERATIONS.json',rows)
            receipt.update(completed_iterations=iteration,last_completed_candidate=row,
                           selected_candidate_iteration=best_index,physical_and_dense_cost=dict(model.cost))
            dump(receipt_path,receipt)
            if iteration==1 or iteration%5==0:
                print(json.dumps(dict(iteration=iteration,Phi=row['physical_Phi'],
                      physical_residual=row['physical_relative_residual'],
                      completed_residual=row['completed_full_data_relative_residual'],
                      locality_offdiag_relative=row['offdiagonal_relative_to_full'],best_iteration=best_index)),flush=True)
        np.savez_compressed(out/'ALL_PHYSICAL_CANDIDATES.npz',m=np.asarray(materials),
                            iterations=np.arange(len(materials)),
                            physical_Phi=np.array([row['physical_Phi']for row in rows]),
                            physical_relative_residual=np.array([row['physical_relative_residual']for row in rows]))
        np.savez_compressed(out/'FINAL_LIFTING_STATE.npz',T_completed=T,T_physical=Tphysical,V_before_projection=V,m_last=m)
        receipt['completion_cost']=dict(model.cost)
        receipt['completion_wall_seconds']=time.perf_counter()-started
        seed=materials[best_index]
        np.savez_compressed(out/'SELECTED_BEFORE_CORRECTOR.npz',m=seed)
        receipt.update(status='CORRECTOR_RUNNING',corrector_started=True,
                       selected_before_corrector=rows[best_index],selected_candidate_iteration=best_index)
        dump(receipt_path,receipt)
        # The single selected physical candidate receives the exact common R4 control.
        final,corr=refine(model,Y,seed,maximum=cfg['corrector_maxiter'])
        dump(out/'CORRECTOR.json',corr)
        final_metric=metric(model,Y,final)
        value,gradient=model.data_value_gradient(final,Y)
        projected=final-np.clip(final-gradient,*cfg['material_box'])
        state=model.state(final)
        np.savez_compressed(out/'prediction.npz',m=final,Y_prediction=state['Y'],
                            gradient=gradient,projected_gradient=projected)
        receipt.update(status='SEALED_DEVELOPMENT_PREDICTION',finished_utc=utc(),
                       prediction='prediction.npz',physical_objective=float(value),
                       physical_relative_residual=final_metric['relative_residual'],
                       projected_gradient_inf=float(la.norm(projected,np.inf)),
                       corrector_stop=corr['stop'],corrector_cost=corr['cost'],
                       corrector_iterations=corr['iterations'],corrector_evaluations=corr['evaluations'],
                       physical_candidate_count=len(materials),corrector_calls=1,
                       final_material_distance_from_selected=float(la.norm(final-seed)))
    except Exception as exc:
        receipt.update(status='FAILED_PAID_CONTROL',error=repr(exc),traceback=traceback.format_exc())
        receipt['failures'].append(dict(error=repr(exc),completed_iterations=receipt['completed_iterations'],utc=utc()))
        raise
    finally:
        peak=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        receipt.update(physical_and_dense_cost=dict(model.cost)if model else {},
                       operation_wall_seconds=dict(model.timings)if model else {},
                       operation_cpu_seconds=dict(model._physics.cpu_timings)if model else {},
                       wall_seconds=time.perf_counter()-started,cpu_seconds=time.process_time()-cpu,
                       process_peak_memory_bytes=int(peak if sys.platform=='darwin'else peak*1024),
                       cost_semantics='LU includes original physical LU and explicitly tagged nonlocal locality_LU; dense T uses retained physical LU; all n=576 RHS counted',
                       timer_scope='public load,SVD,30cycles,31candidate persistence,onecorrector,finalgradient;imports and finalreceiptwrite excluded',
                       FFT=0,GPU=0)
        dump(receipt_path,receipt)
    print(json.dumps({k:receipt[k]for k in ('status','completed_iterations','selected_candidate_iteration',
                     'physical_objective','physical_relative_residual','projected_gradient_inf','wall_seconds',
                     'physical_and_dense_cost')}),flush=True)


if __name__=='__main__':
    with threadpool_limits(limits=1):
        run()
