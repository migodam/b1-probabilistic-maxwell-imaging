"""One frozen extended mature native CIE control; public observations only."""
from pathlib import Path
from collections import Counter
from contextlib import contextmanager
from datetime import datetime, timezone
import json
import sys
import time
import traceback

START_WALL, START_CPU = time.perf_counter(), time.process_time()
ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parents[1]
OUT = ROOT/'results/native_extended'
OUT.mkdir(parents=True, exist_ok=True)
if (OUT/'RECEIPT.json').exists():
    raise FileExistsError('Preserve extended native attempt')
startup = dict(status='RUNNING', executable=sys.executable, argv=sys.argv,
               started_utc=datetime.now(timezone.utc).isoformat(), truth_access='NONE')
(OUT/'STARTUP.json').write_text(json.dumps(startup, indent=2))
try:
    import numpy as np
    from threadpoolctl import threadpool_limits, threadpool_info
    sys.path.insert(0, str(REPO))
    from research.research_fbg_r3_experiments.src.physics import MaxwellModel
    from research.research_fbg_r3_experiments.src import cie_tsom_baseline as native
    from research.research_max_struct_r2.experiments.run_blind import ALLOW, dump, compact, metric
    from research.research_max_struct_r4.experiments.cooperative_paths import refine
    startup.update(status='IMPORTED', wall_seconds=time.perf_counter()-START_WALL,
                   cpu_seconds=time.process_time()-START_CPU)
except Exception as exc:
    startup.update(status='FAILED', error=repr(exc), traceback=traceback.format_exc(),
                   wall_seconds=time.perf_counter()-START_WALL,
                   cpu_seconds=time.process_time()-START_CPU)
    (OUT/'STARTUP.json').write_text(json.dumps(startup, indent=2))
    raise
(OUT/'STARTUP.json').write_text(json.dumps(startup, indent=2))


@contextmanager
def observe_readouts(model, cfg, progress):
    """Forwarding observation only: no changed numerical arguments or results."""
    original = native._readout
    boundaries = np.cumsum([s['maxiter'] for s in cfg['stages']])
    def readout(J, Q, *args, **kwargs):
        # Keep the actual mutable native Counter available even if a solve fails.
        progress['native_cost_reference'] = args[-1]
        m, meta = original(J, Q, *args, **kwargs)
        progress['readouts'] += 1
        count = progress['readouts']
        stage = int(np.searchsorted(boundaries, count, side='left'))
        progress['source_mode_iterations'] += model.ns*(2*cfg['stages'][stage]['mf'])**2
        if count % cfg['checkpoint_every'] == 0:
            path = OUT/'checkpoints'/f'iteration_{count:05d}.npz'
            path.parent.mkdir(parents=True, exist_ok=True)
            np.savez_compressed(path, m=m, J=J, Q=Q,
                                global_iteration=count, stage=stage)
            if count in boundaries:
                np.savez_compressed(OUT/f'stage_{stage}.npz', m=m, J=J, Q=Q)
            record = dict(global_iteration=count, stage=stage,
                source_mode_iterations=progress['source_mode_iterations'],
                wall_seconds=time.perf_counter()-START_WALL,
                cpu_seconds=time.process_time()-START_CPU,
                native_counts_before_iteration_footer=dict(args[-1]),
                physical_model_cost=dict(model.cost),
                material_min=float(m.min()), material_max=float(m.max()),
                path=str(path.relative_to(ROOT)))
            progress['checkpoints'].append(record)
            dump(OUT/'CHECKPOINTS.json', progress['checkpoints'])
            print(json.dumps(dict(event='NATIVE_CHECKPOINT', global_iteration=count,
                stage=stage, wall_seconds=record['wall_seconds'])), flush=True)
        return m, meta
    native._readout = readout
    try:
        yield
    finally:
        native._readout = original


def run():
    cfg = json.loads((ROOT/'configs/native_extended.json').read_text())
    dump(OUT/'FROZEN_CONFIG.json', cfg)
    model = None
    progress = dict(readouts=0, source_mode_iterations=0, checkpoints=[])
    rec = dict(status='RUNNING', configuration=cfg, truth_access='NONE', startup=startup,
        classification=cfg['classification'], failed_attempts=[],
        original_source='research/research_fbg_r3_experiments/src/cie_tsom_baseline.py',
        prior_evaluation='same reg.evaluate values/gradients; unused dense IRLS matrix omitted',
        numerical_source_modified=False, threads=threadpool_info(),
        derivative_check='NOT_NEEDED: unchanged known algorithm, no new derivative',
        cost_scope='startup/import/public load/full auxiliary/physical stage validations/common corrector/prediction/artifact writes; no operator assembly',
        equal_compute_assertion=False)
    dump(OUT/'RECEIPT.json', rec)
    try:
        with np.load(REPO/f"research/research_max_struct_r2/data_public/{cfg['case_id']}.npz",
                     allow_pickle=False) as f:
            if not set(f.files) <= ALLOW:
                raise ValueError('Unexpected public key')
            public = {k:f[k] for k in f.files}
        model = MaxwellModel(public)
        if tuple(model.shape) != (24,24) or float(public['host']) != cfg['host']:
            raise ValueError('Frozen public host/grid mismatch')
        Y = public['Y']
        m0 = np.full(model.n, cfg['host'])
        numerical = {k:cfg[k] for k in ('stages', 'dpcs_rtol', 'sourcewise_cg',
            'coefficient_tol', 'miniter', 'material_inner_maxiter', 'material_inner_ftol',
            'material_inner_gtol', 'material_box')}
        numerical.update(method='cie', regularizer=lambda m:model.reg.evaluate(m)[:2],
                         physical_sigma=float(np.linalg.norm(Y)))
        aux_cpu = time.process_time()
        with observe_readouts(model, cfg, progress):
            aux = native.run_cie(model, Y, m0, numerical)
        rec['auxiliary_cpu_seconds'] = time.process_time()-aux_cpu
        dump(OUT/'AUXILIARY.json', compact(aux))
        selected_index = min(range(len(aux['candidates'])),
                             key=lambda i:aux['candidates'][i]['objective'])
        seed = aux['candidates'][selected_index]
        np.savez_compressed(OUT/'ALL_PHYSICAL_CANDIDATES.npz',
            m=np.asarray([c['m'] for c in aux['candidates']]),
            original_Phi=np.asarray([c['objective'] for c in aux['candidates']]),
            stage=np.asarray([c['stage'] for c in aux['candidates']]))
        np.savez_compressed(OUT/'SELECTED_BEFORE_CORRECTOR.npz', m=seed['m'])
        rec.update(auxiliary_status=aux['status'], auxiliary_cost=aux['cost'],
            auxiliary_selected_stage=seed['stage'], auxiliary_selected_original_Phi=seed['objective'],
            auxiliary_candidate_scores=[{k:v for k,v in c.items() if k!='m'} for c in aux['candidates']])
        dump(OUT/'RECEIPT.json', rec)
        print(json.dumps(dict(event='AUXILIARY_DONE', status=aux['status'],
            selected_stage=seed['stage'], original_Phi=seed['objective'])), flush=True)
        final, corr = refine(model, Y, seed['m'], maximum=cfg['corrector_maxiter'])
        dump(OUT/'CORRECTOR.json', corr)
        final_metric = metric(model, Y, final)
        # refine retains its seed and only accepts original-Phi improvements.
        best_m = final if final_metric['objective'] <= seed['objective'] else seed['m']
        value, gradient = model.data_value_gradient(best_m, Y)
        state = model.state(best_m)
        pg = best_m-np.clip(best_m-gradient, *cfg['material_box'])
        np.savez_compressed(OUT/'prediction.npz', m=best_m, Y_prediction=state['Y'],
                            j=state['j'], Q=state['Q'])
        iterations = [x for x in aux['history'] if 'global_iteration' in x]
        stages = []
        for k in range(len(cfg['stages'])):
            rows = [x for x in iterations if x['stage']==k]
            stages.append(dict(stage=k, iterations=len(rows), mf=cfg['stages'][k]['mf'],
                fourier_columns=rows[-1]['fourier_columns'], apcs_rank=rows[-1]['apcs_rank'],
                source_mode_iterations=sum(model.ns*x['fourier_columns'] for x in rows),
                material_inner_iterations=sum(x.get('material_inner_nit',0) for x in rows)))
        rec.update(status='SEALED_DEVELOPMENT_PREDICTION',
            physical_objective=float(value),
            physical_relative_residual=float(np.linalg.norm(state['Y']-Y)/np.linalg.norm(Y)),
            physical_whitened_residual=float(np.linalg.norm(state['Y']-Y)/float(public['sigma'])),
            noise_radius99=float(public['noise_radius']),
            physical_linear_residual=float(state['state_relative_residual']),
            projected_gradient_inf=float(np.linalg.norm(pg,np.inf)),
            material_min=float(best_m.min()), material_max=float(best_m.max()),
            corrector_stop=corr['stop'], corrector_iterations=corr['iterations'],
            actual_auxiliary_iterations=len(iterations), stages=stages,
            source_mode_iterations=sum(x['source_mode_iterations'] for x in stages),
            source_mode_iterations_definition='sum over completed updates of ns times native Fourier coefficient columns; work proxy, not elapsed-time equivalence or FLOPs',
            dense_current_core_complex_multiply_add_estimate=(5*model.n+4*model.nr)*sum(x['source_mode_iterations'] for x in stages),
            dense_current_core_estimate_scope='nine per-update coefficient matrix products including two current reconstructions and repeated receiver residual in history; excludes setup, elementwise work, material inner work, SVD, validation and corrector',
            prediction='results/native_extended/prediction.npz')
    except Exception as exc:
        rec.update(status='FAILED', error=repr(exc), traceback=traceback.format_exc())
        rec['failed_attempts'].append(dict(error=repr(exc), included_in_total=True))
        raise
    finally:
        rec.update(wall_seconds=time.perf_counter()-START_WALL,
            cpu_seconds=time.process_time()-START_CPU,
            observed_readouts=progress['readouts'],
            observed_source_mode_iterations=progress['source_mode_iterations'],
            native_cost_last_observed=dict(progress.get('native_cost_reference',{})),
            physical_model_cost=dict(model.cost) if model else {},
            physical_operation_wall=dict(model.timings) if model else {},
            physical_operation_cpu=dict(model._physics.cpu_timings) if model else {},
            finished_utc=datetime.now(timezone.utc).isoformat())
        dump(OUT/'RECEIPT.json', rec)
    print(json.dumps({k:rec[k] for k in ('status','wall_seconds','cpu_seconds',
        'physical_objective','physical_relative_residual','physical_model_cost')}), flush=True)


if __name__=='__main__':
    with threadpool_limits(limits=1):
        run()
