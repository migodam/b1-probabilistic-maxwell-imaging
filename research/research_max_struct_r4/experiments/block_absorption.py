"""Single frozen public curvedbar development run, no truth reads or tuning."""
from pathlib import Path
from collections import Counter
from datetime import datetime, timezone
import json
import sys
import time
import traceback
import argparse

PROCESS_WALL = time.perf_counter()
PROCESS_CPU = time.process_time()
ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parents[1]
parser=argparse.ArgumentParser()
parser.add_argument('--variant',choices=['block_absorption','self_absorption','block_absorption_cached'],default='block_absorption')
VARIANT=parser.parse_args().variant
OUT = ROOT / 'results' / VARIANT
OUT.mkdir(parents=True, exist_ok=True)
if (OUT / 'RECEIPT.json').exists():
    raise FileExistsError('Preserve recorded block-absorption attempt')
STARTUP = dict(status='RUNNING', started_utc=datetime.now(timezone.utc).isoformat(),
               executable=sys.executable, argv=sys.argv, truth_access='NONE')
(OUT / 'STARTUP.json').write_text(json.dumps(STARTUP, indent=2))
try:
    import numpy as np
    from threadpoolctl import threadpool_limits, threadpool_info
    sys.path.insert(0, str(REPO))
    from research.research_fbg_r3_experiments.src.physics import MaxwellModel
    from research.research_max_struct_r2.experiments.run_blind import (
        ALLOW, compact, dump, metric, save_stage)
    from research.research_max_struct_r4.experiments.cooperative_paths import refine
    from research.research_max_struct_r4.src.cie_block_absorption import (
        run_block_cie, BlockAbsorption, block_material_value_gradient)
    STARTUP.update(status='IMPORTED', wall_seconds=time.perf_counter()-PROCESS_WALL,
                   cpu_seconds=time.process_time()-PROCESS_CPU)
except Exception as exc:
    STARTUP.update(status='FAILED', error=repr(exc), traceback=traceback.format_exc(),
                   wall_seconds=time.perf_counter()-PROCESS_WALL,
                   cpu_seconds=time.process_time()-PROCESS_CPU)
    (OUT / 'STARTUP.json').write_text(json.dumps(STARTUP, indent=2))
    raise
(OUT / 'STARTUP.json').write_text(json.dumps(STARTUP, indent=2))


def wiring_check(model, cfg):
    wall, cpu = time.perf_counter(), time.process_time()
    counts = Counter()
    record = dict(status='RUNNING', truth_access='NONE', directions=1,
                  physical_model_cost={}, finite_difference_points=2,
                  scope='new block-material derivative wiring only; no identity suite')
    dump(OUT / 'GRADIENT_WIRING.json', record)
    try:
        blocks = BlockAbsorption(model.G, model.shape, counts)
        m = .2 + .01*np.sin(np.arange(model.n)+1)
        J = .1*model.E
        Q = model.E + model.G@J
        counts['G_columns'] += model.ns
        weights = 1/np.maximum(np.sum(abs(model.E)**2, axis=0), 1e-300)
        Os = np.diag(weights)
        direction = np.cos(np.arange(model.n)+.5)
        direction /= np.linalg.norm(direction)
        prior = lambda x: model.reg.evaluate(x)[:2]
        fun = lambda x: block_material_value_gradient(
            blocks, J, Q, x, .2, Os, 1., prior, counts)
        value, grad = fun(m)
        h = cfg['gradient_check_h']
        finite = (fun(m+h*direction)[0]-fun(m-h*direction)[0])/(2*h)
        analytic = float(grad@direction)
        error = abs(analytic-finite)/max(abs(analytic), abs(finite), 1e-12)
        passed = bool(np.isfinite(error) and error <= cfg['gradient_check_relative_tolerance'])
        record.update(status='PASS' if passed else 'FAIL', value=value, beta=.2, h=h,
                      analytic_directional=analytic, finite_difference_directional=finite,
                      relative_error=error, tolerance=cfg['gradient_check_relative_tolerance'])
        if not passed:
            raise AssertionError('New block-material gradient wiring failed')
    except Exception as exc:
        record.update(error=repr(exc), traceback=traceback.format_exc())
        if record['status']=='RUNNING':
            record['status']='FAILED'
        raise
    finally:
        record.update(wall_seconds=time.perf_counter()-wall,
                      cpu_seconds=time.process_time()-cpu, operation_counts=dict(counts))
        dump(OUT / 'GRADIENT_WIRING.json', record)
    return record


def run():
    cfg = json.loads((ROOT / 'configs' / (VARIANT+'.json')).read_text())
    dump(OUT / 'FROZEN_CONFIG.json', cfg)
    model = None
    receipt = dict(status='RUNNING', run_id=VARIANT, configuration=cfg,
        classification=cfg['classification'], truth_access='NONE', source='experiments/block_absorption.py',
        original_source='research/research_max_struct_r2/src/cie_global.py',
        failed_attempts=[], startup=STARTUP, original_prior_hessian_constructed=False,
        prior_evaluation='same model.reg.evaluate values and gradients; unused dense Hessian omitted',
        novelty='no originality assertion; block absorption and mature CIE are known',
        block_cost_semantics=('self mode uses only elementwise complex division; each division charged, no small LU'
            if cfg.get('gin_mode','block4')=='self' else
            ('cached inverse mode charges actual4x4 factor+4RHS solve once per cache miss; subsequent forward/adjoint matrix applies separately'
             if cfg.get('block_solver','direct')=='cached_inverse' else
             'each NumPy batched solve refactors all144 4x4 systems; every factor, solve, RHS counted')),
        threads=threadpool_info(), process_import_and_serialization_included=True)
    dump(OUT / 'RECEIPT.json', receipt)
    try:
        with np.load(REPO / f"research/research_max_struct_r2/data_public/{cfg['case_id']}.npz",
                     allow_pickle=False) as f:
            if not set(f.files) <= ALLOW:
                raise ValueError('Unexpected public key; refuse payload')
            public = {k: f[k] for k in f.files}
        model = MaxwellModel(public)
        if tuple(model.shape)!=(24,24) or float(public['host'])!=cfg['host']:
            raise ValueError('Frozen grid/host mismatch')
        Y = public['Y']
        m0 = np.full(model.n, cfg['host'])
        if cfg.get('gin_mode','block4')=='self' or cfg.get('block_solver','direct')=='cached_inverse':
            receipt['gradient_check']=dict(status='NOT_RERUN',
                reason='same exact generic material derivative; authorized implementation/control mode, no additional derivative suite',
                existing_check='results/block_absorption/GRADIENT_WIRING.json')
        else:
            receipt['gradient_check'] = wiring_check(model, cfg)
        initial = metric(model, Y, m0)
        prior = lambda x: model.reg.evaluate(x)[:2]
        numerical = dict(stages=cfg['stages'], method='cie', sourcewise_cg=False,
            gin_mode=cfg.get('gin_mode','block4'),
            block_solver=cfg.get('block_solver','direct'),
            cache_material_value_gradient=bool(cfg.get('cache_material_value_gradient',False)),
            dpcs_rtol=cfg['dpcs_rtol'], coefficient_tol=cfg['coefficient_tol'], miniter=cfg['miniter'],
            material_inner_maxiter=cfg['material_inner_maxiter'], material_box=cfg['material_box'],
            material_inner_ftol=cfg['material_inner_ftol'], material_inner_gtol=cfg['material_inner_gtol'],
            regularizer=prior, physical_sigma=float(np.linalg.norm(Y)),
            receipt_path=str(OUT/'AUXILIARY_EXECUTION.json'))
        def stage_callback(st):
            save_stage(OUT/f"stage_{st['stage']}.npz", st)
            print(json.dumps(dict(event='STAGE_SAVED', stage=st['stage'])), flush=True)
        numerical['stage_callback'] = stage_callback
        aux = run_block_cie(model, Y, m0, numerical)
        dump(OUT/'AUXILIARY.json', compact(aux))
        scored = [metric(model, Y, c['m']) for c in aux['candidates']]
        seed = min(scored, key=lambda c: c['objective'])
        np.savez_compressed(OUT/'SELECTED_BEFORE_CORRECTOR.npz', m=seed['m'])
        np.savez_compressed(OUT/'ALL_PHYSICAL_CANDIDATES.npz',
            m=np.array([c['m'] for c in scored]),
            original_Phi=np.array([c['objective'] for c in scored]))
        final, corr = refine(model, Y, seed['m'], maximum=cfg['corrector_maxiter'])
        dump(OUT/'CORRECTOR.json', corr)
        final_metric = metric(model, Y, final)
        best = min([initial, seed, final_metric], key=lambda c: c['objective'])
        state = model.state(best['m'])
        np.savez_compressed(OUT/'prediction.npz', m=best['m'], Y_prediction=state['Y'],
                            j=state['j'], Q=state['Q'])
        control_path = REPO/f"research/research_max_struct_r2/results/blind/{cfg['case_id']}/global/RECEIPT.json"
        control = json.loads(control_path.read_text())
        receipt.update(status='SEALED_DEVELOPMENT_PREDICTION',
            selection='original Phi with original prior; initial+auxiliary best+corrector',
            auxiliary_status=aux['status'], auxiliary_selected_original_Phi=seed['objective'],
            auxiliary_candidate_scores=[{k:v for k,v in c.items() if k!='m'} for c in scored],
            physical_objective=best['objective'], physical_relative_residual=best['relative_residual'],
            physical_whitened_residual=float(np.linalg.norm(state['Y']-Y)/float(public['sigma'])),
            linear_relative_residual=float(state['state_relative_residual']),
            material_min=float(best['m'].min()), material_max=float(best['m'].max()),
            corrector_stop=corr['stop'], corrector_iterations=corr['iterations'],
            auxiliary_operation_counts=aux['execution_receipt']['operation_counts'],
            prediction=f'results/{VARIANT}/prediction.npz',
            saved_R2_global=dict(path=str(control_path.relative_to(REPO)),
                auxiliary_original_Phi=control['auxiliary_selected_physical_objective'],
                final_original_Phi=control['physical_objective'],
                final_relative_residual=control['physical_relative_residual'],
                physical_model_cost=control['physical_model_cost'],
                wall_seconds=control['wall_seconds'], cpu_seconds=control['cpu_seconds'],
                comparability='same auxiliary schedule/initialization/prior; historical unscaled corrector differs from this R4 scaled corrector; no equal-compute assertion'),
            final_parameter_selection='NONE; one configured Gin mode, one case, one run')
    except Exception as exc:
        receipt.update(status='FAILED', error=repr(exc), traceback=traceback.format_exc())
        receipt['failed_attempts'].append(dict(error=repr(exc), included_in_total=True))
        if hasattr(exc, 'cie_execution_receipt'):
            receipt['failed_auxiliary_execution']=exc.cie_execution_receipt
        raise
    finally:
        receipt.update(wall_seconds=time.perf_counter()-PROCESS_WALL,
            cpu_seconds=time.process_time()-PROCESS_CPU,
            physical_model_cost=dict(model.cost) if model else {},
            physical_operation_wall=dict(model.timings) if model else {},
            physical_operation_cpu=dict(model._physics.cpu_timings) if model else {},
            finished_utc=datetime.now(timezone.utc).isoformat(),
            cost_scope='startup/imports/public load/wiring check if explicitly enabled/full1500 updates/physical ranking/corrector/artifact writes; final receipt write excluded',
            historical_operator_assembly='not rerun; public frozen operator load only')
        dump(OUT/'RECEIPT.json', receipt)
    print(json.dumps({k:receipt[k] for k in ['status','wall_seconds','cpu_seconds',
        'physical_objective','physical_relative_residual','physical_model_cost']}), flush=True)


if __name__=='__main__':
    with threadpool_limits(limits=1):
        run()
