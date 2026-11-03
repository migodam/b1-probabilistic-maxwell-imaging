"""One fixed pi-phase control of the auxiliary offdiagonal block Green, CPU1, no GT."""
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
parser.add_argument('--variant',choices=['onehop_phase_pi'],default='onehop_phase_pi')
VARIANT=parser.parse_args().variant
OUT = ROOT / 'results' / VARIANT
OUT.mkdir(parents=True, exist_ok=True)
if (OUT / 'RECEIPT.json').exists():
    raise FileExistsError('Preserve recorded pi-phase-control attempt')
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
    from research.research_max_struct_r4.src import cie_onehop as onehop

    class PhasePiAbsorption(onehop.OneHopAbsorption):
        """Change only the copied auxiliary Goff; model.G is never mutated.

        The inherited apply and exact material gradient both call off_action
        on this same instance. Consequently the sign change reaches all
        current/material forward and adjoint actions.
        """
        def __init__(self,G,shape,cost):
            super().__init__(G,shape,cost)
            self.goff *= -1
            self.cost['phase_pi_setup_complex_sign_changes'] += self.goff.size

    # Local to this runner process; the imported source file is unchanged.
    onehop.OneHopAbsorption=PhasePiAbsorption
    OneHopAbsorption=PhasePiAbsorption
    run_onehop_cie=onehop.run_onehop_cie
    STARTUP.update(status='IMPORTED', wall_seconds=time.perf_counter()-PROCESS_WALL,
                   cpu_seconds=time.process_time()-PROCESS_CPU)
except Exception as exc:
    STARTUP.update(status='FAILED', error=repr(exc), traceback=traceback.format_exc(),
                   wall_seconds=time.perf_counter()-PROCESS_WALL,
                   cpu_seconds=time.process_time()-PROCESS_CPU)
    (OUT / 'STARTUP.json').write_text(json.dumps(STARTUP, indent=2))
    raise
(OUT / 'STARTUP.json').write_text(json.dumps(STARTUP, indent=2))


def run():
    cfg = json.loads((ROOT / 'configs' / (VARIANT+'.json')).read_text())
    dump(OUT / 'FROZEN_CONFIG.json', cfg)
    import shutil
    snapshot=OUT/'executed_source_snapshot';snapshot.mkdir(exist_ok=True)
    for source in (Path(__file__),ROOT/'src/cie_onehop.py',ROOT/'configs/onehop_phase_pi.json'):
        shutil.copy2(source,snapshot/source.name)
    model = None
    receipt = dict(status='RUNNING', run_id=VARIANT, configuration=cfg,
        classification=cfg['classification'], truth_access='NONE', source='experiments/onehop_phase_pi.py',
        original_source='research/research_max_struct_r2/src/cie_global.py',
        failed_attempts=[], startup=STARTUP, original_prior_hessian_constructed=False,
        prior_evaluation='same model.reg.evaluate values and gradients; unused dense Hessian omitted',
        novelty='no originality assertion; single auxiliary phase-sign mechanism control',
        phase_control=dict(auxiliary_Goff_multiplier=-1, physical_G_unchanged=True,
            residual='R_pi=(I-L)u; L=diag(a)Goff_original',
            interpretation='same abs(Goff), graph, diagonal metric, whole-box sufficient bound and arithmetic per evaluation; only the real cross term changes sign',
            cost_scope='same operation forms per evaluation; optimizer trajectory and total evaluation counts may differ'),
        block_cost_semantics='actual scalar divisions and block forward/adjoint RHS; zero small LU/inverse',
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
        bounds_counter=Counter()
        bounds_blocks=OneHopAbsorption(model.G,model.shape,bounds_counter)
        receipt['whole_box_bound']=[bounds_blocks.whole_box_bound(cfg['material_box'],beta) for beta in [.2,1.]]
        receipt['bound_setup_operation_counts']=dict(bounds_counter)
        dump(OUT/'WHOLE_BOX_BOUND.json',receipt['whole_box_bound'])
        print(json.dumps(dict(event='WHOLE_BOX_BOUND',bounds=receipt['whole_box_bound'])),flush=True)
        receipt['gradient_check'] = dict(status='NOT_RERUN',
            reason='only fixed auxiliary Goff signs changed; same generic exact material derivative',
            original_check='results/onehop_absorption/GRADIENT_WIRING.json',
            extra_directions=0,physical_model_cost={})
        dump(OUT/'GRADIENT_WIRING.json',receipt['gradient_check'])
        initial = metric(model, Y, m0)
        prior = lambda x: model.reg.evaluate(x)[:2]
        numerical = dict(stages=cfg['stages'], method='cie', sourcewise_cg=False,
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
        def phase_receipt(execution):
            execution.update(phase_control=receipt['phase_control'],
                residual='u-diag(a)Goff_original u; implemented by one fixed negative auxiliary Goff copy')
            dump(OUT/'AUXILIARY_EXECUTION.json',execution)
        numerical['cost_callback']=phase_receipt
        aux = run_onehop_cie(model, Y, m0, numerical)
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
            final_parameter_selection='NONE; one fixed auxiliary pi-phase sign control, one disclosed development case, one run')
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
            cost_scope='startup/imports/public load/analytic bound/full1500 updates/physical ranking/corrector/artifact writes; no derivative test rerun; final receipt write excluded',
            historical_operator_assembly='not rerun; public frozen operator load only')
        dump(OUT/'RECEIPT.json', receipt)
    print(json.dumps({k:receipt[k] for k in ['status','wall_seconds','cpu_seconds',
        'physical_objective','physical_relative_residual','physical_model_cost']}), flush=True)


if __name__=='__main__':
    with threadpool_limits(limits=1):
        run()
