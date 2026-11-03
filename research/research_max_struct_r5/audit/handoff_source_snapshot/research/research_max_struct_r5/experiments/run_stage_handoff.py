"""Frozen natural-stage state handoff, followed by the original R4 refine."""
from pathlib import Path
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key]='1'
import argparse
from collections import Counter
from datetime import datetime,timezone
import json
import sys
import time
import traceback

WALL,CPU=time.perf_counter(),time.process_time()
ROOT=Path(__file__).resolve().parents[1]
REPO=ROOT.parents[1]
sys.path.insert(0,str(REPO))
from research.research_max_struct_r5.src.trajectory_replay import TrajectoryRecorder,dump
parser=argparse.ArgumentParser()
parser.add_argument('--handoff-id',required=True)
args=parser.parse_args()
config_path=ROOT/'configs/handoff_frozen.json'
cfg=json.loads(config_path.read_text())
arm_cfg=cfg['handoffs'][args.handoff_id]
case_id=cfg['case_id']
public_files=[cfg['public_data']]
numerical_cfg=arm_cfg['target']
out=ROOT/'RESULTS/handoff'/args.handoff_id
out.mkdir(parents=True,exist_ok=False)
dump(out/'FROZEN_CONFIG.json',cfg)
receipt=dict(status='RUNNING',handoff_id=args.handoff_id,case_id=case_id,configuration=cfg,
    configuration_path=str(config_path),public_files=public_files,
    started_utc=datetime.now(timezone.utc).isoformat(),executable=sys.executable,
    truth_access='NONE',scoring='NOT_RUN',classification='DISCLOSED_DEVELOPMENT_NATURAL_STAGE_HANDOFF_NOT_BLIND',
    seed=None,seed_reason='Deterministic inverse solve of unchanged public observation; no random draws.',
    failed_attempts=[])
dump(out/'RECEIPT.json',receipt)
model=observer=None
try:
    import numpy as np
    from threadpoolctl import threadpool_limits,threadpool_info
    from research.research_fbg_r3_experiments.src.physics import MaxwellModel
    from research.research_max_struct_r5.src.cie_stage_handoff import run_stage_handoff_cie
    from research.research_max_struct_r4.experiments.cooperative_paths import refine
    with threadpool_limits(limits=1):
        receipt['threads']=threadpool_info()
        start,ct=time.perf_counter(),time.process_time()
        public={}
        for public_file in public_files:
            with np.load(REPO/public_file,allow_pickle=False) as data:
                forbidden={'truth','alpha','mask','yclean','seed','family'}
                if forbidden & {key.lower() for key in data.files}:
                    raise ValueError('Unexpected private payload key')
                public.update({key:data[key] for key in data.files})
        model=MaxwellModel(public);Y=public['Y']
        with np.load(REPO/arm_cfg['source_stage'],allow_pickle=False) as initial:
            m0=initial['m'].copy();initial_J=initial['J'].copy();initial_A=initial['coefficients'].copy()
        with np.load(REPO/arm_cfg['source_candidates'],allow_pickle=False) as initial:
            host_m=initial['m'][0].copy()
            phi_key='original_Phi' if 'original_Phi' in initial.files else 'objective'
            host_Phi=float(initial[phi_key][0])
            host_relative_residual=float(initial['physical_relative_residual'][0])
        if not np.array_equal(host_m,np.full(model.n,cfg['common']['host'])):
            raise ValueError('Paid historical incumbent is not common initial host')
        receipt['source_stage']=arm_cfg['source_stage']
        receipt['host_Phi_reused_without_new_forward']=dict(original_Phi=host_Phi,
            physical_relative_residual=host_relative_residual,source=arm_cfg['source_candidates'],row=0)
        for key,value in cfg.get('expected_public',{}).items():
            if not np.array_equal(np.asarray(public[key]),np.asarray(value)):
                raise ValueError('Frozen public host/prior/box mismatch: '+key)
        receipt['adaptation_differences']=cfg.get('adaptation_differences',[])
        observer=TrajectoryRecorder(out,model,m0)
        receipt['setup_wall_seconds']=time.perf_counter()-start
        receipt['setup_cpu_seconds']=time.process_time()-ct
        numerical={key:value for key,value in cfg['common'].items() if key!='host'}
        numerical.update(numerical_cfg,method='cie',regularizer=lambda m:model.reg.evaluate(m)[:2],
            physical_sigma=float(np.linalg.norm(Y)),
            stage_callback=observer.stage,receipt_path=str(out/'AUXILIARY_EXECUTION.json'))
        def handoff_callback(state):
            dump(out/'HANDOFF_WIRING.json',state['info'])
            np.savez_compressed(out/'HANDOFF_WIRING.npz',
                **{key:value for key,value in state.items() if key!='info'})
        numerical.update(initial_mf=arm_cfg['initial_mf'],
            initial_fourier_coefficients=initial_A,initial_current=initial_J,
            handoff_callback=handoff_callback)
        start,ct,before=time.perf_counter(),time.process_time(),Counter(model.cost)
        auxiliary=run_stage_handoff_cie(model,Y,m0,numerical)
        receipt['auxiliary_wall_seconds']=time.perf_counter()-start
        receipt['auxiliary_cpu_seconds']=time.process_time()-ct
        receipt['auxiliary_physical_model_cost']=dict(Counter(model.cost)-before)
        receipt['auxiliary_execution_receipt']=auxiliary['execution_receipt']
        receipt['auxiliary_status']=auxiliary['status']
        dump(out/'HISTORY.json',auxiliary['history'])
        dump(out/'AUXILIARY_COST.json',auxiliary['cost'])
        # Resume initial material is the original stage0 endpoint. Reuse only
        # the already-paid common host Phi in addition to newly validated
        # resumed initial material and the two true remaining-stage endpoints.
        host_candidate=dict(m=host_m,stage=-2,reason='reused_paid_original_host',
            objective=host_Phi,physical_relative_residual=host_relative_residual)
        all_candidates=[host_candidate]+auxiliary['candidates']
        selected=min(all_candidates,key=lambda row:row['objective'])
        receipt['handoff_wiring']=auxiliary['handoff_info']
        np.savez_compressed(out/'AUXILIARY_CANDIDATES.npz',
            m=np.array([row['m'] for row in all_candidates]),
            original_Phi=np.array([row['objective'] for row in all_candidates]),
            physical_relative_residual=np.array([row['physical_relative_residual'] for row in all_candidates]))
        np.savez_compressed(out/'SELECTED_BEFORE_CORRECTOR.npz',m=selected['m'])
        receipt['auxiliary_candidate_scores']=[{key:value for key,value in row.items() if key!='m'} for row in all_candidates]
        receipt['auxiliary_selected_stage']=selected['stage']
        receipt['auxiliary_selected_original_Phi']=selected['objective']
        receipt['auxiliary_selected_relative_residual']=selected['physical_relative_residual']
        start,ct,before=time.perf_counter(),time.process_time(),Counter(model.cost)
        refined,corrector=refine(model,Y,selected['m'],maximum=cfg['corrector_maxiter'])
        receipt['corrector_wall_seconds']=time.perf_counter()-start
        receipt['corrector_cpu_seconds']=time.process_time()-ct
        receipt['corrector_physical_model_cost']=dict(Counter(model.cost)-before)
        dump(out/'CORRECTOR.json',corrector)
        receipt['corrector']={key:value for key,value in corrector.items() if key!='history'}
        start,ct,before=time.perf_counter(),time.process_time(),Counter(model.cost)
        refined_phi=model.physical_score(refined,Y)
        best=refined if refined_phi<selected['objective'] else selected['m']
        state=model.state(best);residual=state['Y']-Y
        final_phi=.5*float(np.vdot(residual,residual).real/np.vdot(Y,Y).real)+float(model.reg.evaluate(best)[0])
        np.savez_compressed(out/'prediction.npz',m=best,j=state['j'],Q=state['Q'],Y_prediction=state['Y'])
        receipt['finalization_wall_seconds']=time.perf_counter()-start
        receipt['finalization_cpu_seconds']=time.process_time()-ct
        receipt['finalization_physical_model_cost']=dict(Counter(model.cost)-before)
        receipt.update(status='COMPLETED_DEVELOPMENT_STAGE_HANDOFF',physical_objective=final_phi,
            physical_relative_residual=float(np.linalg.norm(residual)/np.linalg.norm(Y)),
            linear_relative_residual=float(state['state_relative_residual']),
            iterations=sum('iteration' in row for row in auxiliary['history']),
            corrector_improved_original_Phi=bool(refined_phi<selected['objective']),
            selection='minimum original Phi of reused paid host + resumed stage0 material + two remaining stage endpoints + original R4 corrector trajectory incumbent')
except Exception as exc:
    receipt.update(status='FAILED',error=repr(exc),traceback=traceback.format_exc())
    receipt['failed_attempts'].append(dict(error=repr(exc),charged=True))
    if hasattr(exc,'cie_execution_receipt'):
        receipt['failed_auxiliary_execution']=exc.cie_execution_receipt
finally:
    receipt.update(finished_utc=datetime.now(timezone.utc).isoformat(),
        wall_seconds=time.perf_counter()-WALL,cpu_seconds=time.process_time()-CPU,
        physical_model_cost=dict(model.cost) if model else {},
        physical_operation_wall=dict(model.timings) if model else {},
        physical_operation_cpu=dict(model._physics.cpu_timings) if model else {},
        callbacks=observer.receipt() if observer else None,
        cost_scope='Process imports/public+input-stage load/setup/1000 remaining auxiliary updates with original initial+stage physics validation/stage observer/handoff reconstruction/original R4 corrector1000 allowance/final incumbent physics and serialization; paid historical host Phi explicitly reused. Final receipt write excluded; nested times must not be added twice.')
    dump(out/'RECEIPT.json',receipt)
print(json.dumps(dict(handoff_id=args.handoff_id,status=receipt['status'],wall_seconds=receipt['wall_seconds'],cpu_seconds=receipt['cpu_seconds'])),flush=True)
if receipt['status']=='FAILED':raise SystemExit(1)
