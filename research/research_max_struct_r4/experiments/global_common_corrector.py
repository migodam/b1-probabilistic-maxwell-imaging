"""Common corrector from sealed R2 global auxiliary seed; never rerun CIE."""
from pathlib import Path
from datetime import datetime, timezone
import json, sys, time, traceback
START_WALL, START_CPU = time.perf_counter(), time.process_time()
ROOT=Path(__file__).resolve().parents[1]; REPO=ROOT.parents[1]
OUT=ROOT/'results/global_common_corrector'
OUT.mkdir(parents=True,exist_ok=True)
if (OUT/'RECEIPT.json').exists():
    raise FileExistsError('Preserve completed common-corrector disposition')
startup=dict(status='RUNNING',executable=sys.executable,argv=sys.argv,
             started_utc=datetime.now(timezone.utc).isoformat())
(OUT/'STARTUP.json').write_text(json.dumps(startup,indent=2))
try:
    import numpy as np
    from threadpoolctl import threadpool_limits
    sys.path.insert(0,str(REPO))
    from research.research_fbg_r3_experiments.src.physics import MaxwellModel
    from research.research_max_struct_r2.experiments.run_blind import ALLOW,dump,metric
    from research.research_max_struct_r4.experiments.cooperative_paths import refine
    startup.update(status='IMPORTED',wall_seconds=time.perf_counter()-START_WALL,
                   cpu_seconds=time.process_time()-START_CPU)
except Exception as exc:
    startup.update(status='FAILED',error=repr(exc),traceback=traceback.format_exc(),
                   wall_seconds=time.perf_counter()-START_WALL,cpu_seconds=time.process_time()-START_CPU)
    (OUT/'STARTUP.json').write_text(json.dumps(startup,indent=2))
    raise
(OUT/'STARTUP.json').write_text(json.dumps(startup,indent=2))


def run():
    cfg=json.loads((ROOT/'configs/global_common_corrector.json').read_text())
    dump(OUT/'FROZEN_CONFIG.json',cfg)
    rec=dict(status='RUNNING',configuration=cfg,truth_access='NONE',startup=startup,
        classification='GT_DISCLOSED_DEVELOPMENT_CASE_NO_TRUTH_ACCESS_BY_RUNNER',
        auxiliary_rerun=False,failed_attempts=[],historical_cost_recharged=False,
        compared_to='block_absorption uses exactly the same refine function and corrector maxiter')
    dump(OUT/'RECEIPT.json',rec);model=None
    try:
        r2=REPO/'research/research_max_struct_r2'
        old_dir=r2/'results/blind'/cfg['case_id']/'global'
        old=json.loads((old_dir/'RECEIPT.json').read_text())
        aux=json.loads((old_dir/'AUXILIARY.json').read_text())
        selected=min(aux['candidates'],key=lambda x:x['objective'])
        with np.load(r2/'data_public'/f"{cfg['case_id']}.npz",allow_pickle=False) as f:
            if not set(f.files)<=ALLOW:raise ValueError('Unexpected public key')
            public={k:f[k] for k in f.files}
        model=MaxwellModel(public);Y=public['Y']
        if selected['stage']==-1:
            seed=np.full(model.n,float(public['host']))
        else:
            with np.load(old_dir/f"stage_{selected['stage']}.npz",allow_pickle=False) as f:
                seed=f['m'].copy()
        with np.load(r2/'predictions'/cfg['case_id']/'global.npz',allow_pickle=False) as f:
            incumbent=f['m'].copy()
        raw=metric(model,Y,seed);historical=metric(model,Y,incumbent)
        rec.update(selected_auxiliary_candidate=selected,
            seed_loaded_from=f"research/research_max_struct_r2/results/blind/{cfg['case_id']}/global/stage_{selected['stage']}.npz",
            seed_current_original_Phi=raw['objective'],historical_final_original_Phi=historical['objective'],
            historical_global_receipt=str((old_dir/'RECEIPT.json').relative_to(REPO)),
            historical_global_physical_cost=old['physical_model_cost'],
            historical_global_wall_seconds=old['wall_seconds'],
            historical_auxiliary_execution=aux.get('execution_receipt'),
            comparability='same saved mature auxiliary algorithm and common1000 scaled corrector; historical published incumbent also retained; block inner work differs')
        dump(OUT/'RECEIPT.json',rec)
        np.savez_compressed(OUT/'SELECTED_BEFORE_CORRECTOR.npz',m=seed)
        final,corr=refine(model,Y,seed,maximum=cfg['corrector_maxiter'])
        dump(OUT/'CORRECTOR.json',corr)
        scores=[raw,historical,metric(model,Y,final)]
        selected_index=int(np.argmin([x['objective'] for x in scores]))
        best=scores[selected_index];state=model.state(best['m'])
        np.savez_compressed(OUT/'prediction.npz',m=best['m'],Y_prediction=state['Y'],j=state['j'],Q=state['Q'])
        rec.update(status='SEALED_DEVELOPMENT_PREDICTION',physical_objective=best['objective'],
            physical_relative_residual=best['relative_residual'],
            physical_whitened_residual=float(np.linalg.norm(state['Y']-Y)/float(public['sigma'])),
            physical_linear_residual=float(state['state_relative_residual']),
            material_min=float(best['m'].min()),material_max=float(best['m'].max()),
            corrector_stop=corr['stop'],corrector_iterations=corr['iterations'],
            candidate_scores=[{k:v for k,v in x.items() if k!='m'} for x in scores],
            selected_candidate=['saved_auxiliary_seed','published_global_incumbent','common_corrector'][selected_index],
            prediction='results/global_common_corrector/prediction.npz')
    except Exception as exc:
        rec.update(status='FAILED',error=repr(exc),traceback=traceback.format_exc())
        rec['failed_attempts'].append(dict(error=repr(exc),included_in_total=True));raise
    finally:
        rec.update(wall_seconds=time.perf_counter()-START_WALL,cpu_seconds=time.process_time()-START_CPU,
            physical_model_cost=dict(model.cost) if model else {},
            physical_operation_wall=dict(model.timings) if model else {},
            physical_operation_cpu=dict(model._physics.cpu_timings) if model else {},
            finished_utc=datetime.now(timezone.utc).isoformat(),
            cost_scope='new startup/import/load/ranking/commoncorrector/prediction only; historical auxiliary separately identified and not rerun')
        dump(OUT/'RECEIPT.json',rec)
    print(json.dumps({k:rec[k] for k in ['status','wall_seconds','cpu_seconds','physical_objective',
        'physical_relative_residual','physical_model_cost']}),flush=True)


if __name__=='__main__':
    with threadpool_limits(limits=1):run()
