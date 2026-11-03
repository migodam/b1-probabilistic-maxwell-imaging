"""Analyze completed frozen replay runs sequentially; public inputs only."""
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[name]='1'
from pathlib import Path
import argparse
import csv
import json
import sys
import time
import traceback
from datetime import datetime, timezone
ROOT=Path(__file__).resolve().parents[1]
REPO=ROOT.parents[1]
sys.path.insert(0,str(REPO))
import numpy as np
from threadpoolctl import threadpool_limits
from research.research_fbg_r3_experiments.src.physics import MaxwellModel
from research.research_max_struct_r5.src.trajectory_diagnostics import analyze_checkpoint


def sanitize(value):
    if isinstance(value,dict):return {k:sanitize(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [sanitize(x) for x in value]
    if isinstance(value,np.generic):value=value.item()
    if isinstance(value,float) and not np.isfinite(value):return None
    return value


def dump(path,value):
    Path(path).write_text(json.dumps(sanitize(value),indent=2,allow_nan=False))


def load_public(cfg):
    public={}
    for name in cfg['public_files']:
        with np.load(REPO/name,allow_pickle=False) as data:
            if any(k.lower() in {'truth','alpha','mask','yclean','seed','family'} for k in data.files):
                raise ValueError('Unexpected private payload key')
            public.update({k:data[k] for k in data.files})
    return public


def analyze_run(cfg):
    source=ROOT/'RESULTS/trajectory_replay'/cfg['run_id']
    ready=json.loads((source/'RECEIPT.json').read_text())
    if ready['status']!='COMPLETED_AUX_ONLY':raise RuntimeError('Replay is not complete')
    paths=sorted(source.glob('checkpoint_s*_i*.npz'))
    if len(paths)!=21:raise ValueError('Expected exactly 21 frozen checkpoints')
    out=ROOT/'RESULTS/trajectory_diagnostics'/cfg['run_id']
    out.mkdir(parents=True,exist_ok=True)
    previous={}
    if (out/'RECEIPT.json').exists():
        previous=json.loads((out/'RECEIPT.json').read_text())
        if previous.get('status')=='COMPLETED':
            print(json.dumps(dict(event='ALREADY_COMPLETED',run_id=cfg['run_id'])),flush=True);return
    wall,cpu=time.perf_counter(),time.process_time()
    receipt=dict(status='RUNNING',run_id=cfg['run_id'],truth_access='NONE',threads=1,device='CPU',
        started_utc=datetime.now(timezone.utc).isoformat(),failed_attempts=previous.get('failed_attempts',[]),checkpoints=[],
        prior_attempts=previous.get('prior_attempts',[])+([previous] if previous else []),
        configuration=cfg,corrector='NOT_RUN',parameter_tuning=False,seed=None)
    dump(out/'RECEIPT.json',receipt)
    model=None
    try:
        public=load_public(cfg);model=MaxwellModel(public);Y=public['Y']
        rows=[]
        for path in paths:
            if (out/(path.stem+'.json')).exists():
                summary=json.loads((out/(path.stem+'.json')).read_text())
                if summary.get('truth_access')!='NONE':raise ValueError('Invalid prior output')
            else:
                with np.load(path,allow_pickle=False) as data:
                    cp={k:data[k] for k in data.files}
                stage=cfg['numerical']['stages'][int(cp['stage'])]
                beta=float(stage.get('beta_material',stage.get('beta',1.)))
                summary,arrays=analyze_checkpoint(model,Y,cp,beta,cfg['numerical'],endpoint=int(cp['iteration'])==500)
                np.savez_compressed(out/(path.stem+'.npz'),**arrays)
                dump(out/(path.stem+'.json'),summary)
            receipt['checkpoints'].append(path.stem)
            row=dict(checkpoint=path.stem,stage=summary['stage'],iteration=summary['iteration'],
                beta=summary['beta_material'],cross_value=summary['cross_value'],**summary['physical'])
            for kind in ('ordinary','self','plus','minus'):
                ls=summary['frozen_LS'][kind]
                for key in ('norm','cosine_actual_step','cosine_negative_physical_gradient'):
                    row['transport_'+kind+'_'+key]=ls['transport'][key]
            for kind in ('self','plus','minus'):
                for component in ('exact_state','frozen','connection'):
                    for key in ('norm','cosine_actual_step','cosine_negative_physical_gradient'):
                        row['gradient_'+kind+'_'+component+'_'+key]=summary['gradients'][kind][component][key]
            rows.append(row)
            print(json.dumps(dict(event='CHECKPOINT_ANALYZED',run_id=cfg['run_id'],checkpoint=path.stem,Phi=row['Phi_before'])),flush=True)
        with (out/'CHECKPOINTS.csv').open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
        receipt['status']='COMPLETED'
    except Exception as exc:
        receipt.update(status='FAILED',error=repr(exc),traceback=traceback.format_exc())
        receipt['failed_attempts'].append(dict(error=repr(exc),charged=True))
        raise
    finally:
        receipt.update(finished_utc=datetime.now(timezone.utc).isoformat(),
            wall_seconds=time.perf_counter()-wall,cpu_seconds=time.process_time()-cpu,
            physical_model_cost=dict(model.cost) if model else {},
            physical_operation_wall=dict(model.timings) if model else {},
            physical_operation_cpu=dict(model._physics.cpu_timings) if model else {},money=None,
            cost_scope='this process run setup/checkpoint diagnostics/endpoints/counterfactual forwards/output; per-checkpoint auxiliary/physics counts in JSON; existing checkpoint files reused if resumed and those historic costs remain in per-checkpoint JSON')
        dump(out/'RECEIPT.json',receipt)


def analyze_legacy_states():
    manifest=json.loads((ROOT/'configs/trajectory_replay_frozen.json').read_text())
    public=load_public(manifest['runs'][0]);Y=public['Y']
    for arm in ('self_absorption','onehop_phase_pi'):
        source=REPO/'research/research_max_struct_r4/results'/arm
        config=json.loads((source/'FROZEN_CONFIG.json').read_text())
        out=ROOT/'RESULTS/trajectory_diagnostics'/('legacy_'+arm)
        out.mkdir(exist_ok=True,parents=True)
        wall,cpu=time.perf_counter(),time.process_time()
        model=MaxwellModel(public)
        receipt=dict(status='RUNNING',truth_access='NONE',actual_step_available=False,
            source=str(source.relative_to(REPO)),device='CPU',threads=1,stages=[],
            corrector='NOT_RUN',parameter_tuning=False,failed_attempts=[],
            configuration=config)
        try:
            for stage in range(3):
                path=source/f'stage_{stage}.npz'
                with np.load(path,allow_pickle=False) as data:
                    state={key:data[key] for key in data.files}
                cp=dict(m_state=state['m'],J_after_current=state['J'],Q=state['Q'],
                        stage=stage,iteration=500)
                stageconfig=config['stages'][stage]
                beta=float(stageconfig.get('beta_material',stageconfig.get('beta',1.)))
                summary,arrays=analyze_checkpoint(model,Y,cp,beta,config,endpoint=True)
                summary['historical_source']=str(path.relative_to(REPO))
                np.savez_compressed(out/f'stage_{stage}.npz',**arrays)
                dump(out/f'stage_{stage}.json',summary)
                receipt['stages'].append(stage)
                print(json.dumps(dict(event='LEGACY_STATE_ANALYZED',arm=arm,stage=stage)),flush=True)
            receipt['status']='COMPLETED_FIXED_STATES'
        except Exception as exc:
            receipt.update(status='FAILED',error=repr(exc),traceback=traceback.format_exc())
            receipt['failed_attempts'].append(dict(error=repr(exc),charged=True))
            raise
        finally:
            receipt.update(wall_seconds=time.perf_counter()-wall,cpu_seconds=time.process_time()-cpu,
                physical_model_cost=dict(model.cost),physical_operation_wall=dict(model.timings),
                physical_operation_cpu=dict(model._physics.cpu_timings),money=None,
                cost_scope='6 historical fixed-state diagnostic reuse, not auxiliary replay; no actual historical step, no HVP along unavailable step; all new forward/adjoint/Jacobian/readout costs charged')
            dump(out/'RECEIPT.json',receipt)


if __name__=='__main__':
    parser=argparse.ArgumentParser();group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--run-id');group.add_argument('--all-ready',action='store_true');group.add_argument('--legacy-stage-states',action='store_true')
    args=parser.parse_args()
    manifest=json.loads((ROOT/'configs/trajectory_replay_frozen.json').read_text())
    with threadpool_limits(limits=1):
        if args.legacy_stage_states:
            analyze_legacy_states();raise SystemExit(0)
        for cfg in manifest['runs']:
            if args.run_id and cfg['run_id']!=args.run_id:continue
            path=ROOT/'RESULTS/trajectory_replay'/cfg['run_id']/'RECEIPT.json'
            if not path.exists() or json.loads(path.read_text()).get('status')!='COMPLETED_AUX_ONLY':
                if args.run_id:raise RuntimeError('Requested replay not ready')
                continue
            analyze_run(cfg)
