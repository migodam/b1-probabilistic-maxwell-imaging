"""Frozen CPU1 four-arm inversion; opens allowlisted PUBLIC payloads only.

No generator import, no GT, no score, no parameter tuning and no hash checks.
"""
from pathlib import Path
from collections import Counter
from contextlib import contextmanager
from datetime import datetime,timezone
import argparse,json,sys,time,traceback
import numpy as np
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
sys.path.insert(0,str(REPO))
from research.research_fbg_r3_experiments.src.physics import MaxwellModel
from research.research_fbg_r3_experiments.src import cie_tsom_baseline as native
from research.research_fbg_r3_experiments.src.full_maxwell_corrector import corrector
from research.research_max_struct_r2.src.cie_global import run_cie as global_cie
from research.research_max_struct_r2.src.profile_current import run_profile

STAGES=[dict(mf=3,beta=1.,maxiter=500),
        dict(mf=6,beta_material=.2,beta_current=1.,maxiter=500),
        dict(mf=10,beta_material=.2,beta_current=1.,maxiter=500)]
ARMS=['native','global','direct','profile']
ALLOW={'G','E0','S','Y','grid_shape','shape','sigma','noise_radius','host','amplitude',
       'material_bounds','smooth_weight','tv_weight','tv_epsilon','prior_width','coords',
       'sources','receivers','k','radius'}

def utc():return datetime.now(timezone.utc).isoformat()

def clean(x):
    if isinstance(x,dict):return {k:clean(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [clean(v) for v in x]
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,np.generic):return x.item()
    return x

def dump(path,obj):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(clean(obj),indent=2,allow_nan=False))

def compact(result):
    return {k:v for k,v in result.items() if k not in ('m','J','stage_states','last_joint_m')
            and k!='candidates'}|{'candidates':[{k:v for k,v in c.items()if k!='m'}for c in result.get('candidates',[])]}

def save_stage(path,state):
    payload={k:state[k] for k in ('m','J','Q','coefficients') if k in state}
    np.savez_compressed(path,**payload)

@contextmanager
def native_endpoints(directory):
    """Forwarding observation wrapper only; mature numerical code unchanged."""
    original=native._readout;count=[0]
    def readout(J,Q,*args,**kwargs):
        m,meta=original(J,Q,*args,**kwargs);count[0]+=1
        if count[0]%500==0:
            save_stage(directory/f'stage_{count[0]//500-1}.npz',dict(m=m,J=J,Q=Q))
        return m,meta
    native._readout=readout
    try:yield count
    finally:native._readout=original

def metric(model,Y,m):
    state=model.state(m);r=state['Y']-Y;rv=model.reg.evaluate(m)[0]
    return dict(m=np.asarray(m).copy(),objective=float(.5*np.vdot(r,r).real/np.vdot(Y,Y).real+rv),
                relative_residual=float(np.linalg.norm(r)/np.linalg.norm(Y)))

def gradient_check():
    path=ROOT/'audit/PROFILE_GRADIENT_WIRING.json'
    if path.exists():return
    start=time.perf_counter();cpu=time.process_time();model=None
    receipt=dict(status='RUNNING',case='existing_R1_4x4_public_operators',new_directions=1,
                 truth_access='DISCLOSED_CONSTANT_DEVELOPMENT_ROOT_ONLY',hash_checks=False,
                 code='src/profile_current.py::run_profile(gradient_check=True)')
    try:
        with np.load(REPO/'research/research_max_struct_r1/data_public/dynamics_operators_4x4.npz',allow_pickle=False) as f:
            ops={k:f[k] for k in ('G','E0','S','grid_shape')}
        model=MaxwellModel(ops)
        Y=model.state(np.full(model.n,.55))['Y'].copy()
        result=run_profile(model,Y,np.full(model.n,.1),gradient_check=True)
        receipt.update(status='COMPLETED',result=result,validation='one finite-direction wiring check; no new identity suite')
    except Exception as exc:
        receipt.update(status='FAILED',error=repr(exc),traceback=traceback.format_exc());raise
    finally:
        receipt.update(wall_seconds=time.perf_counter()-start,cpu_seconds=time.process_time()-cpu,
                       physical_model_cost=dict(model.cost)if model else {},
                       timer_scope='operator load, disclosed-root forward,3 profiled-current solves; excludes imports/receipt write')
        dump(path,receipt)
    print(json.dumps(dict(event='PROFILE_GRADIENT_WIRING',status=receipt['status'],result=receipt['result'])),flush=True)

def run_one(record,arm):
    ident=record['observation_id'];directory=ROOT/'results/blind'/ident/arm
    receipt_path=directory/'RECEIPT.json';prediction=ROOT/'predictions'/ident/f'{arm}.npz'
    if receipt_path.exists() or prediction.exists():
        print(json.dumps(dict(event='EXISTING_DISPOSITION_PRESERVED',observation_id=ident,arm=arm)),flush=True);return
    directory.mkdir(parents=True,exist_ok=True);prediction.parent.mkdir(parents=True,exist_ok=True)
    wall=time.perf_counter();cpu=time.process_time();model=None
    receipt=dict(status='RUNNING',run_id=ident+'_'+arm,observation_id=ident,arm=arm,started_utc=utc(),
        classification='FRESH_BLIND_PREDICTION_NO_GT',truth_access='NONE',hash_checks=False,
        source='experiments/run_blind.py',public_data='data_public/'+record['path'],
        physical_model='unchanged scalar-TM disk-pulse VIE,complex128,576 independent real cells',
        source_metric_policy='fixed original paper source normalizations for native/global; uniform totalenergy for known profile control',
        initialization='all cells public host=.1',material_box=[.02,20.],prior='unchanged historical smooth1e-4/TV1e-4',
        stages=STAGES if arm in ('native','global')else None,coefficient_tol=0,miniter=50,
        iteration_policy='native/global full1500; direct up to1000 L-BFGS; profile up to250',
        corrector=dict(method='original physical L-BFGS-B',maxiter=1000,maxfun=4000,wallcap=None),
        online_selection='original normalized physical objective+same prior; always include common initial material',
        threads=1,device='CPU',failed_attempts=[],operator_assembly='frozen public load; no reassembly')
    dump(receipt_path,receipt)
    try:
        with np.load(ROOT/'data_public'/record['path'],allow_pickle=False) as f:
            if not set(f.files)<=ALLOW:raise ValueError('Unexpected public key, refuse payload')
            public={k:f[k]for k in f.files}
        model=MaxwellModel(public);Y=public['Y'];m0=np.full(model.n,float(public['host']))
        initial=metric(model,Y,m0);candidates=[initial]
        cfg=dict(stages=STAGES,method='cie',dpcs_rtol=1e-3,coefficient_tol=0.,miniter=50,
                 material_inner_maxiter=60,material_box=(.02,20.),regularizer=model.regularizer,
                 physical_sigma=float(np.linalg.norm(Y)),sourcewise_cg=(arm=='native'))
        if arm=='native':
            with native_endpoints(directory) as calls:aux=native.run_cie(model,Y,m0,cfg)
            receipt['native_forwarding_readout_calls']=calls[0]
            candidates.extend(aux['candidates'])
        elif arm=='global':
            cfg['stage_callback']=lambda st:save_stage(directory/f"stage_{st['stage']}.npz",st)
            aux=global_cie(model,Y,m0,cfg);candidates.extend(aux['candidates'])
        elif arm=='direct':
            aux=corrector(model,Y,m0,dict(method='lbfgs',max_nfev=1000,box=(.02,20.)))
            candidates.append(metric(model,Y,aux['m']))
        elif arm=='profile':
            aux=run_profile(model,Y,m0,maxiter=250,beta=1.)
            candidates.append(metric(model,Y,aux['m']))
            save_stage(directory/'profile_endpoint.npz',dict(m=aux['m'],J=aux['J'],Q=model.E+model.G@aux['J']))
        else:raise ValueError(arm)
        # Candidate keys differ in historical corrector; all ranking is original Phi.
        scored=[metric(model,Y,c['m'])for c in candidates]
        seed=min(scored,key=lambda c:c['objective'])
        receipt['auxiliary_selected_physical_objective']=seed['objective']
        receipt['auxiliary']=compact(aux)
        dump(directory/'AUXILIARY.json',receipt['auxiliary'])
        before_corrector=dict(model.cost)
        corr=corrector(model,Y,seed['m'],dict(method='lbfgs',max_nfev=1000,box=(.02,20.)))
        final_candidates=[initial,seed,metric(model,Y,corr['m'])]
        best=min(final_candidates,key=lambda c:c['objective']);state=model.state(best['m'])
        residual=state['Y']-Y;sigma=float(public['sigma'])
        np.savez_compressed(prediction,m=best['m'],Y_prediction=state['Y'],j=state['j'],Q=state['Q'])
        dump(directory/'CORRECTOR.json',compact(corr))
        receipt.update(status='SEALED_PREDICTION',finished_utc=utc(),prediction=str(prediction.relative_to(ROOT)),
            physical_objective=best['objective'],physical_relative_residual=best['relative_residual'],
            physical_whitened_residual=float(np.linalg.norm(residual)/sigma),
            noise_radius99=float(public['noise_radius']),linear_relative_residual=state['state_relative_residual'],
            material_min=float(best['m'].min()),material_max=float(best['m'].max()),
            corrector_result=compact(corr),corrector_model_cost_before=before_corrector,
            post_auxiliary_material_candidate_count=len(scored))
    except Exception as exc:
        receipt.update(status='FAILED',finished_utc=utc(),error=repr(exc),traceback=traceback.format_exc())
        if hasattr(exc,'cie_execution_receipt'):receipt['failed_cie_cost']=exc.cie_execution_receipt
        receipt['failed_attempts'].append(dict(error=repr(exc),included_in_total=True))
    finally:
        receipt.update(wall_seconds=time.perf_counter()-wall,cpu_seconds=time.process_time()-cpu,
            physical_model_cost=dict(model.cost)if model else {},
            physical_operation_wall=dict(model.timings)if model else {},
            physical_operation_cpu=dict(model._physics.cpu_timings)if model else {},
            cost_scope='load/preparation/auxiliary/physical candidate reranking/corrector/prediction+stage serialization; final receipt write excluded',
            shared_historical_operator_assembly_cost=None,
            historical_cost_reason='frozen operator assembly not rerun; generation separately accounted')
        dump(receipt_path,receipt)
    print(json.dumps({k:receipt.get(k)for k in ['status','run_id','wall_seconds','cpu_seconds','physical_relative_residual','physical_whitened_residual','noise_radius99','physical_model_cost','error']}),flush=True)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--case');parser.add_argument('--arm',choices=ARMS)
    parser.add_argument('--gradient-check-only',action='store_true');args=parser.parse_args()
    with threadpool_limits(limits=1):
        gradient_check()
        if args.gradient_check_only:return
        manifest=json.loads((ROOT/'data_public/manifest.json').read_text())
        for record in manifest['records']:
            if args.case and args.case!=record['observation_id']:continue
            for arm in ([args.arm]if args.arm else ARMS):run_one(record,arm)

if __name__=='__main__':main()
