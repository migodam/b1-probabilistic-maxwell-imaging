"""Frozen CPU1 timing harness; observation-only, no solver/training/GT access.
Run only after the parent orchestrator confirms the physical worker has ended.
"""
import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
import argparse,fcntl,gc,json,platform,resource,sys,time
from pathlib import Path
sys.dont_write_bytecode=True
import numpy as np
import scipy
import torch
from model_actions import Action,HIST,S,ONLINE_KEYS
from p3r2.residual_risk import read_anchor
from p3r2.risk_models import observation_features
OUT=Path(__file__).resolve().parents[1]
METHODS={'StrongCNN':'strongCNN','K32_optimized':'original_K32_frozen','detMatchedQ64':'deterministic_matched','distillW128':'distilled_w128_n32768'}
SEED=2026102601

def read_inputs(manifest=None):
    if manifest is None:
        with np.load(HIST/'prepared_validation.npz',allow_pickle=False) as z:
            y=z['y'][:64].copy();parents=z['parent'][:64].astype(str)
        online=read_anchor(S/'inputs/strong_context_1024_256/anchor_validation.npz')
        assert np.array_equal(parents,online['parent'][:64])
        online={k:online[k][:64].copy() for k in ONLINE_KEYS}
        return y,online,parents,dict(regime='historical exposed VAL',source=str(HIST/'prepared_validation.npz'),truth_materialized=False)
    cfg=json.loads(Path(manifest).read_text());rows=cfg['parents']
    if len(rows)<64:raise ValueError('at least 64 parent contexts required')
    parents=np.array([r['parent'] for r in rows[:64]],str);online={k:[] for k in ONLINE_KEYS}
    for parent in parents:
        with np.load(OUT/'raw/contexts'/f'{parent}.npz',allow_pickle=False) as z:
            if str(z['parent'])!=parent:raise ValueError('parent mismatch')
            for k in ONLINE_KEYS:online[k].append(z[k].copy())
        with np.load(OUT/'raw/online'/f'{parent}.npz',allow_pickle=False) as z:
            if set(z.files)!={'observation','target_id','sigma'}:raise ValueError('strict raw observation allowlist')
            if not np.array_equal(z['observation'],online['observation'][-1]):raise ValueError('raw/context observation mismatch')
    online={k:np.stack(v) for k,v in online.items()};y=observation_features(online['observation'])
    return y,online,parents,dict(regime='new observation-only context',source=str(Path(manifest).resolve()),truth_materialized=False)

def peak_bytes():
    value=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(value if sys.platform=='darwin' else value*1024)

def distribution(raw,batch):
    p=np.quantile(raw,[.1,.5,.9])
    return dict(raw_seconds=raw,p10_seconds=float(p[0]),p50_seconds=float(p[1]),p90_seconds=float(p[2]),p50_seconds_per_parent=float(p[1]/batch),unit='seconds per whole batch; per-parent value separately divided by batch size')

def elapsed(fn):
    start=time.perf_counter();value=fn();return time.perf_counter()-start,value

def provider(act,y,parents):
    # StrongCNN has no separate posterior-moment provider: no artificial zero cost.
    return act.moments(y,parents)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--parents-manifest',type=Path,default=None)
    ap.add_argument('--repeats',type=int,default=30)
    ap.add_argument('--output',type=Path,default=OUT/'audit/timing/CONFIRMATION_TIMING.json')
    ap.add_argument('--physics-ended',action='store_true',help='orchestrator must have verified physical-worker completion')
    args=ap.parse_args()
    if not args.physics_ended:raise RuntimeError('parent verification of physical-worker completion is required')
    if args.repeats<30:raise ValueError('at least 30 repeats required')
    # Fail immediately if a physical worker is active; hold its existing lock for
    # the entire timing run, preventing an accidental concurrent physical launch.
    lockpath=OUT/'raw/physics_worker.lock';lockpath.parent.mkdir(parents=True,exist_ok=True)
    with lockpath.open('a') as lock:
        try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:raise RuntimeError('physical worker is active; timing prohibited')
        torch.set_num_threads(1);torch.set_num_interop_threads(1)
        y,online,parents,input_scope=read_inputs(args.parents_manifest)
        for k,v in online.items():
            if len(v)!=64 or not np.isfinite(v).all():raise ValueError('invalid online data')
        expected={'observation':(64,32,8),'residual':(64,32,8),'q0':(64,66),'chi0':(64,576),'J':(64,512,66)}
        if any(online[k].shape!=s for k,s in expected.items()):raise ValueError('online shape contract')
        initial_peak=peak_bytes();records=[];start=time.perf_counter()
        for batch in (1,64):
            yy=y[:batch];pp=parents[:batch];aa={k:v[:batch] for k,v in online.items()}
            cold={label:dict(load=[],first=[],combined=[],memory=[]) for label in METHODS}
            # Interleave arms within each independent Action-load round. Imports
            # and OS filesystem caches stay warm; no disk-cache purge is claimed.
            for repeat in range(args.repeats):
                labels=list(METHODS);labels=labels[repeat%4:]+labels[:repeat%4]
                for label in labels:
                    gc.collect();load,act=elapsed(lambda:Action(METHODS[label],SEED,device='cpu',quadrature_nodes=64))
                    first,pred=elapsed(lambda:act.predict(yy,aa,pp,batch_size=batch))
                    if pred.shape!=(batch,576) or not np.isfinite(pred).all():raise ValueError('invalid first prediction')
                    row=cold[label];row['load'].append(load);row['first'].append(first);row['combined'].append(load+first);row['memory'].append(peak_bytes());del act,pred
            actions={label:Action(method,SEED,device='cpu',quadrature_nodes=64) for label,method in METHODS.items()}
            steady={label:dict(full=[],provider=[]) for label in METHODS}
            for label,act in actions.items():
                act.predict(yy,aa,pp,batch_size=batch)
                if label!='StrongCNN':provider(act,yy,pp)
            for repeat in range(args.repeats):
                labels=list(METHODS);labels=labels[repeat%4:]+labels[:repeat%4]
                for label in labels:
                    act=actions[label];calls=['full'] if label=='StrongCNN' else (['provider','full'] if repeat%2==0 else ['full','provider'])
                    for kind in calls:
                        duration,value=elapsed((lambda:act.predict(yy,aa,pp,batch_size=batch)) if kind=='full' else (lambda:provider(act,yy,pp)))
                        if not np.isfinite(value).all():raise ValueError('nonfinite inference')
                        steady[label][kind].append(duration)
            for label in METHODS:
                row=dict(arm=label,method=METHODS[label],seed=SEED,batch_size=batch,repeats=args.repeats,parents=pp.tolist(),fresh_action_load=distribution(cold[label]['load'],batch),first_full_inference_after_load=distribution(cold[label]['first'],batch),load_plus_first_inference=distribution(cold[label]['combined'],batch),steady_full_predict=distribution(steady[label]['full'],batch),cumulative_process_peak_RSS_bytes_per_load=cold[label]['memory'],cumulative_process_peak_RSS_bytes=peak_bytes())
                row['steady_provider_only']=dict(status='NOT_APPLICABLE',reason='StrongCNN has no separate posterior-moment provider') if label=='StrongCNN' else dict(status='MEASURED',**distribution(steady[label]['provider'],batch))
                records.append(row)
            del actions;gc.collect()
        res=dict(status='COMPLETE',canonical_backend='CPU1',seed=SEED,seed_scope='single canonical frozen training seed; 30 timing repeats are repeated execution, not independent training-seed evidence; other two seeds have separate functional checks',input_scope=input_scope,records=records,physics_RHS=0,new_GT_read=False,training_updates=0,process_peak_RSS_bytes=peak_bytes(),initial_process_peak_RSS_bytes=initial_peak,wall_seconds=time.perf_counter()-start,optional_device=dict(MPS='AVAILABLE_NOT_RUN' if torch.backends.mps.is_available() else 'NOT_AVAILABLE',CUDA='AVAILABLE_NOT_RUN' if torch.cuda.is_available() else 'NOT_AVAILABLE',canonical_backend_unchanged=True),environment=dict(platform=platform.platform(),python=sys.version,numpy=np.__version__,scipy=scipy.__version__,torch=torch.__version__,torch_threads=torch.get_num_threads(),torch_interop_threads=torch.get_num_interop_threads(),numeric_environment={k:os.environ.get(k) for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS')}),scope_notes=['Fresh Action load timings include checkpoint reads, construction, loading normalization/weights, and strict legacy metadata audit; imports and filesystem caches are warm, so not true process/disk cold-start.','First full inference follows each of 30 independent Action loads. Steady inference follows one explicit warmup per arm.','Full predict includes measurement validation, moment provider where present, frozen online feature construction, head, and material packing; excludes physical acquisition and Native/Probit anchor construction.','Provider only includes full four-channel moment representation and K32 parent-local RNG/decoding; excludes feature construction and risk head.','Process peak RSS is cumulative across the harness, not an isolated per-arm allocation peak.','Historical teacher generation and all training costs remain read-only separate sunk-cost records. Full physics and end-to-end aggregated parent runtime must be reported separately from these inference timings.','All arms use identical first parent or first64 parents, observations, residuals, q0, chi0, J, and parent IDs. New raw-context loading materializes only allowed online fields.','Optimized K32 exact CPU prediction equivalence to unoptimized historical sampling and 2e-6 stored-moment gate are already recorded in CHECKPOINT_REPRODUCTION.json.'])
        args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(dict(status=res['status'],output=str(args.output),wall_seconds=res['wall_seconds'])))
if __name__=='__main__':main()
