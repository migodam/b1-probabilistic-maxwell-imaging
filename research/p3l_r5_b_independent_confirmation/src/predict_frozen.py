"""Seal canonical frozen images before any new truth is evaluated.

Reads original online observations and truth-free native/probit contexts only.
An explicit key allowlist, matched observation identity and checkpoint/source
identities are checked. Never imports the generator or offline evaluator.
"""
import os
for _k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS'):
    os.environ[_k]='1'
from pathlib import Path
import sys,json,time,hashlib,resource,argparse,gc
sys.dont_write_bytecode=True
import numpy as np,torch
from model_actions import Action,ONLINE_KEYS,SEEDS,R4
from p3r2.risk_models import observation_features
from p3r2.fusion import passive_box

OUT=Path(__file__).resolve().parents[1]
METHODS={'physics_probit':'physics','strong_cnn':'strongCNN','k32':'original_K32_frozen',
         'deterministic_matched':'deterministic_matched','distilled_w128':'distilled_w128_n32768',
         'det_frozen':'det_frozen'}
def identity(path):
    path=Path(path);return dict(path=str(path.resolve()),bytes=path.stat().st_size,
                               blake2b=hashlib.blake2b(path.read_bytes()).hexdigest())
def write(path,obj):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n');tmp.replace(path)
def save(path,**arrays):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('xb') as f:np.savez_compressed(f,**arrays)
def check_frozen():
    for row in json.loads((OUT/'configs/FREEZE.json').read_text())['files']:
        if identity(row['path'])['blake2b']!=row['blake2b']:raise RuntimeError('frozen source/config changed: '+row['path'])
    manifest=json.loads((OUT/'audit/model_replay/CHECKPOINT_IDENTITIES_BLAKE2B.json').read_text())
    for row in manifest['checkpoints']:
        if identity(row['path'])['blake2b']!=row['blake2b_512']:
            raise RuntimeError('frozen checkpoint changed: '+row['path'])
    return manifest
def read_batch(parents):
    rows=[]
    for parent in parents:
        with np.load(OUT/'raw/contexts'/f'{parent}.npz',allow_pickle=False) as z:
            if str(z['parent'])!=parent:raise ValueError('parent/context identity')
            a={k:z[k].copy() for k in ONLINE_KEYS}
        with np.load(OUT/'raw/online'/f'{parent}.npz',allow_pickle=False) as z:
            if set(z.files)!={'observation','sigma','target_id'}:raise ValueError('observation key allowlist')
            if not np.array_equal(a['observation'],z['observation']):raise ValueError('context observation mismatch')
        rows.append(a)
    return {k:np.stack([r[k] for r in rows]) for k in ONLINE_KEYS}
def predict_group(group):
    torch.set_num_threads(1);check_frozen()
    if not (OUT/f'raw/ANCHORS_COMPLETE_{group}.json').exists():raise RuntimeError('all full-fidelity anchors required')
    parents=[r['parent'] for r in json.loads((OUT/f'configs/parents_{group}.json').read_text())['parents']]
    dest=OUT/'raw/predictions'/group;dest.mkdir(parents=True,exist_ok=True)
    if (dest/'ONLINE_COMPLETE.json').exists():return
    artifacts=[];costs=[];failures=[];tick=time.perf_counter();cpu=time.process_time()
    for method,label in METHODS.items():
        for seed in ([SEEDS[0]] if method=='physics_probit' else SEEDS):
            file=dest/f'{method}_s{seed}.npz'
            if file.exists():
                # Only exact, completed batch archives may be recovered.
                with np.load(file,allow_pickle=False) as z:
                    if not np.array_equal(z['parent_ids'],parents) or z['prediction_raw'].shape!=(256,576):raise ValueError('partial prediction archive')
                artifacts.append(identity(file));continue
            start=time.perf_counter();a=Action(label,seed,device='cpu');load=time.perf_counter()-start
            pred=[];clipped=[];fallback=[];call_seconds=0.;input_seconds=0.
            geometry=getattr(a,'geometry',None)
            if geometry is None:
                reference=Action('strongCNN',seed,device='cpu');geometry=reference.geometry;del reference
            for pos in range(0,len(parents),64):
                ids=parents[pos:pos+64];input_start=time.perf_counter();online=read_batch(ids);y=observation_features(online['observation']);input_seconds+=time.perf_counter()-input_start
                start=time.perf_counter();flag=np.zeros(len(ids),bool)
                try:
                    result=a.predict(y,online,parents=ids,batch_size=64)
                    if result.shape!=(len(ids),576) or not np.iscomplexobj(result):raise ValueError('prediction shape/complex contract')
                    flag=~np.isfinite(result).all(1)
                    if flag.any():failures.append(dict(method=method,seed=seed,parents=np.asarray(ids)[flag].tolist(),error='nonfinite image; frozen anchor fallback'))
                    result[flag]=online['chi0'][flag]
                except Exception as exc:
                    result=online['chi0'].copy();flag[:]=True;failures.append(dict(method=method,seed=seed,parents=ids,error=repr(exc),fallback='same online anchor, no GT choice'))
                call_seconds+=time.perf_counter()-start
                pred.append(result.copy());clipped.append(passive_box(result,geometry['config'],geometry['roi']));fallback.append(flag)
            save(file,parent_ids=np.asarray(parents),prediction_raw=np.concatenate(pred),prediction_clipped=np.concatenate(clipped),fallback_anchor=np.concatenate(fallback))
            artifacts.append(identity(file));costs.append(dict(method=method,seed=int(seed),checkpoint_load_seconds=load,full_inference_seconds=call_seconds,context_IO_seconds=input_seconds))
            del a,pred,clipped,online;gc.collect()
            print(json.dumps(dict(group=group,method=method,seed=int(seed),status='PREDICTION_PERSISTED_NO_GT')),flush=True)
    check_frozen()
    write(dest/'ONLINE_COMPLETE.json',dict(status='ALL_FROZEN_PREDICTIONS_SEALED',group=group,parents=parents,artifacts=artifacts,
          utc=time.time(),truth_read=False,physics_RHS=0,training_updates=0,canonical='CPU1,batch64',
          models=METHODS,seed_slots=list(SEEDS),failures=failures,costs=costs,wall_seconds=time.perf_counter()-tick,CPU_seconds=time.process_time()-cpu,
          peak_RSS_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)))
def main():
    ap=argparse.ArgumentParser();ap.add_argument('group',choices=['A','B','all']);args=ap.parse_args()
    for group in (['A','B'] if args.group=='all' else [args.group]):predict_group(group)
    seals={group:OUT/'raw/predictions'/group/'ONLINE_COMPLETE.json' for group in ['A','B']}
    if all(f.exists() for f in seals.values()):
        artifacts=[]
        for file in seals.values():artifacts.extend(json.loads(file.read_text())['artifacts'])
        write(OUT/'raw/predictions/ALL_PREDICTIONS_SEALED.json',dict(status='ALL_512_SEALED_BEFORE_FIRST_NEW_GT_SCORE',groups={k:str(f.resolve()) for k,f in seals.items()},
              artifacts=artifacts,utc=time.time(),truth_read=False,new_training_updates=0,total_parents=512))
if __name__=='__main__':main()
