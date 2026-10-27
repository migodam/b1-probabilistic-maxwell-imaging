"""Paired new Full512/Low512 FusionRisk heads; frozen compact128 provider.

Training updates see historical TRAIN512 labels only. CDEV64 is explicitly the
checkpoint-selection set; D and historical test labels are never constructed.
Root runs this only after the complete low-context guard and owns the decision
to spend the training job wall budget. No Target or solver is constructed.
"""
from __future__ import annotations
import os
for _name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS',
              'VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_name]='1'
import argparse,copy,json,sys,time,signal
from pathlib import Path
import numpy as np
import torch
TRACK=Path(__file__).resolve().parent;OUT=TRACK.parent;sys.path.insert(0,str(OUT))
from track_c.frozen_heads import load_action,packed_measurement,predict_contexts,ONLINE_KEYS,SEEDS
from shared.physics import write
from p3r4.risk import feature_image,roi_components,full_material
from p3r4.risk_fusion import FusionRisk

FORMAT='R6_MATCHED_TRAIN512_FULL_LOW_FUSION6418_V1'


class TrainingAccount:
    """Separate zero-physics head-training ledger; permits root CPU physics.

    Whole MPS job wall is a conservative GPU reservation. No unobserved kernel
    time estimate is invented. All fits together have the root-assigned1800s cap.
    """
    def __init__(self,job,device,max_wall):
        self.job,self.device=job,device
        self.path=TRACK/'MATCHED_TRAINING_LEDGER.jsonl'
        previous=[json.loads(line) for line in self.path.read_text().splitlines()] if self.path.exists() else []
        live={row['job'] for row in previous if row['event']=='BEGIN'}-{row['job'] for row in previous if row['event']=='END'}
        if live:raise RuntimeError('unreconciled prior matched training job')
        used=sum(row.get('GPU_reserved_whole_job_seconds',0) for row in previous if row['event']=='END')
        remaining=1800-used if device=='mps' else float(max_wall)
        if remaining<=0:raise RuntimeError('assigned1800s MPS training cap reached')
        self.started,self.cpu=time.perf_counter(),time.process_time()
        self.limit=min(float(max_wall),remaining)
        self._append(dict(event='BEGIN',job=job,device=device,new_RHS=0,max_wall_seconds=self.limit))
        signal.signal(signal.SIGALRM,lambda *_:(_ for _ in ()).throw(RuntimeError('matched training assigned wall cap')))
        signal.setitimer(signal.ITIMER_REAL,self.limit)
    def _append(self,row):
        with self.path.open('a') as stream:stream.write(json.dumps(dict(utc=time.time(),**row))+'\n')
    def close(self,status,error):
        if self.device=='mps':torch.mps.synchronize()
        signal.setitimer(signal.ITIMER_REAL,0)
        wall=time.perf_counter()-self.started
        self._append(dict(event='END',job=self.job,status=status,error=error,device=self.device,
            wall_seconds=wall,CPU_seconds=time.process_time()-self.cpu,new_RHS=0,
            GPU_reserved_whole_job_seconds=wall if self.device=='mps' else 0))


def context_batch(paths,parents):
    rows=[]
    for path,parent in zip(paths,parents):
        with np.load(path,allow_pickle=False) as z:
            if str(z['parent'])!=parent:raise ValueError('matched parent/context identity')
            rows.append({key:z[key].copy() for key in ONLINE_KEYS})
    return {key:np.stack([row[key] for row in rows]) for key in ONLINE_KEYS}


def original_features(action,online,y,moments=None):
    payload=action.p;geometry=action.geometry
    normalization=payload.get('feature_payload',payload)['feature_normalization']
    image=feature_image(y,geometry,online,normalization)
    anchor=roi_components(online['chi0'],geometry['roi'])
    if moments is None:
        raw=action.moments(y)
        moments=((raw-payload['moment_mean'])/payload['moment_scale']).astype('float32')
    return image,np.asarray(moments,'float32'),anchor


def require_complete(inputs,low_root,low_config):
    manifest=json.loads((inputs/'PARENTS.json').read_text())
    parents=[str(row['parent']) for row in manifest['parents']]
    if manifest['group']!='HISTORICAL_TRAIN512' or len(parents)!=512:raise ValueError('Original TRAIN512 subset required')
    marker=json.loads((low_root/'ANCHORS_COMPLETE.json').read_text())
    if marker.get('status')!='COMPLETE' or marker.get('n')!=512 or marker.get('parents')!=parents or marker.get('config')!=low_config:
        raise RuntimeError('HOLD: complete matching512 low contexts required before any head fit')
    for parent in parents:
        if not (low_root/parent/low_config/'CONTEXT.npz').exists():raise RuntimeError('missing low context')
    return parents


def train(inputs,low_root,output,low_config='native48_probit32',max_wall=1800,systems=('full','low'),device='cpu'):
    torch.set_num_threads(1)
    parents=[row['parent'] for row in json.loads((inputs/'PARENTS.json').read_text())['parents']]
    if len(parents)!=512:raise ValueError('TRAIN512 required')
    if 'low' in systems:require_complete(inputs,low_root,low_config)
    native_root=TRACK/'native96_only_contexts/TRAIN512'
    if 'native' in systems:require_complete(inputs,native_root,'native96_only')
    if device not in ('cpu','mps') or device=='mps' and not torch.backends.mps.is_available():raise ValueError('available CPU/MPS backend required')
    if device=='mps':torch.mps.set_per_process_memory_fraction(.10)
    devparents=[str(row['parent']) for row in json.loads((OUT/'configs/parents_CDEV.json').read_text())['parents']]
    if len(devparents)!=64 or set(devparents)&set(parents):raise ValueError('independent CDEV selection set')
    with np.load(inputs/'offline_training/LABELS_TRAIN512.npz',allow_pickle=False) as z:
        if z['parent'].astype(str).tolist()!=parents:raise ValueError('TRAIN label order')
        train_truth=z['material'].copy()
    # CDEV labels are allowed solely for the declared checkpoint selection.
    dev_truth=[]
    for parent in devparents:
        with np.load(OUT/'raw/offline_truth'/f'{parent}.npz',allow_pickle=False) as z:dev_truth.append(z['chi'].copy())
    dev_truth=np.array(dev_truth)
    online_train={}
    if 'low' in systems:online_train['low']=context_batch([low_root/p/low_config/'CONTEXT.npz' for p in parents],parents)
    if 'native' in systems:online_train['native']=context_batch([native_root/p/'native96_only/CONTEXT.npz' for p in parents],parents)
    dev={}
    for arm in systems:
        root=TRACK/'native96_only_contexts/CDEV' if arm=='native' else TRACK/'contexts/CDEV'
        config='native96_only' if arm=='native' else ('native96_probit64' if arm=='full' else low_config)
        dev[arm]=context_batch([root/p/config/'CONTEXT.npz' for p in devparents],devparents)
    with np.load(inputs/f'FULL_TEACHER_FEATURES_s{SEEDS[0]}.npz',allow_pickle=False) as z:y=z['y'].copy()
    dev_y=np.stack([packed_measurement(obs) for obs in next(iter(dev.values()))['observation']])
    output.mkdir(parents=True,exist_ok=True)
    preflight=dict(status='READY_TRAIN512_LOW_CONTEXTS_COMPLETE',train_parents=512,CDEV_selection_parents=64,
        low_config=low_config,arms=list(systems),seeds=list(SEEDS),updates=4000,batch_size=64,device=device,
        optimizer=dict(name='AdamW',lr=5e-4,weight_decay=1e-3,gradient_clip_norm=5.),
        capacity=6418,new_provider_updates=0,new_physics_RHS=0,D_GT_read=False,
        limitation='CDEV selects checkpoints and final candidate; independent locked validation still required',
        original_reference_training_n=1024,new_both_arms_training_n=512)
    protocol_path=output/('TRAINING_PROTOCOL_'+'_'.join(systems)+'.json')
    if protocol_path.exists():raise FileExistsError('training invocation protocol already exists')
    write(protocol_path,preflight)
    records=[]
    for seed in SEEDS:
        action=load_action('distilled_w128_n32768',seed)
        geometry=action.geometry
        with np.load(inputs/f'FULL_TEACHER_FEATURES_s{seed}.npz',allow_pickle=False) as z:
            if z['parent'].astype(str).tolist()!=parents or not np.array_equal(z['y'],y):raise ValueError('same TRAIN parents/observations for Full/Low')
            full_features=(z['image'].copy(),z['moments'].copy(),z['anchor'].copy())
        features_by_arm={'full':full_features}
        for arm,online in online_train.items():
            if not np.array_equal(np.stack([packed_measurement(obs) for obs in online['observation']]),y):raise ValueError('same TRAIN observations across arms')
            features_by_arm[arm]=original_features(action,online,y,full_features[1])
        devmoments=action.moments(dev_y)
        devmoments=((devmoments-action.p['moment_mean'])/action.p['moment_scale']).astype('float32')
        dev_features={arm:original_features(action,online,dev_y,devmoments) for arm,online in dev.items()}
        train_target=torch.tensor(roi_components(train_truth,geometry['roi']),device=device)
        dev_target=roi_components(dev_truth,geometry['roi'])
        for arm in systems:
            features=features_by_arm[arm]
            folder=output/f'{arm}_s{seed}'
            if (folder/'RESULT.json').exists():
                records.append(json.loads((folder/'RESULT.json').read_text()));continue
            if folder.exists():raise RuntimeError('partial fit must be audited; no silent retry')
            folder.mkdir()
            book=TrainingAccount(f'matched512_{arm}_{low_config}_s{seed}',device,max_wall)
            start,cpu=time.perf_counter(),time.process_time();status,error='COMPLETE',None
            try:
                torch.manual_seed(seed);rng=np.random.default_rng(seed)
                model=FusionRisk(geometry).to(device)
                count=sum(p.numel() for p in model.parameters())
                if count!=6418:raise ValueError('Original FusionRisk capacity6418 required')
                optimizer=torch.optim.AdamW(model.parameters(),lr=5e-4,weight_decay=1e-3)
                tensor=tuple(torch.tensor(a,dtype=torch.float32,device=device) for a in features)
                val=tuple(torch.tensor(a,dtype=torch.float32,device=device) for a in dev_features[arm])
                def evaluate():
                    model.eval()
                    with torch.inference_mode():prediction=model(*val).cpu().numpy()
                    return float(np.mean((prediction-dev_target)**2)),prediction
                best,prediction=evaluate();state=copy.deepcopy(model.state_dict());selected=0
                with (folder/'LOSSES.jsonl').open('x') as log:
                    for update in range(1,4001):
                        ids=torch.as_tensor(rng.integers(512,size=64),device=device);model.train();optimizer.zero_grad(set_to_none=True)
                        loss=(model(*(a[ids] for a in tensor))-train_target[ids]).square().mean()
                        if not torch.isfinite(loss):raise FloatingPointError('nonfinite matched TRAIN loss')
                        loss.backward();torch.nn.utils.clip_grad_norm_(model.parameters(),5.);optimizer.step()
                        if update%100==0:
                            score,current=evaluate()
                            log.write(json.dumps(dict(update=update,TRAIN_component_MSE=float(loss.item()),CDEV_component_MSE=score,
                                                       wall_seconds=time.perf_counter()-start))+'\n');log.flush()
                            if score<best:best,prediction,selected,state=score,current,update,copy.deepcopy(model.state_dict())
                if device=='mps':torch.mps.synchronize()
                payload=dict(format=FORMAT,arm=arm,seed=seed,low_config=low_config,state_dict={k:v.cpu() for k,v in state.items()},
                    geometry=geometry,train_parents=parents,CDEV_selection_parents=devparents,selected_update=selected,
                    updates=4000,capacity=6418,provider_method='distilled_w128_n32768',provider_frozen=True,
                    new_provider_updates=0,train_n=512,optimizer=preflight['optimizer'])
                torch.save(payload,folder/'BEST.pt')
                material=full_material(prediction,geometry)
                # These are checkpoint-selection outputs, not sealed independent
                # test predictions. Separate inference below seals frozen outputs.
                np.savez_compressed(folder/'CDEV_CHECKPOINT_SELECTION.npz',parent=np.array(devparents),prediction=material,
                                    risk=np.mean(abs(material-dev_truth)**2,1))
                record=dict(status='COMPLETE',arm=arm,seed=seed,capacity=count,train_n=512,updates=4000,
                    selected_update=selected,CDEV_complex_MSE=float(np.mean(abs(material-dev_truth)**2)),
                    training_wall_seconds=time.perf_counter()-start,training_CPU_seconds=time.process_time()-cpu,
                    provider_updates=0,new_physics_RHS=0,device=device,CDEV_selection_only=True,D_GT_read=False,
                    central_ledger_job=book.job,checkpoint=str(folder/'BEST.pt'))
                write(folder/'RESULT.json',record);records.append(record)
                print(json.dumps(record),flush=True)
            except BaseException as exc:
                status,error='FAILED',repr(exc);write(folder/'FAILURE.json',dict(error=error,new_physics_RHS=0));raise
            finally:book.close(status,error)
        del action
    write(output/('HEADS_COMPLETE_'+'_'.join(systems)+'.json'),dict(status='COMPLETE',records=records,
        total_head_fits=len(records),provider_training_updates=0,new_physics_RHS=0,CDEV_selection_only=True,
        total_recorded_training_wall_seconds=sum(r['training_wall_seconds'] for r in records),
        total_recorded_training_CPU_seconds=sum(r['training_CPU_seconds'] for r in records)))
    return records


class MatchedAction:
    """Original frozen compact128 provider plus one selected new6418 head."""
    def __init__(self,checkpoint):
        payload=torch.load(checkpoint,map_location='cpu',weights_only=False)
        if payload.get('format')!=FORMAT or payload['capacity']!=6418 or payload['train_n']!=512:
            raise ValueError('matched TRAIN512 head contract')
        self.payload=payload
        self.original=load_action('distilled_w128_n32768',payload['seed'])
        self.original.action.head.load_state_dict(payload['state_dict'])
        self.original.action.head.eval();self.geometry=self.original.geometry

    def predict(self,y,online,parents=None,batch_size=64):
        return self.original.predict(y,online,parents=parents,batch_size=batch_size)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,default=TRACK/'matched_inputs_TRAIN512')
    parser.add_argument('--low-context-root',type=Path,default=TRACK/'matched_anchors_TRAIN512')
    parser.add_argument('--low-config',default='native48_probit32')
    parser.add_argument('--output',type=Path,default=TRACK/'matched_heads_TRAIN512')
    parser.add_argument('--max-wall-per-fit',type=float,default=1800)
    parser.add_argument('--systems',nargs='+',choices=['full','low','native'],default=['full','low'])
    parser.add_argument('--device',choices=['cpu','mps'],default='mps' if torch.backends.mps.is_available() else 'cpu')
    args=parser.parse_args()
    if not args.output.resolve().is_relative_to(TRACK.resolve()):raise ValueError('Track C output only')
    train(args.input,args.low_context_root,args.output,args.low_config,args.max_wall_per_fit,args.systems,args.device)


if __name__=='__main__':main()
