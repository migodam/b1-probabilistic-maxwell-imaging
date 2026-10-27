"""Export first512 historical TRAIN observations, full contexts and separate labels.

No new observations, Maxwell solves, test labels or CDEV labels are generated.
The only supervised source is Original R4 joint/train32768.npz. Online archives
are constructed by an explicit allowlist and never include chi_GT/q_GT/x_GT.
"""
from __future__ import annotations
import os
for _name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS',
              'VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_name]='1'
import argparse
import json
from pathlib import Path
import sys
import time
import numpy as np

TRACK=Path(__file__).resolve().parent
OUT=TRACK.parent
REPO=OUT.parents[1]
R4=REPO.parent/'b1-p3l-r4-conditional-transport'
SOURCE=R4/'research/p3l_r4_conditional_transport'
sys.path.insert(0,str(OUT))
sys.path.insert(0,str(R4/'src'))
from track_c.frozen_heads import ONLINE_KEYS,SEEDS,packed_measurement,load_action


def export(destination,n=512,teachers=False):
    if n!=512:raise ValueError('revision1 matched subset is frozen first512 OriginalTRAIN')
    destination=Path(destination).resolve()
    if not destination.is_relative_to(TRACK.resolve()) or destination.exists():
        raise FileExistsError('new Track C export destination required')
    destination.mkdir(parents=True)
    start,cpu=time.perf_counter(),time.process_time()
    anchor_source=SOURCE/'inputs/strong_context_1024_256/anchor_train.npz'
    joint_source=SOURCE/'inputs/joint/train32768.npz'
    with np.load(anchor_source,allow_pickle=False) as z:
        parents=z['parent'][:n].astype(str)
        online={key:z[key][:n].copy() for key in ONLINE_KEYS}
        prediction=z['prediction'][:n].copy()
    if len(parents)!=512 or len(set(parents))!=512:raise ValueError('512 unique original TRAIN parents')
    with np.load(joint_source,allow_pickle=False) as z:
        ids={str(parent):i for i,parent in enumerate(z['parent'])}
        indices=np.array([ids[parent] for parent in parents])
        y=z['y'][indices].copy();truth=z['material'][indices].copy()
        target_id=str(z['target_id'])
    expected_y=np.stack([packed_measurement(obs) for obs in online['observation']])
    if not np.array_equal(y,expected_y):raise ValueError('Original TRAIN observation identity mismatch')
    cfg=json.loads((SOURCE/'configs/TARGET.json').read_text())
    (destination/'online').mkdir();(destination/'full_contexts').mkdir();(destination/'offline_training').mkdir()
    for i,parent in enumerate(parents):
        np.savez_compressed(destination/'online'/f'{parent}.npz',observation=online['observation'][i],
                            sigma=np.array(cfg['noise']['sigma']),parent=np.array(parent),target_id=np.array(target_id))
        context={key:online[key][i] for key in ONLINE_KEYS}
        np.savez_compressed(destination/'full_contexts'/f'{parent}.npz',**context,
                            prediction=prediction[i],parent=np.array(parent))
    np.savez_compressed(destination/'offline_training/LABELS_TRAIN512.npz',parent=parents,material=truth)
    manifest=dict(group='HISTORICAL_TRAIN512',n=512,parents=[dict(parent=str(parent)) for parent in parents],
                  order='first512 of Original strong_context_1024_256/anchor_train',new_physics_RHS=0)
    (destination/'PARENTS.json').write_text(json.dumps(manifest,indent=2)+'\n')
    teacher_records=[]
    if teachers:
        from p3r4.risk import feature_image,roi_components
        for seed in SEEDS:
            tick,cpu0=time.perf_counter(),time.process_time()
            action=load_action('distilled_w128_n32768',seed)
            payload=action.p;normalization=payload.get('feature_payload',payload)['feature_normalization']
            moments=action.moments(y,parents)
            normalized=((moments-payload['moment_mean'])/payload['moment_scale']).astype('float32')
            image=feature_image(y,action.geometry,online,normalization)
            anchor=roi_components(online['chi0'],action.geometry['roi'])
            path=destination/f'FULL_TEACHER_FEATURES_s{seed}.npz'
            np.savez_compressed(path,parent=parents,y=y,image=image,anchor=anchor,moments=normalized)
            teacher_records.append(dict(seed=seed,path=str(path),wall_seconds=time.perf_counter()-tick,
                                        CPU_seconds=time.process_time()-cpu0,new_RHS=0,labels_read_by_teacher=False))
            del action
    receipt=dict(status='EXPORTED_ORIGINAL_TRAIN512',parents=512,original_R5_head_training_n=1024,
        matched_new_full_head_required=True,subset='identical first512 parents for newFull and newLow',
        anchor_source=str(anchor_source),supervised_source=str(joint_source),new_physics_RHS=0,
        online_keys=list(ONLINE_KEYS),online_labels_present=False,labels_destination='offline_training/LABELS_TRAIN512.npz',
        teacher_records=teacher_records,wall_seconds=time.perf_counter()-start,CPU_seconds=time.process_time()-cpu,
        GT_scope='historical TRAIN only; no CDEV/test/D labels',original_full_contexts_reused=True)
    (destination/'EXPORT_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n')
    return receipt


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=TRACK/'matched_inputs_TRAIN512')
    parser.add_argument('--n',type=int,default=512)
    parser.add_argument('--teachers',action='store_true')
    args=parser.parse_args();print(json.dumps(export(args.output,args.n,args.teachers)),flush=True)


if __name__=='__main__':main()
