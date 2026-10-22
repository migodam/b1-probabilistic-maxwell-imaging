"""Dedicated joint-supervision training, never online truth access."""
from pathlib import Path
import sys,argparse,json
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3l.core import OUT,Book,write
from p3l.neural_training import train_density

def main():
    a=argparse.ArgumentParser();a.add_argument('model',choices=['mlp_diag','mlp_full','mlp_mdn','mlp_flow','cnn_full']);a.add_argument('--seed',type=int,default=2026102270);args=a.parse_args()
    if args.seed not in [2026102270,2026102271,2026102272]:raise RuntimeError('undeclared training seed')
    meta=json.loads((OUT/'D/training/PACKING_METADATA.json').read_text())
    if meta['groups']!={'train':512,'validation':128}:raise RuntimeError('frozen supervised dataset missing')
    b=Book('D',f'train_{args.model}_s{args.seed}',max_wall=1260);status='COMPLETE';error=None
    try:
        with np.load(OUT/'D/training/TRAIN.npz',allow_pickle=False) as z:train={k:z[k].copy() for k in ['scalar','image','labels']}
        with np.load(OUT/'D/training/VALIDATION.npz',allow_pickle=False) as z:val={k:z[k].copy() for k in ['scalar','image','labels']}
        result=train_density(args.model,train['scalar'],train['image'],train['labels'],val['scalar'],val['image'],val['labels'],output_dir=OUT/f'D/models/{args.model}/seed{args.seed}',label_coordinates=meta['label_coordinates'],seed=args.seed,device='cpu',supervision='JOINT_SIMULATION',qualified_distillation=False,max_updates=3000,max_wall_seconds=1200)
        write(b.jobdir/'TRAINING_SUMMARY.json',dict(model=args.model,seed=args.seed,supervision='joint simulation; no posterior-reference distillation',new_physics_RHS=0,GPU_seconds=0,metadata=result['metadata'],checkpoint_path=result['checkpoint_path'],logs=result['logs']))
        print(json.dumps(dict(model=args.model,seed=args.seed,status='COMPLETE')),flush=True)
    except BaseException as e:status='ERROR';error=repr(e);raise
    finally:b.close(status,error)

if __name__=='__main__':main()
