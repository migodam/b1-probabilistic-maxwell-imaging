"""Train only the two frozen-architecture channel-source controls, no fresh GT."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'): os.environ[key]='1'
from pathlib import Path
import argparse,json,sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3r2.core import OUT
from p3r2.residual_risk import read_anchor,subset,fit_geometry
from p3r2.risk_training import write_json
from p3r2_mechanism.masked_risk import FAMILIES,mask_contract,train_masked_residual,load_mechanism_predictor


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--train',type=Path,default=OUT/'A/anchor_data/train.npz');p.add_argument('--validation',type=Path,default=OUT/'A/anchor_data/validation.npz')
    p.add_argument('--output',type=Path,default=OUT/'A/mechanism_controls');p.add_argument('--target',type=Path,default=OUT/'configs/TARGET.json')
    p.add_argument('--geometry',type=Path,default=ROOT/'research/p1_material_target_20261018/inputs/geometry.json');p.add_argument('--resume',action='store_true')
    a=p.parse_args()
    for path in (a.train,a.validation):
        if 'locked' in str(path).lower() or 'offline_only' in path.parts: raise ValueError('fresh/locked labels prohibited in mechanism training')
    tr=read_anchor(a.train,True);va=read_anchor(a.validation,True)
    if len(tr['parent'])<512 or len(va['parent'])!=128: raise ValueError('original512train+128validation required')
    tr=subset(tr,512)
    if set(tr['parent'])&set(va['parent']): raise ValueError('train/validation parent overlap')
    cfg=json.loads(a.target.read_text());points=json.loads(a.geometry.read_text())['points'];geometry=fit_geometry(tr,cfg,points,16,12,True)
    a.output.mkdir(parents=True,exist_ok=True)
    contract=dict(version=1,objective='SUPERVISED_RISK',train=str(a.train.resolve()),validation=str(a.validation.resolve()),
        train_count=512,validation_count=128,families=FAMILIES,seeds=[11,22,33],parameter_count=5842,
        epochs=100,batch_size=128,patience=15,max_seconds=600,learning_rate=.0005,weight_decay=.001,
        architecture='Original ResidualNet(local_geometry_cnn), unchanged',optimizer_normalization_selection='Original train_residual unchanged, including validation update0',
        masks={f:mask_contract(f) for f in FAMILIES},physics_calls=0,fresh_GT_read=False,posterior_claim=False)
    write_json(a.output/'RUN_CONTRACT.json',contract);results=[]
    for family in FAMILIES:
        for seed in (11,22,33):
            dest=a.output/'n512'/f'{family}_s{seed}'
            print(json.dumps(dict(event='start',family=family,seed=seed)),flush=True)
            if a.resume and (dest/'RESULT.json').exists():
                load_mechanism_predictor(dest/'BEST.pt');result=json.loads((dest/'RESULT.json').read_text())
                if result.get('mask')!=mask_contract(family) or result.get('seed')!=seed: raise ValueError('resume contract mismatch')
            else:
                if dest.exists(): raise FileExistsError('existing run requires explicit --resume: '+str(dest))
                result=train_masked_residual(family,seed,tr,va,geometry,dest,epochs=100,batch_size=128,patience=15,max_seconds=600,learning_rate=5e-4,weight_decay=1e-3,device='cpu')
            results.append(dict(path=str(dest.resolve()),**{k:v for k,v in result.items() if k!='history'}))
            write_json(a.output/'RESULTS.json',dict(status='COMPLETE' if len(results)==6 else 'RUNNING',results=results,objective='SUPERVISED_RISK',fresh_GT_read=False))
            print(json.dumps(dict(event='complete',family=family,seed=seed,parameter_count=result['parameter_count'],selected_update=result['selected_update'])),flush=True)
    write_json(a.output/'COMPLETE.json',dict(status='COMPLETE',runs=6,families=FAMILIES,seeds=[11,22,33],objective='SUPERVISED_RISK',fresh_GT_read=False))


if __name__=='__main__':main()
