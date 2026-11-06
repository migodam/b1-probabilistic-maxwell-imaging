"""NN repair: add legal paper-acquisition training, never Austria labels."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
import argparse,json,time,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
from research.research_max_struct_r8.experiments.generate_warm_data import parameters,material,FAMILIES
from research.research_max_struct_r8.src.paper_physics import acquisition,make_model
from research.research_max_struct_r8.src.warm_features import FeatureOperator

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--split',choices=['train','validation'],required=True);args=ap.parse_args()
    begin,cpu=time.perf_counter(),time.process_time();seed=2026110620+(args.split=='validation');rng=np.random.default_rng(seed);count=256 if args.split=='train' else 64
    out=ROOT/'data_training/paper_domain';out.mkdir(parents=True,exist_ok=True);path=out/(args.split+'.npz')
    if path.exists():raise FileExistsError(path)
    cache={};X=[];labels=[];records=[];receipts=[]
    for i in range(count):
        f=[400e6,750e6][i%2]
        if f not in cache:
            p=acquisition(32,f,dense=True);cache[f]=(p,make_model(p),FeatureOperator(p))
        p,model,feat=cache[f];par=parameters(rng,FAMILIES[(i//2)%8]);amp=float(rng.choice([1,2.5,3,4,6]));m=material(32,par,amp)-.1
        t,c=time.perf_counter(),time.process_time();before=dict(model.cost);clean=model.forward(m)
        noise=float(rng.choice([.03,.1]));sigma=noise*np.linalg.norm(clean)/np.sqrt(clean.size);Y=clean+sigma/np.sqrt(2)*(rng.normal(size=clean.shape)+1j*rng.normal(size=clean.shape))
        z,_,fc=feat(Y);z[-2:]*=.48
        z=np.concatenate((z,np.full((1,32,32),2*p['k']/(2*np.pi),dtype=np.float32),np.ones((1,32,32),np.float32)))
        X.append(z);labels.append(m.reshape(32,32));records.append(dict(object_id=f'r8_domain_{args.split}_{i:04d}',family=par['family'],parameters=par,amplitude=amp,frequency_Hz=f,noise_relative=noise))
        receipts.append(dict(object_id=records[-1]['object_id'],cpu_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-t,physics_cost={k:model.cost.get(k,0)-before.get(k,0) for k in model.cost},features=fc))
    np.savez_compressed(path,X=np.array(X),m=np.array(labels,dtype=np.float32))
    (out/(args.split+'_manifest.json')).write_text(json.dumps(dict(seed=seed,records=records),indent=2))
    rec=dict(status='GENERATED_NON_AUSTRIA_SUPERVISED_DOMAIN_DATA',split=args.split,seed=seed,objects=count,source_types='20 plane waves,40 receivers radius3.75m,2m ROI,32² equalarea discrete physics',truth_access='TRAIN_OR_VALIDATION_ONLY',test_objects=0,records=receipts,cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-begin)
    (ROOT/'audit'/('generation_domain_'+args.split+'.json')).write_text(json.dumps(rec,indent=2));print(json.dumps({k:rec[k] for k in ('status','split','objects','cpu_seconds','wall_seconds')}),flush=True)

if __name__=='__main__':main()
