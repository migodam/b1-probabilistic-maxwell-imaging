"""Actual supervised training; only train and validation objects are accessible."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[name]='1'
import argparse,json,time,sys,resource
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
import torch
from torch import nn
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
from research.research_max_struct_r8.src.warm_network import TinyWarmStart

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--mode',choices=['born','physics'],required=True);ap.add_argument('--epochs',type=int,default=120);ap.add_argument('--device',default='cpu');args=ap.parse_args()
    begin,cpu=time.perf_counter(),time.process_time();torch.set_num_threads(2);torch.manual_seed(2026110604);rng=np.random.default_rng(2026110604)
    out=ROOT/'RESULTS/training'/args.mode;out.mkdir(parents=True,exist_ok=False)
    with np.load(ROOT/'data_training/train.npz') as z:X=z['X'].copy();m=z['m'].copy()
    with np.load(ROOT/'data_training/validation.npz') as z:V=z['X'].copy();v=z['m'].copy()
    if args.mode=='born':X[:,2:6]=0;V[:,2:6]=0
    mean=X.mean((0,2,3));std=np.maximum(X.std((0,2,3)),.1)
    X=(X-mean[None,:,None,None])/std[None,:,None,None];V=(V-mean[None,:,None,None])/std[None,:,None,None]
    x=torch.tensor(X,device=args.device);y=torch.tensor(m,device=args.device);vx=torch.tensor(V,device=args.device);vy=torch.tensor(v,device=args.device)
    arch=dict(channels=X.shape[1],width=16);net=TinyWarmStart(**arch).to(args.device)
    opt=torch.optim.AdamW(net.parameters(),lr=.002,weight_decay=1e-5)
    scheduler=torch.optim.lr_scheduler.CosineAnnealingLR(opt,T_max=args.epochs,eta_min=.0002)
    scale=torch.maximum(y.amax((1,2))-.1,torch.ones(len(y),device=args.device));vs=torch.maximum(vy.amax((1,2))-.1,torch.ones(len(vy),device=args.device))
    history=[];best=float('inf');bestepoch=None;batch=64
    for epoch in range(args.epochs):
        net.train();losses=[]
        for a in range(0,len(x),batch):
            # One permutation per epoch, fixed independently of labels.
            if a==0:order=rng.permutation(len(x))
            idx=torch.as_tensor(order[a:a+batch],device=args.device)
            pred=net(x[idx]);loss=(((pred-y[idx])/scale[idx,None,None])**2).mean()
            opt.zero_grad();loss.backward();opt.step();losses.append(float(loss.detach().cpu()))
        scheduler.step();net.eval()
        with torch.no_grad():
            vp=torch.cat([net(vx[a:a+batch]) for a in range(0,len(vx),batch)])
            vl=float((((vp-vy)/vs[:,None,None])**2).mean().cpu())
        history.append(dict(epoch=epoch+1,train_normalized_MSE=float(np.mean(losses)),validation_normalized_MSE=vl,cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-begin))
        if vl<best:
            best=vl;bestepoch=epoch+1
            torch.save(dict(state_dict={k:t.detach().cpu() for k,t in net.state_dict().items()},architecture=arch,mean=mean.tolist(),std=std.tolist(),feature_mode=args.mode,seed=2026110604,selected_epoch=bestepoch,selection='Minimum independent-material validation normalized MSE; no test data',training_objects=len(x),validation_objects=len(vx)),out/'weights.pt')
        if (epoch+1)%10==0:print(json.dumps(history[-1]),flush=True)
        (out/'history.json').write_text(json.dumps(history,indent=2))
    receipt=dict(status='TRAINED_AND_VALIDATION_SELECTED',feature_mode=args.mode,parameters=sum(p.numel() for p in net.parameters()),architecture=arch,epochs=args.epochs,batch_size=batch,
        best_epoch=bestepoch,validation_normalized_MSE=best,training_objects=len(x),validation_objects=len(vx),test_objects_read=0,
        optimizer='AdamW lr.002 cosine to.0002 weight_decay1e-5',seed=2026110604,device=args.device,torch_version=torch.__version__,cpu_threads=2,
        cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-begin,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        LU=0,complex_RHS=0,truth_access='TRAIN_AND_VALIDATION_LABELS_ONLY',generation_cost_receipts=['audit/generation_train.json','audit/generation_validation.json'],finished_utc=datetime.now(timezone.utc).isoformat())
    (out/'RECEIPT.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt),flush=True)

if __name__=='__main__':main()
