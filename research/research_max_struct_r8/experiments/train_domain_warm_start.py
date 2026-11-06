"""One physically motivated NN repair: explicit acquisition features and data.

The architecture depth/width stay fixed. No primary test or Austria GT labels.
"""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
import json,time,sys,resource
from pathlib import Path
import numpy as np
import torch
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
from research.research_max_struct_r8.src.warm_network import TinyWarmStart

def main():
    begin,cpu=time.perf_counter(),time.process_time();torch.set_num_threads(1);torch.manual_seed(2026110622);rng=np.random.default_rng(2026110622)
    out=ROOT/'RESULTS/training/domain_augmented';out.mkdir(parents=True,exist_ok=False)
    arrays=[]
    for domain in ['original','paper']:
        for split in ['train','validation']:
            path=ROOT/'data_training'/((split+'.npz') if domain=='original' else ('paper_domain/'+split+'.npz'))
            with np.load(path) as z:X=z['X'].copy();m=z['m'].copy()
            if domain=='original':X=np.concatenate((X,np.full((len(X),1,24,24),.96,np.float32),np.zeros((len(X),1,24,24),np.float32)),axis=1)
            arrays.append((X,m))
    mean=(arrays[0][0].mean((0,2,3))+arrays[2][0].mean((0,2,3)))/2
    second=(np.mean(arrays[0][0]**2,(0,2,3))+np.mean(arrays[2][0]**2,(0,2,3)))/2
    std=np.sqrt(np.maximum(second-mean**2,.01));tensors=[]
    for X,m in arrays:
        x=torch.tensor((X-mean[None,:,None,None])/std[None,:,None,None]);y=torch.tensor(m)
        scale=torch.maximum(y.amax((1,2))-(.1 if m.shape[-1]==24 else 0),torch.ones(len(y)))
        tensors.append((x,y,scale))
    arch=dict(channels=12,width=16);net=TinyWarmStart(**arch)
    # Same small architecture, learned weights retained except new feature columns.
    old=torch.load(ROOT/'RESULTS/training/physics/weights.pt',map_location='cpu',weights_only=False)['state_dict'];state=net.state_dict()
    for key in state:
        if key=='enc.0.weight':state[key].zero_();state[key][:,:10]=old[key]
        else:state[key]=old[key]
    net.load_state_dict(state);opt=torch.optim.AdamW(net.parameters(),lr=.0008,weight_decay=1e-5);epochs=80;history=[];best=float('inf');bestep=None
    for epoch in range(epochs):
        net.train();losses=[]
        # Balanced domain updates, 16 batches per domain; 1024 old/1024 resampled new objects per epoch.
        for step in range(32):
            x,y,sc=tensors[0 if step%2==0 else 2];idx=rng.choice(len(x),64,replace=len(x)<64)
            pred=net(x[idx]);loss=(((pred-y[idx])/sc[idx,None,None])**2).mean();opt.zero_grad();loss.backward();opt.step();losses.append(float(loss.detach()))
        net.eval();vals=[]
        with torch.no_grad():
            for index in (1,3):
                x,y,sc=tensors[index];p=torch.cat([net(x[a:a+64]) for a in range(0,len(x),64)]);vals.append(float((((p-y)/sc[:,None,None])**2).mean()))
        score=float(np.mean(vals));history.append(dict(epoch=epoch+1,train_MSE=float(np.mean(losses)),validation_original_MSE=vals[0],validation_paper_MSE=vals[1],balanced_validation_MSE=score,cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-begin))
        if score<best:
            best=score;bestep=epoch+1
            torch.save(dict(state_dict=net.state_dict(),architecture=arch,mean=mean.tolist(),std=std.tolist(),feature_mode='physics_domain_augmented',selected_epoch=bestep,selection='Equal weight original/paper-domain independent validation MSE; no Austria/test labels',seed=2026110622,new_features=['electrical_ROI_width_in_wavelengths','plane_wave_source_flag']),out/'weights.pt')
        (out/'history.json').write_text(json.dumps(history,indent=2))
        if (epoch+1)%10==0:print(json.dumps(history[-1]),flush=True)
    rec=dict(status='TRAINED_PHYSICS_DOMAIN_REPAIR',parameters=sum(p.numel() for p in net.parameters()),architecture=arch,epochs=epochs,best_epoch=bestep,balanced_validation_MSE=best,training_objects=1792,validation_objects=256,test_objects_read=0,Austria_GT_labels_used=0,CPU_threads=1,device='CPU',cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-begin,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,LU=0,complex_RHS=0,previous_checkpoint='physics/weights.pt',prior_original_training_cost_not_free=True)
    (out/'RECEIPT.json').write_text(json.dumps(rec,indent=2));print(json.dumps(rec),flush=True)

if __name__=='__main__':main()
