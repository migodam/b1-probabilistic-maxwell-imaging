import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',PYTORCH_MPS_HIGH_WATERMARK_RATIO='0.12',PYTORCH_MPS_LOW_WATERMARK_RATIO='0.1')
import sys,json,time,resource,traceback
from pathlib import Path
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[2];OLD=ROOT.parent/'b1-p3l-r4-conditional-transport';sys.path.insert(0,str(OLD/'src'))
import numpy as np,torch
from p3r4.risk import LocalRisk,roi_components,full_material
from p3r2.residual_risk import read_anchor,online_features,fit_normalization,normalize
from rrr_controls import covariance_fit,prediction,pca_fit,RANKS,RIDGES
ARMS=['raw_y','decoder_only','physical_residual','same_spectrum_random']
def guard():
 if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss>3*1024**3:raise MemoryError('3GB RSS gate')
def feature_fit(arm,rt,rv,seed):
 if arm=='same_spectrum_random':
  rng=np.random.default_rng(seed);p=rt['physical_residual'];pm,ps,pv,s=pca_fit(p);X=rt['raw_y'];V=rv['raw_y'];mean=X.mean(0);std=X.std(0);std[std<1e-8]=1;R=rng.normal(size=(512,128))/np.sqrt(512);Z=(X-mean)/std@R;zm=Z.mean(0);Z-=zm;ev,eV=np.linalg.eigh(Z.T@Z/len(Z));W=(eV/np.sqrt(np.maximum(ev,1e-10)))@eV.T;factor=W@np.diag(s);return Z@factor,(((V-mean)/std@R)-zm)@factor,dict(raw_mean=mean,raw_scale=std,random_projection=R,projection_mean=zm,spectrum_factor=factor,physical_spectrum=s)
 X=rt[arm];V=rv[arm];mean,std,P,s=pca_fit(X);return (X-mean)/std@P,(V-mean)/std@P,dict(raw_mean=mean,raw_scale=std,projection=P,feature_singular_std=s)
def main():
 torch.set_num_threads(1);start=time.monotonic();device='mps' if torch.backends.mps.is_available() else 'cpu'
 if device=='mps':torch.mps.set_per_process_memory_fraction(.10)
 g=torch.load(OUT/'geometry.pt',weights_only=False,map_location='cpu')['geometry'];tr=dict(np.load(OUT/'prepared_train.npz'));va=dict(np.load(OUT/'prepared_validation.npz'));S=OLD/'research/p3l_r4_conditional_transport';a=read_anchor(S/'inputs/strong_context_1024_256/anchor_train.npz');f=online_features(a,g);target=tr['target'];anchor=tr['anchor'];guard();results=[];basefits=[]
 def raw(d,rawphysical=None):
  return dict(raw_y=d['y'].astype(float),decoder_only=np.concatenate((d['anchor'].reshape(len(d['y']),-1),d['q0']),1).astype(float),physical_residual=np.concatenate(((rawphysical if rawphysical is not None else d['image'][:,2:6]).reshape(len(d['y']),-1),d['q0']),1).astype(float))
 rawt=raw(tr,f['image'][:,2:6]);av=read_anchor(S/'inputs/strong_context_1024_256/anchor_validation.npz');fv=online_features(av,g);rawv=raw(va,fv['image'][:,2:6]);del av,fv
 def basefit(trainids,testids,seed,name):
  dest=OUT/'low_rank_correction'/'nested_baselines'/name;dest.mkdir(parents=True,exist_ok=True)
  if (dest/'RESULT.json').exists():
   basefits.append(json.loads((dest/'RESULT.json').read_text()));return np.load(dest/'heldout.npz')['prediction']
  timer=time.monotonic();torch.manual_seed(seed);rng=np.random.default_rng(seed);norm=fit_normalization({k:v[trainids] for k,v in f.items()});nf=normalize(f,norm);image=np.concatenate((nf['image'],np.broadcast_to(nf['base'][:,-2:,None,None],(1024,2,24,24))),1).astype('float32');x=torch.tensor(image,device=device);b=torch.tensor(anchor,device=device);t=torch.tensor(target,device=device);model=LocalRisk('local_cnn',g).to(device);opt=torch.optim.AdamW(model.parameters(),lr=5e-4,weight_decay=1e-3)
  with (dest/'losses.jsonl').open('w') as log:
   for u in range(1,4001):
    ids=torch.as_tensor(rng.choice(trainids,size=64,replace=True),device=device);model.train();opt.zero_grad(set_to_none=True);loss=(model(x[ids],b[ids])-t[ids]).square().mean();loss.backward();torch.nn.utils.clip_grad_norm_(model.parameters(),5);opt.step()
    if u%100==0:log.write(json.dumps(dict(update=u,train_component_MSE=float(loss.item()),wall_seconds=time.monotonic()-timer))+'\n');log.flush();guard()
  model.eval();ps=[]
  with torch.inference_mode():
   for begin in range(0,len(testids),64):ids=testids[begin:begin+64];ps.append(model(x[ids],b[ids]).cpu().numpy())
  ps=np.concatenate(ps);np.savez_compressed(dest/'heldout.npz',parent=tr['parent'][testids],prediction=ps);torch.save(dict(state_dict={k:v.cpu() for k,v in model.state_dict().items()},normalization=norm,geometry=g,training_parents=tr['parent'][trainids],seed=seed,updates=4000),dest/'BEST.pt');r=dict(status='COMPLETE',name=name,seed=seed,train_size=len(trainids),heldout_size=len(testids),train_heldout_parent_overlap=0,parameters=sum(p.numel() for p in model.parameters()),updates=4000,normalizer_fit='ownTRAINonly',device=device,GPU_seconds=time.monotonic()-timer if device=='mps' else 0,wall_seconds=time.monotonic()-timer,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,physics_RHS=0);(dest/'RESULT.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True);basefits.append(r);return ps
 for seed,oldseed in zip([2026102601,2026102602,2026102603],[2026102501,2026102502,2026102503]):
  rng=np.random.default_rng(seed);outer=np.empty(1024,int);outer[rng.permutation(1024)]=np.arange(1024)%2;outerpred=np.zeros_like(target);nestedpred={arm:{(ridge,k):np.zeros_like(target) for ridge in RIDGES for k in RANKS} for arm in ARMS}
  for fold in [0,1]:
   trainids=np.flatnonzero(outer!=fold);testids=np.flatnonzero(outer==fold);op=basefit(trainids,testids,seed+fold*10,f's{seed}_outer{fold}');outerpred[testids]=op;inner=np.empty(512,int);inner[rng.permutation(512)]=np.arange(512)%2;innerpred=np.zeros_like(target[trainids])
   for infold in [0,1]:
    keep=trainids[inner!=infold];ev=trainids[inner==infold];ip=basefit(keep,ev,seed+fold*10+infold+1,f's{seed}_outer{fold}_inner{infold}');innerpred[inner==infold]=ip
   U=(target[trainids]-innerpred).reshape(512,-1).astype(float)/24;rt={k:v[trainids] for k,v in rawt.items()};rv={k:v[testids] for k,v in rawt.items()}
   for arm in ARMS:
    X,V,feat=feature_fit(arm,rt,rv,seed+fold*100)
    for ridge in RIDGES:
     fit=covariance_fit(X,U,ridge)
     for k in RANKS:nestedpred[arm][(ridge,k)][testids]=op+prediction(fit,V,k).reshape(op.shape)*24
  U=(target-outerpred).reshape(1024,-1).astype(float)/24;finalbase=va[f'strong_s{oldseed}'];base_oof=np.sum((target-outerpred).reshape(1024,-1)**2,1)/576
  for arm in ARMS:
   dest=OUT/'low_rank_correction'/f'nested_{arm}_s{seed}';dest.mkdir(parents=True,exist_ok=True);X,V,feat=feature_fit(arm,rawt,rawv,seed);best=None;candidates=[]
   for ridge in RIDGES:
    fit=covariance_fit(X,U,ridge)
    for k in RANKS:
     ps=finalbase+prediction(fit,V,k).reshape(finalbase.shape)*24;full=full_material(ps,g);vr=np.mean(abs(full-va['truth'])**2,1);npred=nestedpred[arm][(ridge,k)];oof=np.sum((target-npred).reshape(1024,-1)**2,1)/576;row=dict(rank=k,ridge=ridge,validation_full_MSE=float(vr.mean()),nested_baseline_OOF_MSE=float(base_oof.mean()),nested_corrected_OOF_MSE=float(oof.mean()),train_singular_values=fit[4].tolist());candidates.append(row)
     if best is None or row['validation_full_MSE']<best['validation_full_MSE']:best=row;bestfull=full;bestfit=fit;bestOOF=oof;bestOOFpred=npred
   xm,um,W,P,s,Q=bestfit;k=best['rank'];np.savez_compressed(dest/'BEST.npz',input_mean=xm,output_mean=um,input_factor=W@Q[:k].T,output_factor=(P[:,:k]*s[:k]).T,rank=k,ridge=best['ridge'],arm=arm,seed=seed,baseline_old_seed=oldseed,**feat);np.savez_compressed(dest/'validation.npz',parent=va['parent'],prediction=bestfull,risk=np.mean(abs(bestfull-va['truth'])**2,1));np.savez_compressed(dest/'nested_OOF.npz',parent=tr['parent'],prediction=full_material(bestOOFpred,g),risk=bestOOF,baseline_risk=base_oof,outer_fold=outer)
   r=dict(status='COMPLETE',arm=arm,seed=seed,baseline_old_seed=oldseed,selected=best,candidates=candidates,nested_train_fold_sizes=[512,256],normalizers='baselineownTRAINonly; RRRouterTRAINonlyPCAmeanstd; fixedoldyonlygeometry',finalresidual_regime='halfTRAINbaselineOOFresidual1024; fullfrozenTRAINbaselineonVAL, declareddata-regime mismatch',baseline_selection='4000fixedupdates; no heldoutGTselection',feature_seed_semantics='deterministicRRRpairedbaselinefamily; randomcontrol+foldassignmentsseeded',physics_RHS=0)
   (dest/'RESULT.json').write_text(json.dumps(r,indent=2)+'\n');results.append(r);print(json.dumps(dict(arm=arm,seed=seed,selected={k:v for k,v in best.items() if k!='train_singular_values'})),flush=True)
 (OUT/'NESTED_RRR_SUMMARY.json').write_text(json.dumps(dict(status='COMPLETE',records=[dict(arm=r['arm'],seed=r['seed'],selected=r['selected']) for r in results],baseline_fits=basefits,GPU_seconds=sum(r['GPU_seconds'] for r in basefits),wall_seconds=time.monotonic()-start,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,physics_RHS=0),indent=2)+'\n')
if __name__=='__main__':
 status='FAILED';error=None;start=time.monotonic()
 try:main();status='COMPLETE'
 except BaseException as e:error=repr(e);traceback.print_exc();raise
 finally:(OUT/'NESTED_RRR_JOB_CLOSED.json').write_text(json.dumps(dict(status=status,error=error,wall_seconds=time.monotonic()-start,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,physics_RHS=0),indent=2)+'\n')
