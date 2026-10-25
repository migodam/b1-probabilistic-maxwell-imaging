"""Train-only reduced rank regression of predictable material residuals.
Train input covariance whitening, output-input crosscov SVD (never GT-PCA).
Five-fold OOF correction diagnostic and VAL-only ridge/rank selection.
"""
import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
import sys,json,time,resource
from pathlib import Path
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[2];OLD=ROOT.parent/'b1-p3l-r4-conditional-transport';sys.path.insert(0,str(OLD/'src'))
import numpy as np,torch
from p3r4.risk import full_material
from p3r2.fusion import passive_box
RANKS=[2,4,8,16];RIDGES=[.001,.01,.1,1.]

def covariance_fit(X,U,ridge):
 xm=X.mean(0);um=U.mean(0);x=X-xm;u=U-um;B=x.T@x/len(x);lam=ridge*np.trace(B)/len(B);val,V=np.linalg.eigh(B);W=(V/np.sqrt(np.maximum(val+lam,1e-15)))@V.T;T=(u.T@x/len(x))@W;P,s,Q=np.linalg.svd(T,full_matrices=False)
 return xm,um,W,P[:,:16],s,Q[:16]
def prediction(fit,X,k):
 xm,um,W,P,s,Q=fit;return (X-xm)@W@Q[:k].T@(P[:,:k]*s[:k]).T+um

def pca_fit(X,n=128):
 mean=X.mean(0);std=X.std(0);std=np.where(std<1e-8,1,std);x=(X-mean)/std;U,s,V=np.linalg.svd(x,full_matrices=False);V=V[:n].T;return mean,std,V,s[:n]/np.sqrt(len(X))
def main():
 torch.set_num_threads(1);start=time.monotonic();g=torch.load(OUT/'geometry.pt',weights_only=False,map_location='cpu')['geometry'];tr=dict(np.load(OUT/'prepared_train.npz'));va=dict(np.load(OUT/'prepared_validation.npz'));records=[]
 def raw(data,arm):
  if arm=='raw_y':return data['y'].astype(float)
  if arm=='decoder_only':return np.concatenate([data['anchor'].reshape(len(data['y']),-1),data['q0']],axis=1).astype(float)
  if arm=='physical_residual':return np.concatenate([data['image'][:,2:6].reshape(len(data['y']),-1),data['q0']],axis=1).astype(float)
  raise ValueError(arm)
 for seed,oldseed in zip([2026102601,2026102602,2026102603],[2026102501,2026102502,2026102503]):
  rng=np.random.default_rng(seed);fold=np.empty(len(tr['y']),int);order=rng.permutation(len(fold));fold[order]=np.arange(len(fold))%5
  Bt=tr[f'strong_s{oldseed}'];Bv=va[f'strong_s{oldseed}'];U=(tr['target']-Bt).reshape(len(Bt),-1).astype(float)/24
  phys=raw(tr,'physical_residual');pm,ps,pV,pSpectrum=pca_fit(phys);Xphys=(phys-pm)/ps@pV
  for arm in ['raw_y','decoder_only','physical_residual','same_spectrum_random']:
   path=OUT/'low_rank_correction'/f'{arm}_s{seed}';path.mkdir(parents=True,exist_ok=True)
   if (path/'RESULT.json').exists():records.append(json.loads((path/'RESULT.json').read_text()));continue
   timer=time.monotonic()
   if arm=='same_spectrum_random':
    rt=tr['y'].astype(float);rv=va['y'].astype(float);rmean=rt.mean(0);rstd=rt.std(0);rstd[rstd<1e-8]=1;R=rng.normal(size=(512,128))/np.sqrt(512);Z=(rt-rmean)/rstd@R;zmean=Z.mean(0);Z-=zmean;ev,V=np.linalg.eigh(Z.T@Z/len(Z));W=(V/np.sqrt(np.maximum(ev,1e-10)))@V.T;factor=W@np.diag(pSpectrum);Xt=Z@factor;Xv=(((rv-rmean)/rstd@R)-zmean)@factor
    feat=dict(raw_mean=rmean,raw_scale=rstd,random_projection=R,projection_mean=zmean,spectrum_factor=factor,physical_spectrum=pSpectrum)
   else:
    rt=raw(tr,arm);rv=raw(va,arm);mean,std,V,s=pca_fit(rt);Xt=(rt-mean)/std@V;Xv=(rv-mean)/std@V;feat=dict(raw_mean=mean,raw_scale=std,projection=V,feature_singular_std=s)
   candidates=[];best=None;bestpred=None;bestfit=None
   for ridge in RIDGES:
    fit=covariance_fit(Xt,U,ridge);oof={k:np.zeros_like(U) for k in RANKS}
    for f in range(5):
     keep=fold!=f;ff=covariance_fit(Xt[keep],U[keep],ridge)
     for k in RANKS:oof[k][~keep]=prediction(ff,Xt[~keep],k)
    for k in RANKS:
     correction=prediction(fit,Xv,k)*24;pred=Bv+correction.reshape(Bv.shape);full=full_material(pred,g);risk=np.mean(abs(full-va['truth'])**2,1);cr=np.mean(abs(passive_box(full,g['config'],g['roi'])-va['truth'])**2,1)
     row=dict(rank=k,ridge=ridge,validation_full_mse=float(risk.mean()),validation_clipped_mse=float(cr.mean()),OOF_baseline_full_mse=float(np.mean(np.sum(U**2,1))),OOF_corrected_full_mse=float(np.mean(np.sum((U-oof[k])**2,1))),train_singular_values=fit[4].tolist(),train_topk_singular_energy=float(np.sum(fit[4][:k]**2)))
     candidates.append(row)
     if best is None or row['validation_full_mse']<best['validation_full_mse']:best=row;bestpred=full;bestfit=fit
   xm,um,W,P,s,Q=bestfit;k=best['rank'];coeff=W@Q[:k].T;output=(P[:,:k]*s[:k]).T
   np.savez_compressed(path/'BEST.npz',input_mean=xm,output_mean=um,input_factor=coeff,output_factor=output,rank=k,ridge=best['ridge'],arm=arm,seed=seed,baseline_old_seed=oldseed,**feat)
   np.savez_compressed(path/'validation.npz',parent=va['parent'],prediction=bestpred,risk=np.mean(abs(bestpred-va['truth'])**2,1))
   r=dict(status='COMPLETE',arm=arm,seed=seed,baseline_seed=oldseed,method='train crosscov RRR, not GT covariance PCA',crossfit='5fold OOF correction fit; frozen R4 baseline itself trained on original TRAIN1024',normalization='TRAIN-only means scales PCA; no labels for input frame',train_size=1024,validation_size=256,selected=best,candidates=candidates,linear_parameter_count=int(coeff.size+output.size+xm.size+um.size),train_wall_seconds=time.monotonic()-timer,physics_RHS=0,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
   (path/'RESULT.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:v for k,v in r.items() if k not in ['candidates','selected']}),flush=True);records.append(r)
 (OUT/'RRR_SUMMARY.json').write_text(json.dumps(dict(status='COMPLETE',records=[dict(arm=r['arm'],seed=r['seed'],selected=r['selected'],wall_seconds=r['train_wall_seconds']) for r in records],wall_seconds=time.monotonic()-start,physics_RHS=0),indent=2)+'\n')
if __name__=='__main__':main()
