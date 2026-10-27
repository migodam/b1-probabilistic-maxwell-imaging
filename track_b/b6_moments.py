"""Bounded genuine conditional-moment regression; cached TRAIN/VAL physics only."""
import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
import time,json,resource,csv,datetime
from pathlib import Path
import numpy as np
from scipy.special import ndtr
from scipy.linalg import solve
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'research/p3l_r6_task_optimal_maxwell/track_b/b6_moments';R4=ROOT.parent/'b1-p3l-r4-conditional-transport/research/p3l_r4_conditional_transport';R5=ROOT.parent/'b1-p3l-r5-parallel-ab/research/p3l_r5_parallel_ab/stage_b_material_risk';R2=ROOT.parent/'b1-p3l-r2-risk-posterior/research/p3l_r2_risk_posterior_exploration';SIGMA=.00119211798514519
START_CPU=time.process_time();START_WALL=time.perf_counter();FITS=[]
def guard():
 if time.process_time()-START_CPU>108:raise TimeoutError('108 CPU seconds conservative guard, total requested120')
 if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss>4_000_000_000:raise MemoryError('4GB pilot guard')
def write(name,x):(OUT/name).write_text(json.dumps(x,indent=2,default=lambda v:v.item() if isinstance(v,np.generic) else v.tolist())+'\n')
def pack(a):return np.concatenate((a.real.transpose(0,2,1),a.imag.transpose(0,2,1)),2).reshape(len(a),512)
def saturation(s):
 t=np.clip((s+1)/2,0,1);return t*t*t*(10+t*(-15+6*t))
def projections():
 P=np.zeros((512,12));P[[0,32,4*64+16,4*64+32+16],range(4)]=1
 for src in range(8):
  for rec in range(32):
   P[src*64+rec,4]=1/16;P[src*64+32+rec,5]=1/16;P[src*64+rec,6]=(1 if rec<16 else -1)/16;P[src*64+32+rec,7]=(1 if src<4 else -1)/16
 P[:,8:]=np.linalg.qr(np.random.default_rng(202610270627).normal(size=(512,4)))[0];return P
D=json.loads((R5/'ASSET_INVENTORY.json').read_text())['decoder'];ROI=np.array(D['roi']);CFG=D['config'];A=np.array(D['A']);HOST=complex(*CFG['host_chi'])

def load_assets(P):
 cache=OUT/'cached_projection_labels.npz'
 if cache.exists():return dict(np.load(cache)),read_existing()
 data={};records=[]
 for fold,filename in [('train','train32768.npz'),('validation','validation.npz')]:
  with np.load(R4/'inputs/joint'/filename) as j:source=Path(j['source'].item());parent=j['parent']
  with np.load(source) as z:obs=z['observation'];chi=z['chi'];q=z['q'];assert np.array_equal(z['parent'],parent)
  x=np.sqrt(2)*pack(obs);targets=np.column_stack((chi[:,ROI].real.mean(1),chi[:,ROI].imag.mean(1),np.mean(chi[:,ROI].real-HOST.real>.02,1)));F=np.empty((len(parent),12));R=np.empty_like(F);maxobs=0.;maxq=0.;maxnoise=0.;files=0
  for k,name in enumerate(parent):
   if not str(name).startswith('R2_'+fold+'_'):raise RuntimeError('TRAIN/VAL allowlist only')
   path=R2/'inputs/offline_only'/f'{name}.npz'
   if not path.exists():raise FileNotFoundError('MISSING_NOT_RUN no cached noiseless asset: '+str(path))
   with np.load(path) as z:prediction=z['prediction'];noise=z['noise'];qq=z['q']
   maxobs=max(maxobs,float(np.max(abs(prediction+noise-obs[k]))));maxq=max(maxq,float(np.max(abs(qq-q[k]))));fp=np.sqrt(2)*pack(prediction[None])[0];rp=-np.sqrt(2)*pack(noise[None])[0];F[k]=fp@P;R[k]=rp@P;maxnoise=max(maxnoise,float(np.max(abs((fp-x[k])-rp))));files+=1
   if k%512==0:guard()
  assert maxobs==0 and maxq==0
  data.update({fold+'_x':x,fold+'_task':targets,fold+'_f':F,fold+'_r':R,fold+'_parent':parent})
  records.append(dict(fold=fold,n=len(parent),source=str(source),cached_prediction_files=files,source_observation_dtype=str(obs.dtype),direct_prediction_dtype='complex128',max_prediction_plus_noise_observation_error=maxobs,max_cached_q_source_q_error=maxq,max_FminusY_plus_noise_roundoff=maxnoise,read_scope='allowlisted TRAIN/VAL parents only; direct cached prediction+noise, not reconstruction or regeneratedphysics'))
  del obs,chi,q;guard()
 assert not set(data['train_parent'])&set(data['validation_parent']);np.savez_compressed(cache,**data);write('ASSET_RECOVERY.json',dict(status='DIRECT_CACHED_NOISELESS_AVAILABLE',records=records,physics_RHS=0,old_test_or_Dtruth_read=False));return data,records

def read_existing():return json.loads((OUT/'ASSET_RECOVERY.json').read_text())['records']
class Features:
 def __init__(self,kind,x):
  self.kind=kind;self.mean=x.mean(0);self.scale=x.std(0);self.scale[self.scale<1e-10]=1;u=(x-self.mean)/self.scale;rng=np.random.default_rng(202610270628)
  if kind=='raw_y_rbf':
   self.centers=u[np.linspace(0,len(u)-1,127,dtype=int)].copy();self.width=512.
  elif kind=='random_projection':self.proj=np.linalg.qr(rng.normal(size=(512,127)))[0]
  elif kind=='decoder_only':
   self.code_proj=rng.normal(size=(512,66))/np.sqrt(512);code=u@self.code_proj;self.code_scale=code.std(0);self.code_scale[self.code_scale<1e-10]=1;self.material_proj=np.linalg.qr(np.random.default_rng(202610270629).normal(size=(512,127)))[0]
  else:raise ValueError(kind)
  z=self.raw(x);self.zm=z.mean(0);self.zs=z.std(0);self.zs[self.zs<1e-10]=1
 def raw(self,x):
  u=(x-self.mean)/self.scale
  if self.kind=='raw_y_rbf':return np.exp(-np.maximum(np.sum(u*u,1)[:,None]+np.sum(self.centers*self.centers,1)[None]-2*u@self.centers.T,0)/(2*self.width))
  if self.kind=='random_projection':return u@self.proj
  code=(u@self.code_proj)/self.code_scale;field=CFG['field_mean']+code[:,:64]@A[ROI].T;q=saturation(field/CFG['phase_transition_halfwidth']);probit=ndtr(code[:,64:]);real=HOST.real+(CFG['inclusion_increment_real_interval'][0]+np.ptp(CFG['inclusion_increment_real_interval'])*probit[:,0,None])*q;imag=HOST.imag+(CFG['inclusion_increment_imag_interval'][0]+np.ptp(CFG['inclusion_increment_imag_interval'])*probit[:,1,None])*q
  return np.concatenate((real,imag),1)@self.material_proj
 def transform(self,x):return (self.raw(x)-self.zm)/self.zs

def fit(features,x,Y,label):
 tick=time.process_time();Z=features.transform(x);ym=Y.mean(0);M=Z.T@Z/len(Z)+.01*np.eye(127);cond=float(np.linalg.cond(M))
 if not np.isfinite(cond) or cond>1e10:raise np.linalg.LinAlgError('ill-conditioned; no retuning')
 W=solve(M,Z.T@(Y-ym)/len(Z),assume_a='pos');FITS.append(dict(label=label,n=len(x),outputs=Y.shape[1],features=128,ridge=.01,condition=cond,CPU_seconds=time.process_time()-tick));guard();return W,ym

def predict(model,x):
 f,W,ym=model;return f.transform(x)@W+ym

def stats_mse(pred,target,scale):return dict(MSE_per_component=np.mean((pred-target)**2,0),normalized_MSE_per_component=np.mean((pred-target)**2,0)/np.maximum(scale,1e-30),mean_normalized_MSE=float(np.mean(np.mean((pred-target)**2,0)/np.maximum(scale,1e-30))))
def main():
 OUT.mkdir(parents=True,exist_ok=True);P=projections();write('PLAN_FROZEN.json',dict(frozen_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),tasks=['ROImean_material_real','ROImean_material_imag','ROI_fraction_of_cells_with_real_material_increment_gt_.02'],response_coordinate='sqrt2 source-major packed complex response; isotropic coordinate noisevariance sigma²',sigma=SIGMA,projection_sets=['canonical_receivers4','structured_source_receiver_sums4','seeded_random_orthonormal4'],inputs=['raw_y_rbf127 TRAINlandmarks,width512 afterTRAINpercoordinate standardization','random_projection127 fixed orthonormal afterTRAINstandardization','decoder_only: random66D measurementcode→originalC2/probitdecoder→fixed127 materialprojections; uncalibratedpseudo-latent, noanchor/J/latentGT'],features_including_intercept=128,ridge=.01,selection='none; noVALtuning',fold='parentindex parity; eachdirection baselinea0 trainedotherhalf, momentfitownhalf; allfoldnormalizers fitowntrainingonly',moment_targets=['a=Etask','b=Eprojected_noiselessf','c=Etask_times_f','bres=E(f-y)','cres=E((task-a0)*(f-y))'],covariance_estimates=['c-a*b','cres-(a-a0)*bres'],identity_scope='Cov/sigma² identity onlyfortrueposterior; learnedFDgap isregressionconsistency,not biasproof',FD='first128historicalVAL,12fixedprojectiondirections,steps sigma*1e-3 andsigma*1e-4',CPU_limit_seconds=120,physics_RHS=0,noMPS=True))
 data,assets=load_assets(P);tx=data['train_x'];vx=data['validation_x'];tt=data['train_task'];vt=data['validation_task'];tf=data['train_f'];vf=data['validation_f'];tr=data['train_r'];vr=data['validation_r'];fold=np.array([int(str(x).rsplit('_',1)[1])%2 for x in data['train_parent']]);np.savez_compressed(OUT/'FIXED_PROJECTIONS_AND_FOLDS.npz',P=P,train_parent=data['train_parent'],fold=fold,validation_parent=data['validation_parent'])
 summaries=[];predictions={};fdrecords=[]
 for kind in ['raw_y_rbf','random_projection','decoder_only']:
  models=[];foldcov=[];foldres=[];allpred=[];msefold=[]
  for k in [0,1]:
   ids=np.flatnonzero(fold==k);other=np.flatnonzero(fold!=k);fb=Features(kind,tx[other]);W0,m0=fit(fb,tx[other],tt[other],kind+f'_fold{k}_a0_otherhalf');base=(fb,W0,m0);a0t=predict(base,tx[ids]);a0v=predict(base,vx)
   targets=np.concatenate((tt[ids],tf[ids],(tt[ids,:,None]*tf[ids,None,:]).reshape(len(ids),-1),tr[ids],((tt[ids]-a0t)[:,:,None]*tr[ids,None,:]).reshape(len(ids),-1)),1)
   validtargets=np.concatenate((vt,vf,(vt[:,:,None]*vf[:,None,:]).reshape(len(vx),-1),vr,((vt-a0v)[:,:,None]*vr[:,None,:]).reshape(len(vx),-1)),1)
   features=Features(kind,tx[ids]);W,ym=fit(features,tx[ids],targets,kind+f'_fold{k}_moments');model=(features,W,ym);pr=predict(model,vx);a,b,c,br,cr=pr[:,:3],pr[:,3:15],pr[:,15:51].reshape(-1,3,12),pr[:,51:63],pr[:,63:].reshape(-1,3,12);cov=c-a[:,:,None]*b[:,None,:];cres=cr-(a-a0v)[:,:,None]*br[:,None,:];foldcov.append(cov);foldres.append(cres);models.append(model);allpred.append(pr)
   parts={'a':slice(0,3),'b':slice(3,15),'c':slice(15,51),'bres':slice(51,63),'cres':slice(63,99)}
   msefold.append(dict(fold=k,records={label:stats_mse(pr[:,sl],validtargets[:,sl],np.var(targets[:,sl],0)) for label,sl in parts.items()}))
   np.savez_compressed(OUT/f'{kind}_fold{k}_ridge.npz',weights=W,intercept=ym,a0_weights=W0,a0_intercept=m0,validation_prediction=pr,validation_target=validtargets,raw_covariance=cov,centered_covariance=cres)
   guard()
  avg=np.mean(allpred,0);raw=np.mean(foldcov,0);centered=np.mean(foldres,0);derivatives=[]
  for eps in [1e-3,1e-4]:
   der=[];x=vx[:128];delta=SIGMA*eps
   for j in range(12):
    plus=sum(predict(m,x+delta*P[:,j])[:,:3] for m in models)/2;minus=sum(predict(m,x-delta*P[:,j])[:,:3] for m in models)/2;der.append((plus-minus)/(2*delta));guard()
   der=np.stack(der,2);derivatives.append(der)
   for setidx,setname in enumerate(['canonical','structured','random']):
    sl=slice(4*setidx,4*setidx+4);fd=der[:,:,sl];rc=raw[:128,:,sl]/SIGMA**2;cc=centered[:128,:,sl]/SIGMA**2
    fdrecords.append(dict(input=kind,projection_set=setname,step_over_sigma=eps,FD_RMS=float(np.sqrt(np.mean(fd**2))),raw_cov_over_sigma2_RMS=float(np.sqrt(np.mean(rc**2))),centered_cov_over_sigma2_RMS=float(np.sqrt(np.mean(cc**2))),raw_consistency_RMSE=float(np.sqrt(np.mean((rc-fd)**2))),centered_consistency_RMSE=float(np.sqrt(np.mean((cc-fd)**2))),raw_relative_gap=float(np.linalg.norm(rc-fd)/max(np.linalg.norm(fd),1e-30)),centered_relative_gap=float(np.linalg.norm(cc-fd)/max(np.linalg.norm(fd),1e-30)),scope='independentVAL128 regressionconsistency only, no posterior derivative groundtruth'))
  fdstepgap=float(np.linalg.norm(derivatives[0]-derivatives[1])/max(np.linalg.norm(derivatives[1]),1e-30));rawterms=np.abs(avg[:,15:51].reshape(-1,3,12))+np.abs(avg[:,:3,None]*avg[:,None,3:15]);condition=rawterms/np.maximum(np.abs(raw),1e-30)
  summary=dict(input=kind,proper_moment_MSE_by_fold=msefold,ensemble_task_MSE=np.mean((avg[:,:3]-vt)**2,0),ensemble_task_normalized_MSE=np.mean((avg[:,:3]-vt)**2,0)/tt.var(0),raw_covariance_cancellation_median=float(np.median(condition)),raw_covariance_cancellation_p95=float(np.quantile(condition,.95)),raw_centered_covariance_RMS_disagreement=float(np.sqrt(np.mean((raw-centered)**2))),FD_step_relative_difference=fdstepgap)
  summaries.append(summary);predictions[kind+'_raw_cov']=raw;predictions[kind+'_centered_cov']=centered;predictions[kind+'_FD']=derivatives[1];write('PARTIAL.json',dict(status='PARTIAL',records=summaries));guard()
 np.savez_compressed(OUT/'COVARIANCE_AND_FD.npz',**predictions)
 write('RESULTS.json',dict(status='BOUNDED_PILOT_COMPLETE_REGRESSION_CONSISTENCY_ONLY',records=summaries,FD=fdrecords,noise_amplification=1/SIGMA**2,max_asset_subtraction_roundoff=max(r['max_FminusY_plus_noise_roundoff'] for r in assets),roundoff_over_sigma2=max(r['max_FminusY_plus_noise_roundoff'] for r in assets)/SIGMA**2,no_posterior_derivative_or_bias_certificate=True,physics_RHS=0,old_test_or_Dtruth_read=False))
 with (OUT/'DIRECTIONAL_CONSISTENCY.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(fdrecords[0]));w.writeheader();w.writerows(fdrecords)
 write('JOB.json',dict(status='COMPLETE',CPU_seconds=time.process_time()-START_CPU,wall_seconds=time.perf_counter()-START_WALL,CPU_threads=1,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,linear_fit_count=len(FITS),fits=FITS,physics_RHS=0,MPS_seconds=0))
 print(json.dumps(dict(status='COMPLETE',CPU_seconds=time.process_time()-START_CPU,records=[dict(input=r['input'],normalized_task_MSE=r['ensemble_task_normalized_MSE'].tolist(),FD_step_gap=r['FD_step_relative_difference']) for r in summaries]),indent=2))
if __name__=='__main__':
 try:main()
 except BaseException as e:
  OUT.mkdir(parents=True,exist_ok=True);write('JOB.json',dict(status='MISSING_NOT_RUN' if isinstance(e,FileNotFoundError) else 'STOPPED',error=repr(e),CPU_seconds=time.process_time()-START_CPU,wall_seconds=time.perf_counter()-START_WALL,linear_fit_count=len(FITS),fits=FITS,physics_RHS=0));raise
