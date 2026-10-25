import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',PYTORCH_MPS_HIGH_WATERMARK_RATIO='0.12',PYTORCH_MPS_LOW_WATERMARK_RATIO='0.1')
import sys,json,time,copy,resource,traceback
from pathlib import Path
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[2];OLD=ROOT.parent/'b1-p3l-r4-conditional-transport';sys.path.insert(0,str(OLD/'src'))
import numpy as np,torch
from p3r4.risk_fusion import FusionRisk,read_moments
from p3r4.risk import full_material
from p3r2.fusion import passive_box

def main():
 torch.set_num_threads(1);device='mps' if torch.backends.mps.is_available() else 'cpu'
 if device=='mps':torch.mps.set_per_process_memory_fraction(.10)
 g=torch.load(OUT/'geometry.pt',weights_only=False,map_location='cpu')['geometry'];tr=dict(np.load(OUT/'prepared_train.npz'));va=dict(np.load(OUT/'prepared_validation.npz'));S=OLD/'research/p3l_r4_conditional_transport'
 if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss>3*1024**3:raise MemoryError('preallocation RSS gate')
 # Selection sees VAL256 only. Labels are absent from all moment generation routines.
 for arm in ['deterministic_matched','direct_moments_matched']:
  if arm=='deterministic_matched':mt,mv=tr['det'],va['det']
  else:mt=read_moments(S/'inputs/fusion_field_moments/train',tr['parent']);mv=read_moments(S/'inputs/fusion_field_moments/validation',va['parent'])
  mean=mt.mean((0,2,3),keepdims=True,dtype=np.float64).astype('float32');scale=mt.std((0,2,3),keepdims=True,dtype=np.float64).astype('float32');scale[scale<1e-8]=1;mt=(mt-mean)/scale;mv=(mv-mean)/scale
  ti=torch.tensor(tr['image'],device=device);vi=torch.tensor(va['image'],device=device);tm=torch.tensor(mt,device=device);vm=torch.tensor(mv,device=device);ta=torch.tensor(tr['anchor'],device=device);vanc=torch.tensor(va['anchor'],device=device);tt=torch.tensor(tr['target'],device=device);vt=torch.tensor(va['target'],device=device)
  for seed in [2026102601,2026102602,2026102603]:
   path=OUT/'experiments'/f'{arm}_s{seed}';path.mkdir(parents=True,exist_ok=True)
   if (path/'RESULT.json').exists():continue
   start=time.monotonic();cpu=time.process_time();torch.manual_seed(seed);rng=np.random.default_rng(seed);model=FusionRisk(g).to(device);opt=torch.optim.AdamW(model.parameters(),lr=5e-4,weight_decay=1e-3);history=[]
   def evaluate():
    model.eval();preds=[]
    with torch.inference_mode():
     for b in range(0,len(vi),64):preds.append(model(vi[b:b+64],vm[b:b+64],vanc[b:b+64]).cpu().numpy())
    preds=np.concatenate(preds);return float(np.mean((preds-va['target'])**2)),preds
   best,pred=evaluate();state=copy.deepcopy(model.state_dict());bestupdate=0
   with (path/'losses.jsonl').open('w') as log:
    for u in range(1,4001):
     ids=torch.as_tensor(rng.integers(len(ti),size=64),device=device);model.train();opt.zero_grad(set_to_none=True);loss=(model(ti[ids],tm[ids],ta[ids])-tt[ids]).square().mean()
     if not torch.isfinite(loss):raise RuntimeError('nonfinite training loss')
     loss.backward();torch.nn.utils.clip_grad_norm_(model.parameters(),5.);opt.step()
     if u%100==0:
      score,current=evaluate();row=dict(update=u,train_component_mse=float(loss.item()),validation_component_mse=score,wall_seconds=time.monotonic()-start);log.write(json.dumps(row)+'\n');log.flush();history.append(row)
      if score<best:best=score;pred=current;bestupdate=u;state=copy.deepcopy(model.state_dict())
      if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss>3*1024**3:raise MemoryError('3GB resource gate')
   payload=dict(arm=arm,seed=seed,state_dict={k:v.cpu() for k,v in state.items()},geometry=g,moment_mean=mean,moment_scale=scale,selected_update=bestupdate,feature_payload=torch.load(OUT/'geometry.pt',weights_only=False,map_location='cpu')['reference_head_payload'],source_root=str(S),training_parents=tr['parent'],optimizer=dict(lr=5e-4,weight_decay=1e-3,batch_size=64,updates=4000),physics_RHS=0)
   torch.save(payload,path/'BEST.pt');full=full_material(pred,g);clipped=passive_box(full,g['config'],g['roi']);risk=np.mean(abs(full-va['truth'])**2,1);cr=np.mean(abs(clipped-va['truth'])**2,1);np.savez_compressed(path/'validation.npz',parent=va['parent'],prediction=full,risk=risk,clipped_prediction=clipped,clipped_risk=cr)
   r=dict(status='COMPLETE',arm=arm,seed=seed,device=device,parameter_count=sum(p.numel() for p in model.parameters()),provider_parameter_count=682943 if arm=='deterministic_matched' else None,train_size=1024,validation_size=256,selected_update=bestupdate,updates=4000,validation_full_complex_mse=float(risk.mean()),validation_clipped_full_complex_mse=float(cr.mean()),wall_seconds=time.monotonic()-start,cpu_seconds=time.process_time()-cpu,gpu_seconds=time.monotonic()-start if device=='mps' else 0,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,physics_RHS=0,history=history)
   (path/'RESULT.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:v for k,v in r.items() if k!='history'}),flush=True)
  del ti,vi,tm,vm,ta,vanc,tt,vt
  if device=='mps':torch.mps.empty_cache()
if __name__=='__main__':
 status='FAILED';error=None;start=time.monotonic()
 try:main();status='COMPLETE'
 except BaseException as e:error=repr(e);traceback.print_exc();raise
 finally:(OUT/('MATCHED_JOB_CLOSED_'+str(time.time_ns())+'.json')).write_text(json.dumps(dict(status=status,error=error,wall_seconds=time.monotonic()-start,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,physics_RHS=0),indent=2)+'\n')
