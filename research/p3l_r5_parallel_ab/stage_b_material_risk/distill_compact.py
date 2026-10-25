import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',PYTORCH_MPS_HIGH_WATERMARK_RATIO='0.12',PYTORCH_MPS_LOW_WATERMARK_RATIO='0.1')
import sys,json,time,copy,resource,traceback,argparse
from pathlib import Path
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[2];OLD=ROOT.parent/'b1-p3l-r4-conditional-transport';sys.path.insert(0,str(OLD/'src'))
import numpy as np,torch
from torch import nn
from material_moments import raw_gaussian,diagonal_image
from p3r4.density import load_density
from p3r4.risk_fusion import FusionRisk
from p3r4.risk import full_material
from p3r2.fusion import passive_box
class Compact(nn.Module):
 def __init__(self):
  super().__init__();self.net=nn.Sequential(nn.Linear(512,64),nn.SiLU(),nn.Linear(64,64),nn.SiLU(),nn.Linear(64,1024));nn.init.zeros_(self.net[-1].weight);nn.init.zeros_(self.net[-1].bias)
 def forward(self,x):return self.net(x).reshape(-1,4,256)
def guard():
 if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss>3*1024**3:raise MemoryError('3GB RSS gate')
def teacher(n):
 path=OUT/f'teacher_train{n}.npz'
 if path.exists():return dict(np.load(path))
 start=time.monotonic();S=OLD/'research/p3l_r4_conditional_transport';model=load_density(S/'models/density/extended_diag_gaussian_s2026102501/BEST.pt');decoder=json.loads((OUT/'ASSET_INVENTORY.json').read_text())['decoder'];g=torch.load(OUT/'geometry.pt',weights_only=False,map_location='cpu')['geometry'];roi=g['roi']
 with np.load(S/'inputs/joint/train32768.npz') as z:y=z['y'][:n];parents=z['parent'][:n]
 rows=[]
 for b in range(0,n,32):
  m,C=raw_gaussian(model,y[b:b+32]);image,_=diagonal_image(m,C,decoder,64);rows.append(image.reshape(-1,4,576)[:,:,roi]);guard()
 t=np.concatenate(rows);np.savez_compressed(path,y=y,teacher=t,parent=parents);(OUT/f'TEACHER_GENERATION_{n}.json').write_text(json.dumps(dict(status='COMPLETE',train=n,wall_seconds=time.monotonic()-start,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,physics_RHS=0,truth_labels_read=False,scope='allowlisted y,parent only; raw jointarchive heldmaterial labels notmaterialized'),indent=2)+'\n');return dict(y=y,teacher=t,parent=parents)
def main(n):
 torch.set_num_threads(1);start=time.monotonic();device='mps' if torch.backends.mps.is_available() else 'cpu'
 if device=='mps':torch.mps.set_per_process_memory_fraction(.1)
 g=torch.load(OUT/'geometry.pt',weights_only=False,map_location='cpu')['geometry'];tr=dict(np.load(OUT/'prepared_train.npz'));va=dict(np.load(OUT/'prepared_validation.npz'));t=teacher(n);ym=t['y'].mean(0);ys=t['y'].std(0);ys[ys<1e-8]=1;tm=t['teacher'].mean((0,2),keepdims=True);ts=t['teacher'].std((0,2),keepdims=True);ts[ts<1e-8]=1;vy=va['y'];vt=va['det'].reshape(-1,4,576)[:,:,g['roi']]
 X=torch.tensor((t['y']-ym)/ys,device=device);Y=torch.tensor((t['teacher']-tm)/ts,device=device);VX=torch.tensor((vy-ym)/ys,device=device);VY=torch.tensor((vt-tm)/ts,device=device);guard();results=[]
 for seed in [2026102601,2026102602,2026102603]:
  dest=OUT/'learned_context'/f'distilled_n{n}_s{seed}';dest.mkdir(parents=True,exist_ok=True)
  if (dest/'RESULT.json').exists():results.append(json.loads((dest/'RESULT.json').read_text()));continue
  torch.manual_seed(seed);rng=np.random.default_rng(seed);model=Compact().to(device);opt=torch.optim.AdamW(model.parameters(),lr=5e-4,weight_decay=1e-3);best=float('inf');bestup=0;timer=time.monotonic()
  with (dest/'provider_losses.jsonl').open('w') as f:
   for u in range(1,4001):
    ids=torch.as_tensor(rng.integers(n,size=128),device=device);model.train();opt.zero_grad(set_to_none=True);loss=(model(X[ids])-Y[ids]).square().mean();loss.backward();torch.nn.utils.clip_grad_norm_(model.parameters(),5);opt.step()
    if u%100==0:
     model.eval()
     with torch.inference_mode():val=float((model(VX)-VY).square().mean().item())
     f.write(json.dumps(dict(update=u,train_teacher_MSE=float(loss.item()),validation_teacher_MSE=val,wall_seconds=time.monotonic()-timer))+'\n');f.flush();guard()
     if val<best:best=val;bestup=u;state=copy.deepcopy(model.state_dict())
  providersecs=time.monotonic()-timer;model.load_state_dict(state);model.eval();predmom=[]
  for dd in [tr,va]:
   pp=[]
   with torch.inference_mode():
    for b in range(0,len(dd['y']),128):pp.append(model(torch.tensor((dd['y'][b:b+128]-ym)/ys,device=device)).cpu().numpy()*ts+tm)
   roiarr=np.concatenate(pp);roiarr[:,2:]=np.maximum(roiarr[:,2:],0);im=np.zeros((len(roiarr),4,576),np.float32);im[:,:2]=np.asarray(g['config']['host_chi'])[None,:,None];im[:,:,g['roi']]=roiarr;predmom.append(im.reshape(-1,4,24,24))
  torch.save(dict(state_dict={k:v.cpu() for k,v in state.items()},y_mean=ym,y_scale=ys,target_mean=tm,target_scale=ts,geometry=g,training_parents=t['parent'],seed=seed,train_size=n,teacher_validation_MSE=best,selected_update=bestup,parameter_count=sum(p.numel() for p in model.parameters()),format='R5_compact_distilled_gaussian_moments_v1',physics_RHS=0),dest/'PROVIDER.pt')
  mt,mv=predmom;mm=mt.mean((0,2,3),keepdims=True,dtype=float).astype('float32');ms=mt.std((0,2,3),keepdims=True,dtype=float).astype('float32');ms[ms<1e-8]=1;mt=(mt-mm)/ms;mv=(mv-mm)/ms;head=FusionRisk(g).to(device);torch.manual_seed(seed);head=FusionRisk(g).to(device);opt=torch.optim.AdamW(head.parameters(),lr=5e-4,weight_decay=1e-3);timer=time.monotonic();ti=torch.tensor(tr['image'],device=device);vi=torch.tensor(va['image'],device=device);mt=torch.tensor(mt,device=device);mv=torch.tensor(mv,device=device);ta=torch.tensor(tr['anchor'],device=device);a=torch.tensor(va['anchor'],device=device);target=torch.tensor(tr['target'],device=device)
  def evaluate():
   head.eval();ps=[]
   with torch.inference_mode():
    for b in range(0,len(vi),64):ps.append(head(vi[b:b+64],mv[b:b+64],a[b:b+64]).cpu().numpy())
   ps=np.concatenate(ps);return float(np.mean((ps-va['target'])**2)),ps
  bval,bpred=evaluate();bs=copy.deepcopy(head.state_dict());bu=0
  with (dest/'head_losses.jsonl').open('w') as f:
   for u in range(1,4001):
    ids=torch.as_tensor(rng.integers(1024,size=64),device=device);head.train();opt.zero_grad(set_to_none=True);loss=(head(ti[ids],mt[ids],ta[ids])-target[ids]).square().mean();loss.backward();torch.nn.utils.clip_grad_norm_(head.parameters(),5);opt.step()
    if u%100==0:
     score,ps=evaluate();f.write(json.dumps(dict(update=u,train_material_component_MSE=float(loss.item()),validation_material_component_MSE=score,wall_seconds=time.monotonic()-timer))+'\n');f.flush();guard()
     if score<bval:bval=score;bpred=ps;bs=copy.deepcopy(head.state_dict());bu=u
  torch.save(dict(state_dict={k:v.cpu() for k,v in bs.items()},geometry=g,moment_mean=mm,moment_scale=ms,provider_path=str(dest/'PROVIDER.pt'),seed=seed,selected_update=bu,feature_payload=torch.load(OUT/'geometry.pt',weights_only=False,map_location='cpu')['reference_head_payload']),dest/'HEAD.pt');full=full_material(bpred,g);cr=passive_box(full,g['config'],g['roi']);risk=np.mean(abs(full-va['truth'])**2,1);np.savez_compressed(dest/'validation.npz',parent=va['parent'],prediction=full,risk=risk)
  r=dict(status='COMPLETE',revision='same4threpresentationarm_compact_distillation',seed=seed,provider_train=n,head_train=1024,validation=256,provider_parameters=sum(p.numel() for p in model.parameters()),head_parameters=6418,provider_selected_update=bestup,head_selected_update=bu,teacher_validation_normalized_MSE=best,validation_full_complex_MSE=float(risk.mean()),validation_clipped_MSE=float(np.mean(abs(cr-va['truth'])**2)),provider_training_seconds=providersecs,head_training_seconds=time.monotonic()-timer,device=device,GPU_seconds=providersecs+time.monotonic()-timer if device=='mps' else 0,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,physics_RHS=0)
  (dest/'RESULT.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True);results.append(r);del ti,vi,mt,mv,ta,a,target
  if device=='mps':torch.mps.empty_cache()
 (OUT/f'DISTILL_SUMMARY_n{n}.json').write_text(json.dumps(dict(status='COMPLETE',records=results,mean_validation_full_MSE=float(np.mean([r['validation_full_complex_MSE'] for r in results])),wall_seconds=time.monotonic()-start,physics_RHS=0),indent=2)+'\n')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--train-size',type=int,choices=[4096,32768],default=4096);args=p.parse_args();status='FAILED';error=None;start=time.monotonic()
 try:main(args.train_size);status='COMPLETE'
 except BaseException as e:error=repr(e);traceback.print_exc();raise
 finally:(OUT/f'DISTILL_JOB_CLOSED_n{args.train_size}.json').write_text(json.dumps(dict(status=status,error=error,controlled_wall_seconds=time.monotonic()-start,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,physics_RHS=0),indent=2)+'\n')
