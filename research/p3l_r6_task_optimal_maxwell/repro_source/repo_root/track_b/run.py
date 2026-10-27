"""R6 finite-task provider learning; immutable historical TRAIN/VAL only, zero physics."""
import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',PYTORCH_MPS_HIGH_WATERMARK_RATIO='0.12',PYTORCH_MPS_LOW_WATERMARK_RATIO='0.1')
import sys,json,time,copy,resource,traceback,argparse
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'research/p3l_r6_task_optimal_maxwell/track_b';HIST=ROOT.parent/'b1-p3l-r5-parallel-ab/research/p3l_r5_parallel_ab/stage_b_material_risk';OLD=ROOT.parent/'b1-p3l-r4-conditional-transport';sys.path.insert(0,str(OLD/'src'))
import numpy as np,torch
from torch import nn
from p3r4.risk_fusion import FusionRisk
from p3r4.risk import full_material
SEEDS=[2026102601,2026102602,2026102603];ARMS=['feature_control1024','direct_gt','feature_action','gt_action']
class Compact(nn.Module):
 def __init__(self):
  super().__init__();self.net=nn.Sequential(nn.Linear(512,128),nn.SiLU(),nn.Linear(128,128),nn.SiLU(),nn.Linear(128,1024));nn.init.zeros_(self.net[-1].weight);nn.init.zeros_(self.net[-1].bias)
 def forward(self,x):return self.net(x).reshape(-1,4,256)
def write(path,obj):path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(obj,indent=2)+'\n')
def guard(start):
 rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
 if rss>5_000_000_000:raise MemoryError('5GB RSS cap')
 if time.monotonic()-start>3500:raise TimeoutError('conservative <1 device-hour wall cap')
def setup(device):
 torch.set_num_threads(1)
 if device=='mps':torch.mps.set_per_process_memory_fraction(.10)
 p=torch.load(HIST/'geometry.pt',weights_only=False,map_location='cpu');g=p['geometry'];tr=dict(np.load(HIST/'prepared_train.npz'));va=dict(np.load(HIST/'prepared_validation.npz'));t=dict(np.load(HIST/'teacher_train32768.npz'))
 with np.load(OLD/'research/p3l_r4_conditional_transport/inputs/joint/train32768.npz') as z:
  assert np.array_equal(z['parent'],t['parent']);assert np.array_equal(z['y'],t['y']);joint_material=z['material']
 ids=np.array([{str(x):i for i,x in enumerate(t['parent'])}[str(v)] for v in tr['parent']]);assert np.array_equal(t['y'][ids],tr['y']);assert not set(tr['parent'])&set(va['parent']);assert np.array_equal(joint_material[ids],tr['truth'])
 del joint_material
 ym=t['y'].mean(0);ys=t['y'].std(0);ys[ys<1e-8]=1;tm=t['teacher'].mean((0,2),keepdims=True);ts=t['teacher'].std((0,2),keepdims=True);ts[ts<1e-8]=1
 arrays={'X':(tr['y']-ym)/ys,'VX':(va['y']-ym)/ys,'Y':(t['teacher'][ids]-tm)/ts,'VY':(va['det'].reshape(-1,4,576)[:,:,g['roi']]-tm)/ts,'image':tr['image'],'vimage':va['image'],'anchor':tr['anchor'],'vanchor':va['anchor'],'target':tr['target'],'vtarget':va['target']}
 d={k:torch.tensor(v,dtype=torch.float32,device=device) for k,v in arrays.items()}
 meta=dict(train_joint=len(t['y']),train_action=len(tr['y']),validation=len(va['y']),train_ids=ids.tolist(),train_parents=tr['parent'].tolist(),validation_parents=va['parent'].tolist(),parameter_count=sum(v.numel() for v in Compact().parameters()),device=device,source_paths=[str(HIST/x) for x in ['teacher_train32768.npz','prepared_train.npz','prepared_validation.npz','geometry.pt']],physics_RHS=0,no_locked_or_old_test_read=True,normalization='y and teacher moment statistics historical TRAIN32768, no GT beyond historical TRAIN/VAL')
 return g,tr,va,t,p,d,ym,ys,tm,ts,meta

def embed(roi,g):
 im=torch.zeros((len(roi),4,576),device=roi.device);im[:,:2]=torch.as_tensor(g['config']['host_chi'],device=roi.device)[None,:,None];im[:,:,g['roi']]=torch.cat((roi[:,:2],roi[:,2:].clamp(min=0)),1);return im.reshape(-1,4,24,24)
def frozen_head(seed,g,device):
 s=seed-100;path=OLD/f'research/p3l_r4_conditional_transport/models/risk/riskfusion_diag_gaussian_s{s}/BEST.pt';p=torch.load(path,weights_only=False,map_location='cpu');h=FusionRisk(g).to(device);h.load_state_dict(p['state_dict']);h.eval()
 for w in h.parameters():w.requires_grad_(False)
 return h,torch.tensor(p['moment_mean'],device=device),torch.tensor(p['moment_scale'],device=device),str(path)
def infer(h,im,mom,anc):
 with torch.no_grad():return torch.cat([h(im[b:b+64],mom[b:b+64],anc[b:b+64]) for b in range(0,len(im),64)])
def risks(pred,va,g):
 full=full_material(pred,g);risk=np.mean(abs(full-va['truth'])**2,1);return full,risk

def run_arm(arm,seed,data,jobstart,device):
 g,tr,va,t,p,d,ym,ys,tm,ts,meta=data;dest=OUT/'models'/f'{arm}_s{seed}';dest.mkdir(parents=True,exist_ok=True)
 if (dest/'RESULT.json').exists():return json.loads((dest/'RESULT.json').read_text())
 start=time.monotonic();torch.manual_seed(seed);rng=np.random.default_rng(seed);model=Compact().to(device);h,mm,ms,hpath=frozen_head(seed,g,device);tm_t=torch.tensor(tm,device=device);ts_t=torch.tensor(ts,device=device);opt=torch.optim.AdamW(model.parameters(),lr=5e-4,weight_decay=1e-3)
 with torch.no_grad():ta= infer(h,d['image'],(torch.tensor(tr['det'],device=device)-mm)/ms,d['anchor']);va_teacher=infer(h,d['vimage'],(torch.tensor(va['det'],device=device)-mm)/ms,d['vanchor'])
 action_scale=float(d['target'].var(dim=0,unbiased=False).mean().item())
 def provider_metrics():
  model.eval()
  with torch.no_grad():
   roi=model(d['VX'])*ts_t+tm_t;mom=(embed(roi,g)-mm)/ms;pred=infer(h,d['vimage'],mom,d['vanchor']);feat=float(((roi-torch.tensor(va['det'].reshape(-1,4,576)[:,:,g['roi']],device=device))/ts_t).square().mean().item());score=float((pred-d['vtarget']).square().mean().item());ad=float((pred-va_teacher).square().mean().item())
  return score,feat,ad,pred.cpu().numpy()
 best,_,_,bp=provider_metrics();bu=0;bs=copy.deepcopy(model.state_dict());first=best
 with (dest/'provider_losses.jsonl').open('w') as f:
  f.write(json.dumps(dict(update=0,validation_material_component_MSE=best))+'\n')
  for u in range(1,4001):
   ids=torch.as_tensor(rng.integers(1024,size=64),device=device);model.train();opt.zero_grad(set_to_none=True);raw=model(d['X'][ids]);feat=(raw-d['Y'][ids]).square().mean();roi=raw*ts_t+tm_t
   if arm=='feature_control1024':loss=feat;gt=act=torch.zeros((),device=device)
   else:
    pred=h(d['image'][ids],(embed(roi,g)-mm)/ms,d['anchor'][ids]);gt=(pred-d['target'][ids]).square().mean();act=(pred-ta[ids]).square().mean()
    loss={'direct_gt':gt/action_scale,'feature_action':feat+act/action_scale,'gt_action':(gt+act)/action_scale}[arm]
   loss.backward();torch.nn.utils.clip_grad_norm_(model.parameters(),5);opt.step()
   if u%100==0:
    val,vf,ad,pp=provider_metrics();f.write(json.dumps(dict(update=u,loss=float(loss.item()),train_feature_MSE=float(feat.item()),train_GT_action_MSE=float(gt.item()),train_teacher_action_MSE=float(act.item()),validation_material_component_MSE=val,validation_teacher_feature_MSE=vf,validation_teacher_action_distortion=ad,wall_seconds=time.monotonic()-start))+'\n');f.flush();guard(jobstart)
    if val<best:best=val;bu=u;bs=copy.deepcopy(model.state_dict());bp=pp
 providerseconds=time.monotonic()-start;model.load_state_dict(bs);model.eval();mom=[]
 with torch.no_grad():
  for x in [d['X'],d['VX']]:mom.append(embed(model(x)*ts_t+tm_t,g).detach())
 full,fr=risks(bp,va,g);np.savez_compressed(dest/'frozen_validation.npz',parent=va['parent'],prediction=full,risk=fr)
 torch.save(dict(state_dict={k:v.cpu() for k,v in bs.items()},y_mean=ym,y_scale=ys,target_mean=tm,target_scale=ts,geometry=g,seed=seed,selected_update=bu,arm=arm,format='R6_functional_compact_v1',parameter_count=214272,training_parents=tr['parent'],normalization_parents=t['parent'],frozen_head_path=hpath,physics_RHS=0),dest/'PROVIDER.pt')
 mt,mv=mom;hm=mt.mean((0,2,3),keepdim=True);hs=mt.std((0,2,3),keepdim=True,unbiased=False);hs=torch.where(hs<1e-8,torch.ones_like(hs),hs);mt=(mt-hm)/hs;mv=(mv-hm)/hs
 torch.manual_seed(seed);head=FusionRisk(g).to(device);opt=torch.optim.AdamW(head.parameters(),lr=5e-4,weight_decay=1e-3);rng=np.random.default_rng(seed);ht=time.monotonic()
 def evaluate():
  pp=infer(head,d['vimage'],mv,d['vanchor']);return float((pp-d['vtarget']).square().mean().item()),pp.cpu().numpy()
 besth,bpred=evaluate();bhs=copy.deepcopy(head.state_dict());bhu=0
 with (dest/'head_losses.jsonl').open('w') as f:
  f.write(json.dumps(dict(update=0,validation_material_component_MSE=besth))+'\n')
  for u in range(1,4001):
   ids=torch.as_tensor(rng.integers(1024,size=64),device=device);head.train();opt.zero_grad(set_to_none=True);loss=(head(d['image'][ids],mt[ids],d['anchor'][ids])-d['target'][ids]).square().mean();loss.backward();torch.nn.utils.clip_grad_norm_(head.parameters(),5);opt.step()
   if u%100==0:
    val,pp=evaluate();f.write(json.dumps(dict(update=u,train_material_component_MSE=float(loss.item()),validation_material_component_MSE=val,wall_seconds=time.monotonic()-ht))+'\n');f.flush();guard(jobstart)
    if val<besth:besth=val;bpred=pp;bhs=copy.deepcopy(head.state_dict());bhu=u
 headseconds=time.monotonic()-ht
 torch.save(dict(state_dict={k:v.cpu() for k,v in bhs.items()},geometry=g,moment_mean=hm.cpu().numpy(),moment_scale=hs.cpu().numpy(),seed=seed,selected_update=bhu,provider_path=str(dest/'PROVIDER.pt'),feature_payload=p['reference_head_payload'],arm=arm,format='R6_functional_matched_head_v1'),dest/'HEAD.pt')
 full,mr=risks(bpred,va,g);np.savez_compressed(dest/'validation.npz',parent=va['parent'],prediction=full,risk=mr)
 result=dict(status='COMPLETE',arm=arm,seed=seed,train=1024,validation=256,provider_parameters=214272,head_parameters=6418,provider_updates=4000,head_updates=4000,provider_selected_update=bu,head_selected_update=bhu,initial_provider_material_component_MSE=first,frozen_validation_full_complex_MSE=float(fr.mean()),matched_validation_full_complex_MSE=float(mr.mean()),action_loss_scale=action_scale,provider_training_seconds=providerseconds,head_training_seconds=headseconds,wall_seconds=time.monotonic()-start,device=device,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,MPS_current_allocated_bytes=torch.mps.current_allocated_memory() if device=='mps' else 0,MPS_driver_allocated_bytes=torch.mps.driver_allocated_memory() if device=='mps' else 0,physics_RHS=0)
 write(dest/'RESULT.json',result);write(OUT/'jobs'/f'{arm}_s{seed}.json',result);print(json.dumps(result),flush=True);return result

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--prepare-only',action='store_true');parser.add_argument('--arms',nargs='+',default=ARMS);args=parser.parse_args();device='mps' if torch.backends.mps.is_available() else 'cpu';start=time.monotonic();data=setup(device);meta=data[-1]
 write(OUT/'READINESS.json',dict(status='READY',**meta,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
 plan=dict(status='FROZEN_BEFORE_INITIAL_TRAINING',arms=args.arms,seeds=SEEDS,architecture='512-128-128-1024 SiLU, zero output init',provider_updates=4000,matched_head_updates=4000,batch_size=64,optimizer='AdamW lr5e-4 wd1e-3 clip5',provider_initialization='fresh same seed for every arm',frozen_head='corresponding historical Gaussian fusion head seed provider_seed-100',feature_control1024='normalized teacher feature MSE only',direct_gt='GT material head-output MSE / train target component variance',feature_action='normalized feature MSE + frozen teacher action MSE / same scale',gt_action='(GT material head-output MSE + frozen teacher action MSE) / same scale',feature_data='paired TRAIN1024; original teacher normalization from TRAIN32768 disclosed',same_context_schedule=True,selection='provider by VAL256 frozen-head material MSE incl update0; matched head by VAL256 material MSE incl update0 every100 updates',historical_original_baseline='immutable width128 n32768 feature-distill checkpoints retained separately; different provider training count/selection from causal1024 comparison',resource_cap=dict(MPS_device_seconds=3600,conservative_wall_gate=3500,RSS_bytes=5_000_000_000,CPU_threads=1),physics_RHS=0,posterior_claim=False)
 if not (OUT/'PLAN_FROZEN.json').exists():write(OUT/'PLAN_FROZEN.json',plan)
 if args.prepare_only:return
 records=[];status='FAILED';err=None
 try:
  for arm in args.arms:
   for seed in SEEDS:
    guard(start);records.append(run_arm(arm,seed,data,start,device));torch.mps.empty_cache() if device=='mps' else None
  status='COMPLETE'
 except BaseException as e:err=repr(e);raise
 finally:
  write(OUT/'TRAINING_JOB.json',dict(status=status,error=err,wall_seconds=time.monotonic()-start,records=records,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,device=device,physics_RHS=0))
if __name__=='__main__':main()
