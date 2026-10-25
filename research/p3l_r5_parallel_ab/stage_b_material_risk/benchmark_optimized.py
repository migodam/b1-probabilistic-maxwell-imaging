import sys,json,time,resource,os,platform
from pathlib import Path
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[2];OLD=ROOT.parent/'b1-p3l-r4-conditional-transport';sys.path.insert(0,str(OLD/'src'))
import numpy as np,torch
from distill_compact import Compact
from distill_compact_w128 import Compact as Compact128
from material_moments import raw_gaussian,diagonal_image
from p3r4.density import load_density
from p3r4.probability_eval import decode_gaussianized
from p3r4.risk_fusion import FusionRisk,parent_seed
from p3r4.fusion_provider_compat import load_fusion_compatible
from p3r4.risk import feature_image,roi_components
from p3r2.residual_risk import read_anchor
S=OLD/'research/p3l_r4_conditional_transport';torch.set_num_threads(1);d=dict(np.load(OUT/'prepared_validation.npz'));anchor=read_anchor(S/'inputs/strong_context_1024_256/anchor_validation.npz');decoder=json.loads((OUT/'ASSET_INVENTORY.json').read_text())['decoder'];model=load_density(S/'models/density/extended_diag_gaussian_s2026102501/BEST.pt');original=load_fusion_compatible(S/'models/risk/riskfusion_diag_gaussian_s2026102501/BEST.pt');rows=[];equivalence=[]
def K32(y,parents):
 z=[]
 for parent in parents:
  gen=torch.Generator(device='cpu').manual_seed(parent_seed(parent,20261025121));z.append(torch.cat([torch.randn((8,66),generator=gen) for _ in range(4)]))
 z=torch.stack(z)
 with torch.inference_mode():
  cy=torch.from_numpy(model._y(y));m,ls=model.model.head(model.model.encoder(cy)).chunk(2,-1);sx=(m[:,None]+torch.exp(ls.clamp(-5,3))[:,None]*z).numpy()
 st=model.stats;x=(sx*st['x_scale']+st['x_mean'])@st['frame_inverse'];chi,_,_=decode_gaussianized(x.reshape(-1,66),decoder);chi=chi.reshape(len(y),32,576);return np.stack((chi.real.mean(1),chi.imag.mean(1),chi.real.var(1),chi.imag.var(1)),1).reshape(-1,4,24,24).astype('float32')
for start in range(0,256,64):equivalence.append(float(np.max(abs(K32(d['y'][start:start+64],d['parent'][start:start+64])-d['K32'][start:start+64]))))
for label in ['optimized_original_K32','deterministic_q16','distilled_w64_n4096','distilled_w64_n32768','distilled_w128_n32768']:
 cold=time.monotonic()
 if label.startswith('distilled_'):
  width=128 if 'w128' in label else 64;n=4096 if label.endswith('4096') else 32768;folder=f'distilled_w128_n{n}_s2026102601' if width==128 else f'distilled_n{n}_s2026102601';dest=OUT/'learned_context'/folder
  if not (dest/'PROVIDER.pt').exists():continue
  cp=torch.load(dest/'PROVIDER.pt',map_location='cpu',weights_only=False);cm=Compact128() if width==128 else Compact();cm.load_state_dict(cp['state_dict']);cm.eval();p=torch.load(dest/'HEAD.pt',map_location='cpu',weights_only=False);net=FusionRisk(p['geometry']);net.load_state_dict(p['state_dict']);net.eval()
  def mom(y,parents):
   with torch.inference_mode():rr=cm(torch.tensor((y-cp['y_mean'])/cp['y_scale'])).numpy()*cp['target_scale']+cp['target_mean']
   rr[:,2:]=np.maximum(rr[:,2:],0);im=np.zeros((len(y),4,576),np.float32);im[:,:2]=np.asarray(cp['geometry']['config']['host_chi'])[None,:,None];im[:,:,cp['geometry']['roi']]=rr;return im.reshape(-1,4,24,24)
 else:
  if label=='optimized_original_K32':
   model=load_density(S/'models/density/extended_diag_gaussian_s2026102501/BEST.pt');original=load_fusion_compatible(S/'models/risk/riskfusion_diag_gaussian_s2026102501/BEST.pt');p,net=original.p,original.model;mom=K32
  else:
   model=load_density(S/'models/density/extended_diag_gaussian_s2026102501/BEST.pt')
   p=torch.load(OUT/'experiments/deterministic_matched_s2026102601/BEST.pt',map_location='cpu',weights_only=False);net=FusionRisk(p['geometry']);net.load_state_dict(p['state_dict']);net.eval()
   def mom(y,parents):m,C=raw_gaussian(model,y);return diagonal_image(m,C,decoder,16)[0]
 loadtime=time.monotonic()-cold;featurep=p.get('feature_payload',p);g=p['geometry']
 for batch in [1,64]:
  y=d['y'][:batch];parents=d['parent'][:batch];anc={k:anchor[k][:batch] for k in ['observation','residual','q0','chi0','J']};times=[];pt=[]
  def fullcall():
   m=mom(y,parents);image=feature_image(y,g,anc,featurep['feature_normalization']);a=roi_components(anc['chi0'],g['roi']);m=(m-p['moment_mean'])/p['moment_scale']
   with torch.inference_mode():return net(torch.tensor(image),torch.tensor(m),torch.tensor(a)).numpy()
  start=time.monotonic();fullcall();first=time.monotonic()-start
  for rep in range(30):
   start=time.monotonic();mom(y,parents);pt.append(time.monotonic()-start);start=time.monotonic();fullcall();times.append(time.monotonic()-start)
  rows.append(dict(arm=label,batch=batch,checkpoint_load_seconds=loadtime,first_inference_seconds=first,full_warm_median_seconds=float(np.median(times)),full_warm_per_parent_seconds=float(np.median(times)/batch),provider_warm_per_parent_seconds=float(np.median(pt)/batch),provider_p10_seconds=float(np.quantile(pt,.1)),provider_p90_seconds=float(np.quantile(pt,.9)),raw_provider_seconds=pt,raw_full_pipeline_seconds=times,timing_units='secondsperwholebatch; perparent separatelynormalized',full_p10=float(np.quantile(times,.1)),full_p90=float(np.quantile(times,.9))))
res=dict(status='COMPLETE',device='CPU1',environment=dict(platform=platform.platform(),machine=platform.machine(),processor=platform.processor(),logical_CPU=os.cpu_count(),torch=torch.__version__,torch_threads=torch.get_num_threads(),numpy=np.__version__,numeric_environment={k:os.environ.get(k) for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS']}),timing_trials=30,scope='optimized original K32: compute sameGaussianencoder onceperparentbatch, originalchunk8 RNG32draws parentseed preserved; allproviders include fourchannels; fullcost additionally originalallowedfeatureimage+6418head',optimizedK32_vs_storedoriginal4channel_max_abs=max(equivalence),equivalence_tolerance=2e-6,rows=rows,physics_RHS=0,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
(OUT/'OPTIMIZED_COST_BENCHMARK.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res))
