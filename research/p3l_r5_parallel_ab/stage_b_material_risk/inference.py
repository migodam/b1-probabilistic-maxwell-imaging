"""Observation-only inference for frozen Track B finalists.
Usage: action=Action('distilled_w128_n32768',2026102601); chi=action.predict(y,online)
Online must contain observation,residual,q0,chi0,J from unchanged R4 strong anchor.
No truth, physics solve, optimizer update, or label access occurs here.
"""
import sys,json
from pathlib import Path
import numpy as np,torch
OUT=Path(__file__).resolve().parent;REPO=OUT.parents[2];OLD=REPO.parent/'b1-p3l-r4-conditional-transport';sys.path.insert(0,str(OLD/'src'))
from material_moments import raw_gaussian,gaussian_image
from p3r4.risk import feature_image,roi_components,full_material
from p3r4.risk_fusion import FusionRisk,parent_seed
from p3r4.density import load_density
from p3r4.fusion_provider_compat import load_fusion_compatible
from p3r4.probability_eval import decode_gaussianized
from distill_compact_w128 import Compact as Compact128
SEEDS=[2026102601,2026102602,2026102603]
class Action:
 def __init__(self,method='distilled_w128_n32768',seed=2026102601,device='cpu',quadrature_nodes=64):
  if seed not in SEEDS:raise ValueError('frozen training seed required')
  if method not in ['distilled_w128_n32768','deterministic_matched','original_K32_frozen']:raise ValueError('frozen method required')
  self.method=method;self.seed=seed;self.device=device;self.quadrature_nodes=quadrature_nodes;self.decoder=json.loads((OUT/'ASSET_INVENTORY.json').read_text())['decoder'];S=OLD/'research/p3l_r4_conditional_transport'
  if method=='original_K32_frozen':
   oldseed=2026102501+SEEDS.index(seed);self.reference=load_fusion_compatible(S/f'models/risk/riskfusion_diag_gaussian_s{oldseed}/BEST.pt',device);self.p=self.reference.p;self.head=self.reference.model
  else:
   directory=OUT/'learned_context'/f'distilled_w128_n32768_s{seed}' if method.startswith('distilled_') else OUT/'experiments'/f'deterministic_matched_s{seed}'
   self.p=torch.load(directory/('HEAD.pt' if method.startswith('distilled_') else 'BEST.pt'),map_location='cpu',weights_only=False);self.head=FusionRisk(self.p['geometry']).to(device);self.head.load_state_dict(self.p['state_dict']);self.head.eval()
   if method.startswith('distilled_'):
    self.cp=torch.load(directory/'PROVIDER.pt',map_location='cpu',weights_only=False);self.provider=Compact128().to(device);self.provider.load_state_dict(self.cp['state_dict']);self.provider.eval()
  self.geometry=self.p['geometry']
  if not method.startswith('distilled_'):self.provider=load_density(S/'models/density/extended_diag_gaussian_s2026102501/BEST.pt',device)
 def moments(self,y,parents=None):
  if self.method.startswith('distilled_'):
   p=self.cp
   with torch.inference_mode():rr=self.provider(torch.tensor((y-p['y_mean'])/p['y_scale'],device=self.device)).cpu().numpy()*p['target_scale']+p['target_mean']
   rr[:,2:]=np.maximum(rr[:,2:],0);im=np.zeros((len(y),4,576),np.float32);im[:,:2]=np.asarray(self.geometry['config']['host_chi'])[None,:,None];im[:,:,self.geometry['roi']]=rr;return im.reshape(-1,4,24,24)
  if self.method=='deterministic_matched':
   m,C=raw_gaussian(self.provider,y);return gaussian_image(m,C,self.decoder,self.quadrature_nodes)
  if parents is None or len(parents)!=len(y):raise ValueError('parent identifiers required for originalK32 frozen samplingseed')
  out=[]
  for yy,parent in zip(y,parents):
   x=self.provider.sample(yy[None],32,parent_seed(parent,20261025121),batch_size=8)[0];chi,_,_=decode_gaussianized(x,self.decoder);out.append(np.stack((chi.real.mean(0),chi.imag.mean(0),chi.real.var(0),chi.imag.var(0))).reshape(4,24,24))
  return np.array(out,dtype='float32')
 def predict(self,y,online,parents=None,batch_size=64):
  y=np.asarray(y)
  if y.ndim!=2 or y.shape[1]!=512 or np.iscomplexobj(y) or not np.isfinite(y).all():raise ValueError('finite real source-major512 measurement rows')
  allowed={k:online[k] for k in ['observation','residual','q0','chi0','J']};p=self.p;fp=p.get('feature_payload',p);out=[]
  for b in range(0,len(y),batch_size):
   yy=y[b:b+batch_size];aa={k:v[b:b+batch_size] for k,v in allowed.items()};m=self.moments(yy,None if parents is None else parents[b:b+batch_size]);im=feature_image(yy,self.geometry,aa,fp['feature_normalization']);anchor=roi_components(aa['chi0'],self.geometry['roi']);m=((m-p['moment_mean'])/p['moment_scale']).astype('float32')
   with torch.inference_mode():pred=self.head(torch.tensor(im,device=self.device),torch.tensor(m,device=self.device),torch.tensor(anchor,device=self.device)).cpu().numpy()
   out.append(full_material(pred,self.geometry))
  return np.concatenate(out) if out else np.empty((0,576),complex)
