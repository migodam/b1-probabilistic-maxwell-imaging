"""Label-free deployment API. Existing y-only online anchor inputs are required."""
from pathlib import Path
import sys
import numpy as np
import torch
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(Path(__file__).resolve().parent));sys.path.insert(0,str(ROOT.parent/'b1-p3l-r4-conditional-transport/src'))
try:
 from .run import Compact,OUT,HIST,OLD
except ImportError:
 from run import Compact,OUT,HIST,OLD
from p3r4.risk_fusion import FusionRisk
from p3r4.risk import feature_image,roi_components,full_material
class Action:
 def __init__(self,method='feature_action',seed=2026102601,head_mode='matched',device='cpu',provider_seed=None):
  if method not in ['original_feature32768','feature_control1024','direct_gt','feature_action','gt_action']:raise ValueError('unknown frozen candidate')
  if head_mode not in ['matched','frozen']:raise ValueError('matched/frozen required')
  provider_seed=seed if provider_seed is None else provider_seed;self.method=method;self.seed=seed;self.provider_seed=provider_seed;self.device=device;self.head_mode=head_mode
  def folder(s):return HIST/'learned_context'/f'distilled_w128_n32768_s{s}' if method=='original_feature32768' else OUT/'models'/f'{method}_s{s}'
  self.p=torch.load(folder(provider_seed)/'PROVIDER.pt',map_location='cpu',weights_only=False);self.g=self.p['geometry'];self.provider=Compact().to(device);self.provider.load_state_dict(self.p['state_dict']);self.provider.eval()
  if head_mode=='matched':self.hp=torch.load(folder(seed)/'HEAD.pt',map_location='cpu',weights_only=False);self.feature_payload=self.hp['feature_payload']
  else:
   self.hp=torch.load(OLD/f'research/p3l_r4_conditional_transport/models/risk/riskfusion_diag_gaussian_s{seed-100}/BEST.pt',map_location='cpu',weights_only=False);self.feature_payload=self.hp
  self.head=FusionRisk(self.g).to(device);self.head.load_state_dict(self.hp['state_dict']);self.head.eval()
 def moments(self,y,batch_size=64):
  rows=[]
  with torch.no_grad():
   for b in range(0,len(y),batch_size):
    x=torch.tensor((np.asarray(y[b:b+batch_size])-self.p['y_mean'])/self.p['y_scale'],dtype=torch.float32,device=self.device);raw=self.provider(x).cpu().numpy()*self.p['target_scale']+self.p['target_mean'];raw[:,2:]=np.maximum(raw[:,2:],0);im=np.zeros((len(raw),4,576),np.float32);im[:,:2]=np.asarray(self.g['config']['host_chi'])[None,:,None];im[:,:,self.g['roi']]=raw;rows.append(im.reshape(-1,4,24,24))
  return np.concatenate(rows)
 def predict_prepared(self,image,anchor,moments,batch_size=64):
  mom=(moments-self.hp['moment_mean'])/self.hp['moment_scale'];rows=[]
  with torch.no_grad():
   for b in range(0,len(image),batch_size):rows.append(self.head(torch.tensor(image[b:b+batch_size],device=self.device),torch.tensor(mom[b:b+batch_size],device=self.device),torch.tensor(anchor[b:b+batch_size],device=self.device)).cpu().numpy())
  return full_material(np.concatenate(rows),self.g)
 def predict(self,y,online,batch_size=64):
  if any(k not in online for k in ['observation','residual','q0','chi0','J']):raise ValueError('original complete y-only online anchor required')
  y=np.asarray(y,dtype=np.float32)
  if y.ndim!=2 or y.shape[1]!=512:raise ValueError('y shape [N,512] required')
  image=feature_image(y,self.g,online,self.feature_payload['feature_normalization']);anchor=roi_components(online['chi0'],self.g['roi']);return self.predict_prepared(image,anchor,self.moments(y,batch_size),batch_size)
