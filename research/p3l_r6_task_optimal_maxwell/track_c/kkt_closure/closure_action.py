"""Approximate KKT-closed feature algorithm; frozen provider/head, no physics.

Positive stationarity closure for residual y-F is g=z. This substitutes only the
normalized backprojection direction; it never constructs or modifies a Jacobian.
"""
import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
import sys
from pathlib import Path
import numpy as np
import torch
OUT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(OUT));sys.path.insert(0,str(OUT.parents[1].parent/'b1-p3l-r4-conditional-transport/src'))
from track_c.frozen_heads import load_action
from p1_material.prior import saturation_prime
from p3r4.risk import roi_components,full_material

class Action:
 def __init__(self,method='deterministic_matched',seed=2026102601,closure_mode='q',noise_ratio=1.,_base=None):
  if method not in ['deterministic_matched','distilled_w128_n32768']:raise ValueError('unchanged R5 frozen method required')
  if closure_mode not in ['q','zero','negativeq']:raise ValueError('q/zero/negativeq required')
  if not np.isfinite(noise_ratio) or noise_ratio<=0:raise ValueError('known sigma_new/sigma_feature_training must be positive')
  torch.set_num_threads(1);self.method=method;self.seed=seed;self.closure_mode=closure_mode;self.noise_ratio=float(noise_ratio);self.base=load_action(method,seed).action if _base is None else _base;self.p=self.base.p;self.geometry=self.base.geometry;self.fp=self.p.get('feature_payload',self.p);self.normalization=self.fp['feature_normalization'];self.cfg=self.geometry['config'];self.A=np.asarray(self.geometry['A'],np.float32);self.mask=np.zeros(576,np.float32);self.mask[self.geometry['roi']]=1;self.aw=self.A@np.asarray(self.geometry['weak']);self.ast=self.A@np.asarray(self.geometry['strong'])
 def raw_image(self,online,gradient=None):
  q=np.asarray(online['q0'],dtype=np.float32);n=len(q);images=[]
  for i in range(n):
   f=self.cfg['field_mean']+self.A@q[i,:64];df=saturation_prime(f/self.cfg['phase_transition_halfwidth'])/self.cfg['phase_transition_halfwidth']*self.mask
   delta=complex(self.cfg['inclusion_increment_real_interval'][0]+np.diff(self.cfg['inclusion_increment_real_interval'])[0]*q[i,64],self.cfg['inclusion_increment_imag_interval'][0]+np.diff(self.cfg['inclusion_increment_imag_interval'])[0]*q[i,65]);dw=(df[:,None]*self.aw)*abs(delta);ds=(df[:,None]*self.ast)*abs(delta)
   if gradient is not None:g=np.asarray(gradient[i],np.float32)
   elif self.closure_mode=='zero':g=np.zeros(64,np.float32)
   else:g=self.noise_ratio*q[i,:64]*(1 if self.closure_mode=='q' else -1)
   direction=g/max(np.linalg.norm(g),1.);back=delta*df*(self.A@direction);anchor=online['chi0'][i]
   images.append(np.stack((anchor.real,anchor.imag,np.sqrt(np.sum(dw*dw,axis=1)),np.sqrt(np.sum(ds*ds,axis=1)),back.real,back.imag)).reshape(6,24,24))
  return np.asarray(images,np.float32)
 def image(self,online,gradient=None):
  image=self.raw_image(online,gradient);s=self.normalization;image=(image-s['image_mean'])/s['image_scale'];q=np.asarray(online['q0'],np.float32);u=(q[:,-2:]-s['base_mean'][:,-2:])/s['base_scale'][:,-2:];return np.concatenate((image,np.broadcast_to(u[:,:,None,None],(len(q),2,24,24))),1).astype(np.float32)
 def predict(self,y,online,parents=None,batch_size=64):
  y=np.asarray(y)
  if y.dtype!=np.float32 or y.ndim!=2 or y.shape[1]!=512 or not np.isfinite(y).all():raise ValueError('original float32 source-major measurement [N,512] required')
  if np.asarray(online['q0']).shape!=(len(y),66) or np.asarray(online['chi0']).shape!=(len(y),576):raise ValueError('original native-coordinate q0/chi0 shape required')
  rows=[]
  for b in range(0,len(y),batch_size):
   sl=slice(b,b+batch_size);aa={k:online[k][sl] for k in ['q0','chi0']};im=self.image(aa);mom=self.base.moments(y[sl],None if parents is None else parents[sl]);mom=((mom-self.p['moment_mean'])/self.p['moment_scale']).astype(np.float32);anc=roi_components(aa['chi0'],self.geometry['roi'])
   with torch.inference_mode():pred=self.base.head(torch.from_numpy(im),torch.from_numpy(mom),torch.from_numpy(anc)).numpy()
   rows.append(full_material(pred,self.geometry))
  return np.concatenate(rows)
