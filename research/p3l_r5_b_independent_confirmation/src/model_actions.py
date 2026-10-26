"""Frozen observation-only CPU1 actions. No physics solve or truth access.
Historical source and checkpoint files are imported/read in place. K32 preserves
four chunk-8 Gaussian draws for each parent but evaluates the encoder once/batch.
"""
from pathlib import Path
import sys
sys.dont_write_bytecode=True
import numpy as np
import torch
ROOT=Path(__file__).resolve().parents[3]
B1=ROOT.parent
HIST=B1/'b1-p3l-r5-parallel-ab/research/p3l_r5_parallel_ab/stage_b_material_risk'
R4=B1/'b1-p3l-r4-conditional-transport'
S=R4/'research/p3l_r4_conditional_transport'
for path in (R4/'src',HIST):
 if str(path) not in sys.path:sys.path.insert(0,str(path))
from inference import Action as HistoricalAction
from p3r4.risk import RiskActionPredictor
from p3r4.risk_fusion import parent_seed
from p3r4.probability_eval import decode_gaussianized
ONLINE_KEYS=('observation','residual','q0','chi0','J')
SEEDS=(2026102601,2026102602,2026102603)
class Action:
 def __init__(self,method='distilled_w128_n32768',seed=2026102601,device='cpu',quadrature_nodes=64):
  if device!='cpu':raise ValueError('confirmation canonical backend is CPU1')
  if torch.get_num_threads()!=1:torch.set_num_threads(1)
  seed=int(seed)
  if seed in (2026102501,2026102502,2026102503):seed+=100
  if seed not in SEEDS:raise ValueError('frozen seed required')
  aliases={'strongCNN':'strongCNN','strong':'strongCNN','strongCNN3seeds':'strongCNN','original_K32':'original_K32_frozen','K32':'original_K32_frozen','original_K32_frozen3seeds':'original_K32_frozen','deterministic_matched3seeds':'deterministic_matched','distilled_w128_n327683seeds':'distilled_w128_n32768'}
  self.method=aliases.get(method,method);self.seed=seed;self.device=device
  if self.method=='physics':return
  if self.method=='strongCNN':
   with torch.random.fork_rng(devices=[]):self.action=RiskActionPredictor(S/f'models/risk/strongrisk_spectral_s{seed-100}/BEST.pt',device)
   self.geometry=self.action.geometry;self.p=self.action.p;return
  if self.method not in ('original_K32_frozen','deterministic_matched','distilled_w128_n32768','det_frozen'):raise ValueError('unsupported frozen action')
  # Constructors initialize tensors before loading weights; isolate their RNG use.
  with torch.random.fork_rng(devices=[]):
   self.action=HistoricalAction('original_K32_frozen' if self.method=='det_frozen' else self.method,seed,device,quadrature_nodes)
  self.geometry=self.action.geometry;self.p=self.action.p
  if self.method=='original_K32_frozen':self.action.moments=self.moments
  elif self.method=='det_frozen':
   from material_moments import raw_gaussian,gaussian_image
   self.action.moments=lambda y,parents=None:gaussian_image(*raw_gaussian(self.action.provider,y),self.action.decoder,quadrature_nodes)
 def moments(self,y,parents=None):
  if self.method!='original_K32_frozen':return self.action.moments(y,parents)
  if parents is None or len(parents)!=len(y):raise ValueError('K32 needs frozen parent identifiers')
  model=self.action.provider;z=[]
  for parent in parents:
   gen=torch.Generator(device='cpu').manual_seed(parent_seed(parent,20261025121))
   z.append(torch.cat([torch.randn((8,66),generator=gen) for _ in range(4)]))
  with torch.inference_mode():
   cy=torch.from_numpy(model._y(np.asarray(y)));m,ls=model.model.head(model.model.encoder(cy)).chunk(2,-1)
   sx=(m[:,None]+torch.exp(ls.clamp(-5,3))[:,None]*torch.stack(z)).numpy()
  st=model.stats;x=(sx*st['x_scale']+st['x_mean'])@st['frame_inverse']
  chi,_,_=decode_gaussianized(x.reshape(-1,66),self.action.decoder);chi=chi.reshape(len(y),32,576)
  return np.stack((chi.real.mean(1),chi.imag.mean(1),chi.real.var(1),chi.imag.var(1)),1).reshape(-1,4,24,24).astype('float32')
 def predict(self,y,online,parents=None,batch_size=64):
  y=np.asarray(y)
  if y.ndim!=2 or y.shape[1]!=512 or np.iscomplexobj(y) or not np.isfinite(y).all():raise ValueError('finite real source-major512 measurement rows required')
  allowed={k:online[k] for k in ONLINE_KEYS}
  if any(len(v)!=len(y) or not np.isfinite(v).all() for v in allowed.values()):raise ValueError('online length or finite-value mismatch')
  if self.method=='physics':return np.asarray(allowed['chi0']).copy()
  if self.method=='strongCNN':return self.action.predict(y,allowed,batch_size)
  return self.action.predict(y,allowed,parents,batch_size)
