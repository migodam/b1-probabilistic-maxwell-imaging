"""Prepare-only local response/JVP surrogate, distinct from action-error ranker.

Output is the whitened nonlinear secant remainder beyond cached J0*h.
Candidate h is the differentiated input; x0,f0,r0,J0 context is FIXED.
No candidate Full data, GT, T1 or queried J1 enters inference features.
No density, certificate, automatic Full evaluator or optimizer is implied.
"""
import time
import numpy as np
from .core import pack

V_SEED=2026102518
SEEDS=(2026102501,2026102502,2026102503)


def whiten(z,sigma):return pack(np.asarray(z))*np.sqrt(2)/sigma


def unwhiten(w,sigma,receivers=32,sources=8):
 a=np.asarray(w).reshape(sources,2*receivers)
 return (a[:,:receivers]+1j*a[:,receivers:]).T*sigma/np.sqrt(2)


def fixed_directions(dimension):
 rng=np.random.default_rng(np.random.SeedSequence([V_SEED,dimension]))
 return rng.choice(np.array([-1.,1.]),size=(dimension,2))/np.sqrt(dimension)


def assemble_input(x0,h,f0,r0,J0):
 return np.r_[np.asarray(x0),np.asarray(h),np.asarray(f0),np.asarray(r0),np.asarray(J0)@h]


def input_tangent_matrix(J0):
 """d u_raw / d h = [0; I; 0; 0; J0], context held fixed."""
 J0=np.asarray(J0);m,d=J0.shape
 return np.r_[np.zeros((d,d)),np.eye(d),np.zeros((m,d)),np.zeros((m,d)),J0]


def parameter_count(dimension=66,white_dimension=512):
 input_dim=2*dimension+3*white_dimension
 return input_dim*256+256+256*256+256+256*white_dimension+white_dimension+4*256


def build_local_forward_model(dimension,input_mean,input_scale,output_scale):
 """Lazy CPU model definition; preparing data never calls an optimizer."""
 import torch
 from torch import nn
 class LocalResponse(nn.Module):
  def __init__(self):
   super().__init__();m=len(output_scale);input_dim=2*dimension+3*m
   self.net=nn.Sequential(nn.Linear(input_dim,256),nn.LayerNorm(256),nn.SiLU(),
     nn.Linear(256,256),nn.LayerNorm(256),nn.SiLU(),nn.Linear(256,m))
   self.register_buffer('input_mean',torch.as_tensor(input_mean,dtype=torch.float32))
   self.register_buffer('input_scale',torch.as_tensor(input_scale,dtype=torch.float32))
   self.register_buffer('output_scale',torch.as_tensor(output_scale,dtype=torch.float32))
   self.dimension=dimension
  def forward(self,h,x0,f0,r0,J0):
   # Context detachment is deliberate; only candidate delta is differentiated.
   x0=x0.detach();f0=f0.detach();r0=r0.detach();J0=J0.detach()
   linear=torch.bmm(J0,h.unsqueeze(-1)).squeeze(-1)
   raw=torch.cat((x0,h,f0,r0,linear),dim=-1)
   zero=torch.cat((x0,torch.zeros_like(h),f0,r0,torch.zeros_like(f0)),dim=-1)
   return self.output_scale*(self.net((raw-self.input_mean)/self.input_scale)-self.net((zero-self.input_mean)/self.input_scale))
  def prediction(self,h,x0,f0,r0,J0):
   return f0+torch.bmm(J0,h.unsqueeze(-1)).squeeze(-1)+self(h,x0,f0,r0,J0)
 model=LocalResponse()
 if sum(p.numel() for p in model.parameters())>800000:raise ValueError('forward/JVP model0.8M capacity cap')
 return model


def residual_jvp(model,h,context,direction,*,create_graph=False):
 import torch
 function=lambda delta:model(delta,context['x0'],context['f0'],context['r0'],context['J0'])
 _,jvp=torch.autograd.functional.jvp(function,h,direction,create_graph=create_graph,strict=False)
 return jvp


def value_jvp_loss(model,h,context,value_target,jvp_target,V,*,derivative_weight):
 """Equal normalized value/JVP weights for DINO; derivative_weight=0 control.

    V is unit raw-coordinate Rademacher d×2, not normalized-input directions.
    Model performs input/output chain rules before comparison in white units.
 """
 import torch
 prediction=model(h,context['x0'],context['f0'],context['r0'],context['J0'])
 value=(((prediction-value_target)/model.output_scale)**2).mean()
 derivative=torch.zeros((),device=h.device)
 if derivative_weight:
  errors=[]
  for k in range(2):
   direction=V[:,k][None,:].expand_as(h)
   predicted=residual_jvp(model,h,context,direction,create_graph=True)
   errors.append((((predicted-jvp_target[:,:,k])/model.output_scale)**2).mean())
  derivative=torch.stack(errors).mean()
 return value+float(derivative_weight)*derivative,dict(value=value,derivative=derivative)


def timed_single_prediction(model,h,context):
 """Single candidate assembly+normalization+two-MLP anchored forward time.

    No label/candidate-Full argument exists. Caller separately counts cold model
    load, anchor acquisition and any later 8-RHS Full fallback.
 """
 import torch
 started=time.perf_counter()
 with torch.no_grad():prediction=model.prediction(h,context['x0'],context['f0'],context['r0'],context['J0'])
 return prediction,time.perf_counter()-started
