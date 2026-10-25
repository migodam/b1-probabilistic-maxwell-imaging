"""Deterministic pushforward of the original C2 decoder under any Gaussian law.
Population marginal variance (ddof0 limit); no posterior correctness assertion.
"""
import numpy as np
from scipy.special import ndtr, owens_t
from numpy.polynomial.legendre import leggauss
from numpy.polynomial import Polynomial
from functools import lru_cache

@lru_cache(None)
def nodes(q):return leggauss(q)
def H(s):
 t=np.clip((np.asarray(s)+1)/2,0,1);return t*t*t*(10+t*(-15+6*t))
def phi(z):return np.exp(-.5*np.asarray(z)**2)/np.sqrt(2*np.pi)
def stable_cdf_difference(a,b):
 a,b=np.broadcast_arrays(a,b);return np.where(a>=0,ndtr(-a)-ndtr(-b),ndtr(b)-ndtr(a))
def bvn_same(t,rho):
 return ndtr(t)-2*owens_t(t,np.sqrt(np.maximum(0,(1-rho)/(1+rho))))
def bvn(t,u,rho,q=64):
 t,u,rho=np.broadcast_arrays(t,u,rho);g,w=nodes(q);angle=np.arcsin(np.clip(rho,-1,1));theta=angle[...,None]*(g+1)/2;r=np.sin(theta)
 # rho=sin(theta) removes the 1/sqrt(1-rho^2) endpoint singularity.
 exponent=-(t[...,None]**2-2*r*t[...,None]*u[...,None]+u[...,None]**2)/(2*np.cos(theta)**2)
 correction=angle/(4*np.pi)*np.sum(np.exp(np.minimum(exponent,0))*w,-1)
 result=ndtr(t)*ndtr(u)+correction
 # Near-singular conditional properties use adaptive integration, preserving
 # arbitrary correlation instead of trusting a fixed endpoint quadrature.
 near=np.abs(rho)>.95
 if np.any(near):
  from scipy.integrate import quad
  result=np.array(result,copy=True);flat=result.reshape(-1);tt=t.reshape(-1);uu=u.reshape(-1);aa=angle.reshape(-1)
  for i in np.flatnonzero(near.reshape(-1)):
   f=lambda th:np.exp(-(tt[i]**2-2*np.sin(th)*tt[i]*uu[i]+uu[i]**2)/(2*np.cos(th)**2))/(2*np.pi)
   flat[i]=ndtr(tt[i])*ndtr(uu[i])+quad(f,0,aa[i],epsabs=2e-13,epsrel=2e-13,limit=300)[0]
 return np.clip(result,0,1)
def phase_moments(ms,vs,q=64):
 """Seam-split integration, exact plateau tails; stable for extreme/degenerate Gaussians."""
 ms,vs=np.broadcast_arrays(np.asarray(ms,float),np.asarray(vs,float));sd=np.sqrt(np.maximum(vs,1e-300))
 a=(-1-ms)/sd;b=(1-ms)/sd;lo=np.maximum(a,-12);hi=np.minimum(b,12);width=np.maximum(hi-lo,0)
 g,w=nodes(q);z=lo[...,None]+width[...,None]*(g+1)/2;s=ms[...,None]+sd[...,None]*z
 weight=width[...,None]/2*w*phi(z);h=H(s)
 tail=ndtr((ms-1)/sd);e=tail+np.sum(weight*h,-1);e2=tail+np.sum(weight*h*h,-1)
 deg=vs<=1e-24;return np.where(deg,H(ms),e),np.where(deg,H(ms)**2,e2)
def truncated_phase_moments(ms,vs):
 """Polynomial/truncated-normal recurrence, fallback on dangerous cancellation."""
 ms,vs=np.broadcast_arrays(np.asarray(ms,float),np.asarray(vs,float));out=[]
 for m,v in zip(ms.ravel(),vs.ravel()):
  if v<=1e-24:out.append((float(H(m)),float(H(m)**2)));continue
  sd=np.sqrt(v);a=(-1-m)/sd;b=(1-m)/sd
  if sd<.1 or sd>3 or abs(a)>8 or abs(b)>8:out.append(tuple(float(z) for z in phase_moments(m,v)));continue
  I=[float(stable_cdf_difference(a,b)),float(phi(a)-phi(b))]
  for n in range(2,11):I.append(a**(n-1)*phi(a)-b**(n-1)*phi(b)+(n-1)*I[n-2])
  s=Polynomial([m,sd]);p=Polynomial([.5,.9375,0,-.625,0,.1875])(s)
  vals=[]
  for poly in [p,p*p]:vals.append(ndtr(-b)+sum(c*I[k] for k,c in enumerate(poly.coef)))
  if min(vals)<-1e-10 or max(vals)>1+1e-10:vals=[float(z) for z in phase_moments(m,v)]
  out.append(vals)
 return tuple(np.asarray(out)[:,i].reshape(ms.shape) for i in (0,1))
def correlated_cell(ms,vs,mu,Cuu,c,intervals,q=96):
 """Two probit properties jointly coupled to field. Returns means and 2x2 covariance."""
 mu=np.asarray(mu,float);Cuu=np.asarray(Cuu,float);c=np.asarray(c,float);intervals=np.asarray(intervals,float)
 if vs<0 or np.linalg.eigvalsh(np.block([[np.array([[vs]]),c[None]],[c[:,None],Cuu]])).min()<-1e-10:raise ValueError('nonPSD Gaussian')
 if vs<=1e-24:
  z=np.array([0.]);weight=np.array([1.]);s=np.array([ms]);cm=mu[None];V=Cuu
 else:
  sd=np.sqrt(vs);edges=sorted(set([-12.,12.,float(np.clip((-1-ms)/sd,-12,12)),float(np.clip((1-ms)/sd,-12,12))]));g,w=nodes(q)
  z=np.concatenate([a+(b-a)*(g+1)/2 for a,b in zip(edges[:-1],edges[1:])]);weight=np.concatenate([(b-a)/2*w for a,b in zip(edges[:-1],edges[1:])])*phi(z);s=ms+sd*z;cm=mu+z[:,None]*c/sd;V=Cuu-np.outer(c,c)/vs
 V=(V+V.T)/2;vv=np.maximum(np.diag(V),0);t=cm/np.sqrt(1+vv);p=ndtr(t);p2=bvn_same(t,vv/(1+vv));d0=intervals[:,0];d1=intervals[:,1]-d0;ed=d0+d1*p;ed2=d0*d0+2*d0*d1*p+d1*d1*p2
 h=H(s);mean=np.sum(weight[:,None]*h[:,None]*ed,axis=0);sec=np.sum(weight[:,None]*h[:,None]**2*ed2,axis=0)
 rho=V[0,1]/np.sqrt((1+vv[0])*(1+vv[1]));crossp=bvn(t[:,0],t[:,1],np.full(len(t),rho));crossd=d0[0]*d0[1]+d1[0]*d0[1]*p[:,0]+d0[0]*d1[1]*p[:,1]+d1[0]*d1[1]*crossp
 cross=np.sum(weight*h*h*crossd);cov=np.array([[sec[0],cross],[cross,sec[1]]])-np.outer(mean,mean)
 return mean,cov

def raw_gaussian(model,y):
 import torch
 with torch.inference_mode():
  t=torch.from_numpy(model._y(y)).to(model.device);mean,ls=model.model.head(model.model.encoder(t)).chunk(2,-1)
  m=mean.cpu().numpy().astype(float);v=np.exp(2*ls.clamp(-5,3).cpu().numpy().astype(float))
 st=model.stats;F=np.asarray(st.get('frame_inverse',st['frame'].T));m=(m*st['x_scale']+st['x_mean'])@F;factors=np.sqrt(v)*st['x_scale'];C=np.einsum('bi,ij,ik->bjk',factors**2,F,F)
 return m,C

def diagonal_image(m,C,decoder,q=64):
 cfg=decoder['config'];A=np.asarray(decoder['A'])/cfg['phase_transition_halfwidth'];ms=cfg['field_mean']/cfg['phase_transition_halfwidth']+m[:,:64]@A.T
 block=C[:,:64,:64];diag=np.diagonal(block,axis1=1,axis2=2);off=block.copy();idx=np.arange(64);off[:,idx,idx]=0
 if np.count_nonzero(off)==0:vs=diag@(A*A).T
 else:vs=np.einsum('pi,bij,pj->bp',A,block,A,optimize=True)
 cross=np.einsum('pi,bij->bpj',A,C[:,:64,64:]);maxcross=float(np.max(abs(cross)))
 if maxcross>1e-12:raise ValueError('correlated production law requires correlated_cell; diagonal product prohibited')
 e,e2=phase_moments(ms,vs,q);u=m[:,64:];v=np.diagonal(C[:,64:,64:],axis1=1,axis2=2);t=u/np.sqrt(1+v);p=ndtr(t);p2=bvn_same(t,v/(1+v));ints=np.array([cfg['inclusion_increment_real_interval'],cfg['inclusion_increment_imag_interval']]);d0=ints[:,0];d1=ints[:,1]-d0;d=d0+d1*p;d2=d0*d0+2*d0*d1*p+d1*d1*p2
 mean=d[:,:,None]*e[:,None];var=np.maximum(0,d2[:,:,None]*e2[:,None]-mean**2);mask=np.zeros(576);mask[decoder['roi']]=1;mean*=mask;var*=mask;mean+=np.array(cfg['host_chi'])[None,:,None]
 return np.concatenate((mean,var),axis=1).reshape(-1,4,24,24).astype('float32'),maxcross

def gaussian_image(m,C,decoder,q=96,return_joint=False):
 """Full Gaussian original-x pushforward, retaining field/property couplings.
 Four R4 marginal channels plus optional pointwise ReIm crosscovariance.
 This is a pointwise moment representation, never a joint spatial sampler.
 """
 m=np.atleast_2d(np.asarray(m,float));C=np.asarray(C,float)
 if C.ndim==2:C=C[None]
 if C.shape!=(len(m),66,66) or m.shape[1]!=66:raise ValueError('Gaussianoriginalx66 dimensions')
 if not np.isfinite(m).all() or not np.isfinite(C).all() or np.min(np.linalg.eigvalsh(C))<-1e-10:raise ValueError('finitePSDGaussian')
 cfg=decoder['config'];A=np.asarray(decoder['A'])/cfg['phase_transition_halfwidth'];roi=np.asarray(decoder['roi']);ints=np.array([cfg['inclusion_increment_real_interval'],cfg['inclusion_increment_imag_interval']]);cross=np.einsum('pi,bij->bpj',A,C[:,:64,64:]);out=np.zeros((len(m),4,576));out[:,:2]=np.asarray(cfg['host_chi'])[None,:,None];joint=np.zeros((len(m),576))
 if np.max(abs(cross))<1e-12 and not return_joint:return diagonal_image(m,C,decoder,q)[0]
 for b in range(len(m)):
  ms=cfg['field_mean']/cfg['phase_transition_halfwidth']+A@m[b,:64];vs=np.einsum('pi,ij,pj->p',A,C[b,:64,:64],A)
  for p in roi:
   mean,cov=correlated_cell(ms[p],vs[p],m[b,64:],C[b,64:,64:],cross[b,p],ints,q)
   out[b,:2,p]+=mean;out[b,2:,p]=np.maximum(np.diag(cov),0);joint[b,p]=cov[0,1]
 return (out.reshape(-1,4,24,24).astype('float32'),joint.reshape(-1,24,24)) if return_joint else out.reshape(-1,4,24,24).astype('float32')
