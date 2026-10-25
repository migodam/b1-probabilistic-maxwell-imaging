import json,time,resource,sys
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.special import ndtr
from material_moments import *
OUT=Path(__file__).resolve().parent
start=time.monotonic();rng=np.random.default_rng(2026102604);records=[]
for ms,vs in [(0.,1.),(.9,.03),(-.9,.03),(4.,.01),(-4.,.01),(0.,0.),(1.,1e-20),(-1.,1e-20),(8.,100.),(0.,.0001),(-20.,1.)]:
 e,e2=phase_moments(ms,vs,96);rec=truncated_phase_moments(ms,vs)
 ref=[]
 if vs==0:ref=[H(ms),H(ms)**2]
 else:
  sd=np.sqrt(vs);a=(-1-ms)/sd;b=(1-ms)/sd
  for power in [1,2]:
   lo=max(-14,a);hi=min(14,b);v=quad(lambda z:H(ms+sd*z)**power*phi(z),lo,hi,epsabs=2e-13,epsrel=2e-13,limit=300)[0] if hi>lo else 0.
   ref.append(v+ndtr(-b))
 records.append(dict(test='phase',ms=ms,vs=vs,quadrature=[float(e),float(e2)],reference=list(map(float,ref)),quadrature_error=float(np.max(abs(np.array([e,e2])-ref))),recurrence_error=float(np.max(abs(np.array(rec)-ref)))))
intervals=np.array([[.2,1.6],[.01,.4]])
for name,L,mean in [('correlated',np.array([[.9,0,0],[.5,.7,0],[-.2,.25,.6]]),[.2,-.5,.3]),('strong_correlation',np.array([[1,0,0],[.8,.05,0],[-.7,.02,.06]]),[-.8,.4,-.2]),('zero_field',np.diag([0,.4,.5]),[.2,-3,3]),('tail',np.diag([.1,.3,.5]),[-4,0,0]),('near_saturated',np.diag([.04,.3,.5]),[1.02,0,0])]:
 C=L@L.T;mu=np.array(mean);det,cov=correlated_cell(mu[0],C[0,0],mu[1:],C[1:,1:],C[0,1:],intervals,96)
 # Independent adaptive quadrature of conditional property and conditional field.
 v=C[0,0];sd=np.sqrt(v)
 def pref(z,power,j):
  if v==0:s=mu[0];cm=mu[1+j];vv=C[1+j,1+j]
  else:s=mu[0]+sd*z;cm=mu[1+j]+C[0,1+j]/sd*z;vv=max(0,C[1+j,1+j]-C[0,1+j]**2/v)
  d0,dmax=intervals[j];d1=dmax-d0
  def f(u):return (d0+d1*ndtr(cm+np.sqrt(vv)*u))**power*phi(u)
  return H(s)**power*quad(f,-12,12,epsabs=1e-11,epsrel=1e-11)[0]*(phi(z) if v else 1.)
 ref=np.zeros(2);sec=np.zeros(2)
 for j in range(2):
  if v:
   cuts=sorted(set([-12.,12.,float(np.clip((-1-mu[0])/sd,-12,12)),float(np.clip((1-mu[0])/sd,-12,12))]))
   for power,dest in [(1,ref),(2,sec)]:dest[j]=sum(quad(lambda z:pref(z,power,j),a,b,epsabs=3e-11,epsrel=3e-11)[0] for a,b in zip(cuts[:-1],cuts[1:]))
  else:
   ref[j]=pref(0,1,j);sec[j]=pref(0,2,j)
 n=1000000;sm=rng.normal(size=(n,3))@L.T+mu;h=H(sm[:,0]);d=intervals[:,0]+(intervals[:,1]-intervals[:,0])*ndtr(sm[:,1:]);mat=d*h[:,None];mc=mat.mean(0);mcov=np.cov(mat,rowvar=False,ddof=0);se=mat.std(0)/np.sqrt(n)
 prod=mat[:,:,None]*mat[:,None,:];mc_cov_se=prod.std(0)/np.sqrt(n)+abs(det[:,None])*se[None,:]+abs(det[None,:])*se[:,None]
 records.append(dict(MC_covariance_max_z=float(np.max(abs(mcov-cov)/np.maximum(mc_cov_se,1e-16))),test=name,deterministic_mean=det.tolist(),deterministic_covariance=cov.tolist(),reference_mean=ref.tolist(),reference_variance=(sec-ref**2).tolist(),reference_abs_error=float(max(np.max(abs(det-ref)),np.max(abs(np.diag(cov)-(sec-ref**2))))),MC_n=n,MC_mean=mc.tolist(),MC_covariance=mcov.tolist(),MC_mean_max_z=float(np.max(abs(mc-det)/np.maximum(se,1e-16)))))
res=dict(status='PASS' if max(r.get('quadrature_error',0) for r in records)<1e-9 and max(r.get('recurrence_error',0) for r in records)<1e-7 and max(r.get('reference_abs_error',0) for r in records)<1e-8 and max(r.get('MC_mean_max_z',0) for r in records)<6 and max(r.get('MC_covariance_max_z',0) for r in records)<6 else 'FAIL',records=records,wall_seconds=time.monotonic()-start,physics_RHS=0)
(OUT/'MOMENT_VALIDATION.json').write_text(json.dumps(res,indent=2)+'\n');print(res['status'],res['wall_seconds']);print([(r['test'],r.get('quadrature_error'),r.get('reference_abs_error'),r.get('MC_mean_max_z')) for r in records])
