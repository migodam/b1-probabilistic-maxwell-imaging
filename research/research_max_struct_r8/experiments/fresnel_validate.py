from pathlib import Path
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
import sys,json,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT.parents[1]))
import numpy as np
from scipy.special import hankel1,h1vp,jv,jvp
from scipy.linalg import lstsq
from research.research_max_struct_r8.src.fresnel_adapter import prepare
from research.research_fbg_r3_experiments.src.physics import MaxwellModel


def reference(d,cs,a=.015,eps=3.,order=8):
 """Independent TM interface collocation with outgoing cylindrical multipoles.

 Both E_z and normal derivative continuous for nonmagnetic dielectrics.
 This reference has no grid/pixel permittivity, no material inverse tuning.
 """
 k=2*np.pi*2e9/299792458.;ki=k*np.sqrt(eps);modes=np.arange(-order,order+1);nm=len(modes);nc=len(cs)
 ang=2*np.pi*np.arange(64)/64;norm=np.c_[np.cos(ang),np.sin(ang)]
 points=np.concatenate([c+a*norm for c in cs]);normal=np.tile(norm,(nc,1));N=len(points)
 B=np.zeros((2*N,2*nc*nm),complex)
 def wave(c,interior=False):
  v=points-c;dist=np.linalg.norm(v,axis=1);unit=v/dist[:,None];theta=np.arctan2(v[:,1],v[:,0]);phase=np.exp(1j*theta[:,None]*modes)
  fun,der,kw=(jv,jvp,ki) if interior else (hankel1,h1vp,k)
  W=fun(modes[None,:],kw*dist[:,None]);dr=np.sum(unit*normal,axis=1);dt=np.sum(np.c_[-unit[:,1],unit[:,0]]*normal,axis=1)/dist
  return W*phase,phase*(kw*der(modes[None,:],kw*dist[:,None])*dr[:,None]+1j*modes*W*dt[:,None])
 for i,c in enumerate(cs):
  W,D=wave(c);B[:N,i*nm:(i+1)*nm]=W;B[N:,i*nm:(i+1)*nm]=D/k
  W,D=wave(c,True);mask=np.repeat(np.arange(nc)==i,64);cols=slice((nc+i)*nm,(nc+i+1)*nm)
  B[:N,cols]=-W*mask[:,None];B[N:,cols]=-D/k*mask[:,None]
 inc=d['incident_evaluator'](points);der=d['incident_evaluator'](points,normal)
 rhs=-np.r_[inc,der/k];scale=np.linalg.norm(B,axis=0);coef,_,rank,_=lstsq(B/scale,rhs,cond=1e-12);coef/=scale[:,None]
 ra=np.arange(72)*np.pi/36;rx=.76*np.c_[np.cos(ra),np.sin(ra)];pred=np.zeros((72,36),complex)
 for i,c in enumerate(cs):
  v=rx-c;dist=np.linalg.norm(v,axis=1);theta=np.arctan2(v[:,1],v[:,0]);W=hankel1(modes[None,:],k*dist[:,None])*np.exp(1j*theta[:,None]*modes)
  pred+=W@coef[i*nm:(i+1)*nm]
 return pred,dict(interface_relative_residual=float(np.linalg.norm(B@coef-rhs)/np.linalg.norm(rhs)),rank=int(rank),order=order)


def main():
 cpu=time.process_time();wall=time.perf_counter();rows=[]
 for name,cs in [('dielTM_dec8f',[np.array([-.03,0.])]),('twodielTM_8f',[np.array([-.045,0.]),np.array([.045,0.])])]:
  d=prepare(ROOT/'data_local/fresnel'/f'{name}.exp');p,rec=reference(d,cs);mask=np.isfinite(d['Y']);train=d['train_mask'];hold=d['withheld_mask']
  rel=lambda p,mask:float(np.linalg.norm((p-d['Y'])[mask])/np.linalg.norm(d['Y'][mask]))
  rec.update(case=name,reference='Known nominal circular geometry only; not exact material GT',all_measured_relative=rel(p,mask),train_relative=rel(p,train),heldout_relative=rel(p,hold))
  # A free complex gain is only diagnostic; it is never used by inversions.
  alpha=np.vdot(p[train],d['Y'][train])/np.vdot(p[train],p[train]);rec['diagnostic_gain_real']=float(alpha.real);rec['diagnostic_gain_imag']=float(alpha.imag);rec['gain_diagnostic_heldout_relative']=rel(alpha*p,hold)
  rec['discretization']=[]
  for side in [24,36]:
   z=prepare(ROOT/'data_local/fresnel'/f'{name}.exp',side=side);h=.16/side;ctr=z['centers'];frac=np.zeros(len(ctr));sub=(np.arange(12)+.5)/12-.5
   for dx in sub:
    for dy in sub:
     pts=ctr+h*np.array([dx,dy]);inside=np.zeros(len(ctr),bool)
     for c in cs:inside |= np.linalg.norm(pts-c,axis=1)<=.015
     frac+=inside/144
   mm=2*frac;model=MaxwellModel(z['public']);pp=model.forward(mm)
   rec['discretization'].append(dict(side=side,relative_to_independent_reference=float(np.linalg.norm((pp-p)[mask])/np.linalg.norm(p[mask])),relative_to_measured=rel(pp,mask),contrast_integral=float(mm.sum()*h*h)))
  # Verify receiver adjoint and real-material derivative, meaningful for mask adapter.
  model=d['model'];rng=np.random.default_rng(806);J=rng.normal(size=model.E.shape)+1j*rng.normal(size=model.E.shape);R=rng.normal(size=d['Y_train'].shape)+1j*rng.normal(size=d['Y_train'].shape)
  lhs=np.vdot(model.receive(J),R);rhs=np.vdot(J,model.receive_adjoint(R));rec['receiver_adjoint_relative']=float(abs(lhs-rhs)/max(abs(lhs),abs(rhs)))
  m=np.full(model.n,.1);direction=rng.normal(size=model.n);direction/=np.linalg.norm(direction);v,g=model.data_value_gradient(m,d['Y_train']);step=1e-5
  fd=(model.physical_score(m+step*direction,d['Y_train'])-model.physical_score(m-step*direction,d['Y_train']))/(2*step);rec['gradient_directional_analytic']=float(g@direction);rec['gradient_directional_finite_difference']=float(fd)
  rows.append(rec);print(json.dumps(rec),flush=True)
 (ROOT/'RESULTS/fresnel/REFERENCE_VALIDATION.json').write_text(json.dumps(dict(cases=rows,cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall),indent=2)+'\n')
if __name__=='__main__':main()
