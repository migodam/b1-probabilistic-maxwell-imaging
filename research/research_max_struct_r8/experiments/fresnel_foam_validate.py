from fresnel_validate import *
from research.research_max_struct_r8.src.paper_physics import MatrixFreeModel
from scipy.special import h1vp

def nominal(coords):
 # Figure-based nominal positions; no measured-target pose/material fitting.
 # Foam centre0; innerplastic−5mm; externalplastic−60mm. Both axes/locations
 # are nominal diagnostics and not pixelGT (Geffrin2005; Dubois2005 Fig7).
 ctr=np.asarray(coords);foam=np.linalg.norm(ctr,axis=1)<.04;inside=np.linalg.norm(ctr-np.array([-.005,0]),axis=1)<.0155;ext=np.linalg.norm(ctr-np.array([-.060,0]),axis=1)<.0155
 return np.where(inside|ext,2.,np.where(foam,.45,0.))

def reference_foam(d,order=18):
 k=2*np.pi*8e9/299792458.;kf=k*np.sqrt(1.45);kd=k*np.sqrt(3.);nm=2*order+1;modes=np.arange(-order,order+1);cs=[np.array([0.,0.]),np.array([-.060,0.]),np.array([-.005,0.])];rr=[.04,.0155,.0155]
 theta=2*np.pi*np.arange(80)/80;unit=np.c_[np.cos(theta),np.sin(theta)];pts=np.concatenate([c+r*unit for c,r in zip(cs,rr)]);norm=np.tile(unit,(3,1));N=len(pts);B=np.zeros((2*N,6*nm),complex)
 def wave(c,kw,kind):
  v=pts-c;dist=np.linalg.norm(v,axis=1);u=v/dist[:,None];ang=np.arctan2(v[:,1],v[:,0]);phase=np.exp(1j*ang[:,None]*modes);f,df=(hankel1,h1vp) if kind=='H' else (jv,jvp);W=f(modes[None,:],kw*dist[:,None]);dr=np.sum(u*norm,axis=1);dt=np.sum(np.c_[-u[:,1],u[:,0]]*norm,axis=1)/dist
  return W*phase,phase*(kw*df(modes[None,:],kw*dist[:,None])*dr[:,None]+1j*modes*W*dt[:,None])
 # afoamHk,aextHk,bfoamJkf,cinnerHkf,dinnerJkd,eextJkd
 specs=[(cs[0],k,'H',[1,1,0]),(cs[1],k,'H',[1,1,0]),(cs[0],kf,'J',[-1,0,1]),(cs[2],kf,'H',[-1,0,1]),(cs[2],kd,'J',[0,0,-1]),(cs[1],kd,'J',[0,-1,0])]
 for i,(c,kw,kind,signs) in enumerate(specs):
  W,D=wave(c,kw,kind);sg=np.repeat(signs,80);B[:N,i*nm:(i+1)*nm]=sg[:,None]*W;B[N:,i*nm:(i+1)*nm]=sg[:,None]*D/k
 exterior=np.repeat([1,1,0],80);inc=d['incident_evaluator'](pts)*exterior[:,None];der=d['incident_evaluator'](pts,norm)*exterior[:,None];rhs=-np.r_[inc,der/k];scale=np.linalg.norm(B,axis=0);z,_,rank,_=lstsq(B/scale,rhs,cond=1e-12);z/=scale[:,None]
 a=2*np.pi*np.arange(360)/360;rx=1.67*np.c_[np.cos(a),np.sin(a)];pred=np.zeros((360,18),complex)
 for i,c in enumerate(cs[:2]):
  v=rx-c;dist=np.linalg.norm(v,axis=1);ang=np.arctan2(v[:,1],v[:,0]);W=hankel1(modes[None,:],k*dist[:,None])*np.exp(1j*ang[:,None]*modes);pred+=W@z[i*nm:(i+1)*nm]
 return pred,dict(order=order,rank=int(rank),interface_relative_residual=float(np.linalg.norm(B@z-rhs)/np.linalg.norm(rhs)))

def main():
 cpu=time.process_time();wall=time.perf_counter();d=prepare(ROOT/'data_local/fresnel/FoamTwinDielTM.exp',frequency=8.,side=40,width=.2,matrix_free=True);pred,rec=reference_foam(d);mask=np.isfinite(d['Y']);train=d['train_mask'];hold=d['withheld_mask'];gain=np.vdot(pred[train],d['Y'][train])/np.vdot(pred[train],pred[train])
 rec.update(qualification='Nominal literaturegeometrydiagnostic, coordinatesapproxfromprimaryfigures, noexactpixelGT, nooptimizationofpose/material',
  nominal_parameters=dict(foam_center=[0,0],foam_radius=.04,foam_epsilon=1.45,inner_center=[-.005,0],outer_center=[-.060,0],plastic_radius=.0155,plastic_epsilon=3.),
  train_relative=float(np.linalg.norm((pred-d['Y'])[train])/np.linalg.norm(d['Y'][train])),heldout_relative=float(np.linalg.norm((pred-d['Y'])[hold])/np.linalg.norm(d['Y'][hold])),
  diagnostic_train_gain=[float(gain.real),float(gain.imag)],gain_diagnostic_heldout_relative=float(np.linalg.norm((gain*pred-d['Y'])[hold])/np.linalg.norm(d['Y'][hold])),training_channels=int(train.sum()),withheld_channels=int(hold.sum()),discretization=[])
 for side in [40,60]:
  z=d if side==40 else prepare(ROOT/'data_local/fresnel/FoamTwinDielTM.exp',frequency=8.,side=side,width=.2,matrix_free=True);h=.2/side;mm=np.zeros(side*side);sub=(np.arange(8)+.5)/8-.5
  for dx in sub:
   for dy in sub:mm+=nominal(z['centers']+h*np.array([dx,dy]))/64
  model=MatrixFreeModel(z['public']);st=model.state(mm);p=st['Y'];rec['discretization'].append(dict(side=side,relative_to_independent_reference=float(np.linalg.norm((p-pred)[mask])/np.linalg.norm(pred[mask])),relative_to_measurements=float(np.linalg.norm((p-z['Y'])[mask])/np.linalg.norm(z['Y'][mask])),state_relative_residual=st['state_relative_residual'],physical_cost=dict(model.cost)))
 model=d['model'];rng=np.random.default_rng(806);J=rng.normal(size=model.E.shape)+1j*rng.normal(size=model.E.shape);R=rng.normal(size=d['Y_train'].shape)+1j*rng.normal(size=d['Y_train'].shape);lhs=np.vdot(model.receive(J),R);rhs=np.vdot(J,model.receive_adjoint(R));rec['receiver_adjoint_relative']=float(abs(lhs-rhs)/max(abs(lhs),abs(rhs)))
 m=np.full(model.n,.1);direction=rng.normal(size=model.n);direction/=np.linalg.norm(direction);v,g=model.data_value_gradient(m,d['Y_train']);step=1e-4;fd=(model.physical_score(m+step*direction,d['Y_train'])-model.physical_score(m-step*direction,d['Y_train']))/(2*step);rec['gradient_directional_analytic']=float(g@direction);rec['gradient_directional_finite_difference']=float(fd)
 rec.update(cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,incident_calibration=d['calibration']);dump(ROOT/'RESULTS/fresnel/FOAM_REFERENCE_VALIDATION.json',rec);print(json.dumps({k:v for k,v in rec.items() if k!='incident_calibration'}),flush=True)
from fresnel_invert import dump
if __name__=='__main__':main()
