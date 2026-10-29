"""R9 DEV stabilized weak subspaces, frozen R7 target. E0/E1 no physics.
Only this module owns exploration ledger and outputs; history is read-only.
"""
from __future__ import annotations
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'):os.environ[k]='1'
import argparse, importlib.util, json, sys, time, uuid, fcntl, resource
from pathlib import Path
from collections import Counter
from contextlib import contextmanager
import numpy as np
from scipy.linalg import null_space,eigh
from scipy.stats import chi2
P=Path(__file__).resolve().parents[1]
ROOT=P.parents[3]
H7=ROOT/'B1/b1-p3l-r7-material-fibers/research/p3l_r7_material_solution_fibers'
# Resolve from B1 rather than relying on callers' working directory.
H7=P.parents[2]/'b1-p3l-r7-material-fibers/research/p3l_r7_material_solution_fibers'
sys.path.insert(0,str(H7/'src'))
spec=importlib.util.spec_from_file_location('r9_frozen_r7_physics',H7/'src/physics.py')
hist=importlib.util.module_from_spec(spec);sys.modules[spec.name]=hist;spec.loader.exec_module(hist)
from material_jets import real_material, material_derivative, decode, match_distance, feasible_progress,gn_step
OUT=P/'results/exploration';AUD=P/'audit/exploration';CONF=P/'config/exploration_freeze.json'
ARMS=['A_NATIVE_FULLJ','B_FIXED_WEAK','C_DYNAMIC_WEAK','S2_GLOBAL','D_FOURIER','D_WAVELET','E_SOM','F_FROZENJ','R_RANDOM','U_UNCORRECTED']
CAP=24000
IO=Counter();_rawsave=hist.save
def tracked_save(path,**data):
 st=time.perf_counter();IO['array_save_attempts']+=1
 try:
  _rawsave(path,**data);IO['array_save_completed']+=1;IO['array_save_bytes']+=Path(path).stat().st_size
 except BaseException:
  IO['array_save_failures']+=1;raise
 finally:IO['array_save_seconds']+=time.perf_counter()-st
hist.save=tracked_save

def write(p,x):hist.write(p,x)
def save(p,**x):hist.save(p,**x)
def append(p,x):
 p.parent.mkdir(parents=True,exist_ok=True)
 with p.open('a') as f:f.write(json.dumps(hist.plain(x),allow_nan=False,separators=(',',':'))+'\n');f.flush();os.fsync(f.fileno())
def arrays(p):
 with np.load(p,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}
def point(p):
 a=arrays(p);return hist.Point(a['x'],a['q'],hist.EndpointState(a['material'],a['current'],a['field'],a['prediction'],float(a['residual'])),str(p),a.get('receiver_T'),a.get('J'),a.get('Jnative'))
def prior_only():
 # OriginalTarget construction does not load the physics backend.
 class Dummy:finished=False;test_mode=False
 return hist.Target(Dummy()).prior

def frame(prior,x,J):
 t=time.perf_counter();M=real_material(material_derivative(prior,x));U,s,Vh=np.linalg.svd(M,full_matrices=False)
 keep=s>max(s[0]*1e-8,1e-10);T=Vh[keep].T/s[keep][None,:];K=J@T
 ev,W=eigh(K.T@K+T.T@T);W=W[:,:min(4,len(ev))]
 Qw=M@T@W;Qweak,_=np.linalg.qr(T@W);v=T@W[:,0]
 if v[np.argmax(abs(v))]<0:v=-v
 return dict(M=M,T=T,K=K,Qw=Qw,Qweak=Qweak,v=v,rank=int(keep.sum()),material_singular=s,prior_aware_eigenvalues=ev,svd_seconds=time.perf_counter()-t)

def alternate(prior,base,g,kind):
 # Material bases are mapped by least squares back into original x66 support.
 pts=prior.points;xx=pts[:,0]+.5;yy=pts[:,1]+.5;cols=[]
 if kind=='FOURIER':
  for k,l in [(0,0),(1,0),(0,1),(1,1),(2,0),(0,2),(2,1),(1,2)]:
   cols += [np.cos(2*np.pi*(k*xx+l*yy)),np.sin(2*np.pi*(k*xx+l*yy)) if k+l else np.ones(576)]
 else:
  for iy in range(4):
   for ix in range(4):cols.append(((xx>=ix/4)&(xx<(ix+1)/4)&(yy>=iy/4)&(yy<(iy+1)/4)).astype(float))
 B=np.column_stack(cols)*prior.mask[:,None]
 inc=decode(prior,base.x)[1]['inclusion_increment'];B=real_material(inc*B)
 L=np.linalg.lstsq(g['M'],B,rcond=1e-8)[0];projected=g['M']@L
 err=float(np.linalg.norm(projected-B)/np.linalg.norm(B))
 _,s,Vh=np.linalg.svd(projected,full_matrices=False);keep=s>max(s[0]*1e-8,1e-10)
 lift=L@Vh[keep].T/s[keep][None,:]
 ev,W=eigh((base.J@lift).T@(base.J@lift)+lift.T@lift);v=lift@W[:,0]
 if v[np.argmax(abs(v))]<0:v=-v
 return v,dict(projection_relative_error=err,mapped_rank=int(keep.sum()),target_basis_columns=B.shape[1],same_original_latent_support=True,basis_restriction=True)

def som(prior,base,g):
 bg=arrays(hist.R4/'research/p3l_r4_conditional_transport/inputs/background/weak.npz');S,G=bg['S'],bg['G']
 _,sv,Vh=np.linalg.svd(S,full_matrices=False);V=Vh.conj().T
 # Receiver-ambiguous current projected outside observed rank32 row space.
 desired=base.state.current-V@(V.conj().T@base.state.current)
 rhs=desired-base.state.chi[:,None]*(G@desired)
 field=base.state.field;dchi=np.sum(np.conj(field)*rhs,axis=1)/np.maximum(np.sum(abs(field)**2,axis=1),1e-30)
 defect=float(np.linalg.norm(rhs-dchi[:,None]*field)/max(np.linalg.norm(rhs),1e-30))
 v=np.linalg.lstsq(g['M'],real_material(dchi),rcond=1e-8)[0]
 liftdef=float(np.linalg.norm(g['M']@v-real_material(dchi))/max(np.linalg.norm(real_material(dchi)),1e-30))
 info=dict(shared_material_equation_relative_defect=defect,decoder_tangent_lift_relative_defect=liftdef,receiver_singular_values=sv.tolist(),lift_tolerance=1e-6,all_illuminations_shared_material=True)
 if defect>1e-6 or liftdef>1e-6 or np.linalg.norm(g['M']@v)<1e-12:return None,dict(status='MATERIAL_LIFT_FAIL',scope='this chosen receiver-ambiguous current direction only; not all SOM directions',**info)
 v/=np.linalg.norm(g['M']@v)
 if v[np.argmax(abs(v))]<0:v=-v
 return v,dict(status='MATERIAL_LIFT_PASS',**info)

def projector_drift(Qa,Qb):
 c=np.linalg.svd(Qa.T@Qb,compute_uv=False);return float(np.sqrt(max(0.,1.-np.clip(c.min(),0.,1.)**2)))

def e0():
 OUT.mkdir(parents=True,exist_ok=True);AUD.mkdir(parents=True,exist_ok=True);prior=prior_only();rows=[]
 globalsum=np.zeros((66,66));gsources=[]
 for i in range(4):
  pp=H7/'E2/parents'/f'R6_CDEV_{i:04d}'/'ANCHOR.npz';bb=point(pp);qg=frame(prior,bb.x,bb.J)['Qweak'];globalsum+=qg@qg.T;gsources.append(str(pp))
 vals,vec=eigh(globalsum);Qglobal=vec[:,-4:];save(AUD/'S2_GLOBAL_FRAME.npz',Qweak=Qglobal,eigenvalues=vals);write(AUD/'S2_GLOBAL_PROVENANCE.json',dict(sources=gsources,role='original R7 development only; no locked objects',new_RHS=0))
 for idx in range(8):
  parent=f'R7_LOCKED_{idx:03d}';hp=H7/'E2/parents'/parent;b=point(hp/'ANCHOR.npz');g=frame(prior,b.x,b.J)
  # Saved analytic receiver chain identity uses arithmetic, no solver.
  C=np.einsum('rp,ps->rsp',b.T,b.state.field).reshape(256,576)
  d=C@material_derivative(prior,b.x)
  J=np.column_stack([hist.pack(d[:,i].reshape(32,8))*np.sqrt(2)/hist.SIGMA for i in range(66)])
  err=float(np.linalg.norm(J-b.J)/np.linalg.norm(b.J));de=float(np.linalg.norm(decode(prior,b.x)[0]-b.state.chi))
  geom=arrays(hp/'BASE_GEOMETRY.npz');samev=float(np.linalg.norm(g['v']-geom['v']))
  alts={};directions={}
  for kind in ['FOURIER','WAVELET']:
   directions['D_'+kind],alts[kind]=alternate(prior,b,g,kind)
  sv,si=som(prior,b,g)
  directions['E_SOM']=np.zeros(66) if sv is None else sv
  save(AUD/parent/'ONLINE_INPUT.npz',x=b.x,q=b.q,material=b.state.chi,prediction=b.state.data,current=b.state.current,field=b.state.field,residual=b.state.residual,J=b.J,Jnative=b.Jnative,receiver_T=b.T,observation=arrays(hp/'OBSERVATION.npz')['observation'],random_v=geom['random_v'],**directions)
  write(AUD/parent/'BASIS_METADATA.json',dict(D_FOURIER=alts['FOURIER'],D_WAVELET=alts['WAVELET'],E_SOM=si))
  rows.append(dict(parent=parent,role='EXPOSED_DEV_ONLY',anchor=str(hp/'ANCHOR.npz'),observation=str(hp/'OBSERVATION.npz'),geometry=str(hp/'BASE_GEOMETRY.npz'),decoder_error=de,J_relative_error=err,R7_direction_difference=samev,rank=g['rank'],alternates=alts,SOM=si,source_phi=float(np.vdot(b.state.data-arrays(hp/'OBSERVATION.npz')['observation'],b.state.data-arrays(hp/'OBSERVATION.npz')['observation']).real/hist.SIGMA**2)))
  save(AUD/parent/'ANCHOR_FRAME.npz',**g)
 # Exact decoder plateau gauge only: no physics or inferential gauge claim.
 x=np.zeros(66);x[0]=-100.;v=np.zeros(66);v[1]=.01
 dg=float(np.linalg.norm(decode(prior,x)[0]-decode(prior,x+v)[0]));d0=float(np.linalg.norm(material_derivative(prior,x)@v))
 write(AUD/'E0_SOURCE_AUDIT.json',dict(status='PASS' if max(max(r['decoder_error'],r['J_relative_error']) for r in rows)<1e-10 else 'FAIL',new_RHS=0,rows=rows,decoder_plateau_gauge=dict(finite_decoder_change=dg,tangent_material_norm=d0,physics_RHS=0,universal_gauge_claim=False),gates=dict(material_metric='norm(concat(Re dm,Im dm))/24',radius_relative_tolerance=.02,progress_fraction=.5,DeltaPhi_max=2,noise_ball_chi2_df=512,noise_ball_threshold=float(chi2.ppf(.95,512)/2),DeltaPrior_max=10,prior='0.5*x@x in original Gaussianized x66',legal='original prior legal q'),cost_model=dict(state_RHS=8,J_first_RHS=32,corrector_maxstates=3,corrected_arm_with_endpointJ_maxRHS=128,corrected_arms=8,objects=8,radii=2,endpoint_RHS_bound=16384,paths_BC_5nodes_8objects_RHS_bound=3200,uncorrected_RHS_bound=640,probe_RHS_bound=56,total_prospective_RHS=20280,cap=CAP),history_read_only=True,GT_opened=False))
 freeze=dict(version=1,campaign='R9_S_EXPLORATION',frozen_utc=time.time(),exposed_DEV_objects=[r['parent'] for r in rows],radii=[.01,.02],weak_rank=4,arms=ARMS,max_states=3,trust_initial=.5,trust_material_fraction=.5,max_latent_radius=4.,match_tolerance=.02,weak_definition='four smallest eigenvectors of K.T K + T.T T with M T orthonormal, effective material cutoff max(s0*1e-8,1e-10)',B='Qweak=QR(T Wweak) Gaussian-orthonormal; exclude fixed anchor Qweak from latent correction dx',C='exclude recomputed current Gaussian-orthonormal Qweak(x,J(x)), Procrustes aligned to anchor; same initial v as B',A='exact R7 native material-effective corrector',D='fixed 16-column Fourier or4x4Haar scaling material targets; mapped by M pseudoinverse into original x66; weakest prior-aware direction',S2_GLOBAL='original R7 four dev weak Gaussian-projector mean top4; fixed correction complement; deterministic top aggregate vector direction in global Qweak, anchor material-unit normalized',E='one SOM ambiguous current outside receiver row space; require shared all-source material and decoder tangent lift error<=1e-6 else MATERIAL_LIFT_FAIL',paths=['B_FIXED_WEAK','C_DYNAMIC_WEAK'],path_radius=.02,path_nodes=[0,.25,.5,.75,1],path_claim='five nodes sampled along actual accepted correction polyline by material arclength; no continuous certificate',retraction='A,F,R,D exact R7 radial; B,C,S2 decoder-only Newton constrained to respective Qs; A/B not pure basis causal comparison',C_scope='current weak-complement refresh only; same anchor initial direction as B; no multistep direction-coverage claim',J_accounting='algorithm J only when next correction needs current J; Frozen uses baseJ; selected endpoint J audit after algorithm clock+cost snapshot, paid restore8+J32 if necessary; U endpointJ audit separate',prior_and_gates='EXACT R7',GT_online=False,cap_RHS=CAP,costprobe_max_RHS=56,single_worker=True,live_LU_max=1,BLAS_threads=1,dtype='float64/complex128',physics_enable_required=True)
 if CONF.exists():
  old=json.loads(CONF.read_text());freeze['frozen_utc']=old['frozen_utc']
  if freeze!=old:raise RuntimeError('freeze differs; explicit versioning required')
 else:write(CONF,freeze)
 frozen=AUD/'FROZEN_EXPLORATION_SOURCE.py'
 if frozen.exists() and frozen.read_bytes()!=Path(__file__).read_bytes():raise RuntimeError('frozen source mismatch; explicit versioning required')
 if not frozen.exists():frozen.write_bytes(Path(__file__).read_bytes())
 print('E0_COMPLETE zeroRHS',flush=True)

class Book:
 def __init__(self,job,max_wall=7200):
  self.finished=False;self.test_mode=False;self.job=job;self.jobdir=OUT/'jobs'/job;self.root=OUT
  self.context={};self.counts=Counter();self.credit=0;self.active=None;self.max_wall=max_wall
  if not CONF.exists():raise RuntimeError('freeze missing')
  if (AUD/'FROZEN_EXPLORATION_SOURCE.py').read_bytes()!=Path(__file__).read_bytes():raise RuntimeError('source differs from prephysics freeze')
  auth=json.loads((P/'config/exploration_physics_enable.json').read_text())
  if auth.get('enabled') is not True:raise RuntimeError('physics disabled')
  self.cap=min(CAP,int(auth['cap_RHS']));self.lock=(P/'shared/physics_worker.lock').open('a');fcntl.flock(self.lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  if self.jobdir.exists():raise RuntimeError('no implicit job replay')
  self.jobdir.mkdir(parents=True);self.started=time.perf_counter();self.cpu=time.process_time();self.io0=IO.copy();self.total=self.totals()
  self._event('JOB_BEGIN',pid=os.getpid(),max_wall=max_wall)
 def totals(self):
  path=OUT/'RESOURCE_LEDGER.jsonl'
  if not path.exists():return dict(RHS=0,forward_RHS=0,derivative_RHS=0,states=0,LU=0,failures=0,sequence=0)
  raw=path.read_bytes()
  if not raw.endswith(b'\n'):raise RuntimeError('torn ledger')
  return json.loads(raw.splitlines()[-1])['after']
 def _event(self,event,rhs=0,states=0,lu=0,kind=None,**kw):
  if time.perf_counter()-self.started>self.max_wall and rhs:raise RuntimeError('controlled wall cap')
  t=self.total.copy()
  if t['RHS']+rhs>self.cap:raise RuntimeError('exploration RHS cap')
  for k,v in [('RHS',rhs),('states',states),('LU',lu)]:t[k]+=v
  if kind:t[kind+'_RHS']+=rhs
  if event=='PHYSICS_END' and kw.get('error') is not None:t['failures']+=1
  t['sequence']+=1
  append(OUT/'RESOURCE_LEDGER.jsonl',dict(event=event,job=self.job,utc=time.time(),RHS=rhs,states=states,LU=lu,kind=kind,context=self.context,**kw,after=t))
  self.total=t;write(OUT/'RESOURCE_SNAPSHOT.json',t)
 def reserve_state(self,material,*,resolution,purpose):
  if resolution!=24 or self.active is None or self.credit:raise RuntimeError('state contract')
  self._event('STATE_RESERVED',rhs=8,states=1,lu=1,kind='forward',purpose=purpose);self.credit=8
 def reserve_rhs(self,count,role):
  if self.active is None:raise RuntimeError('solve outside physics')
  if role=='full_forward':
   if self.credit!=count:raise RuntimeError('forward credit')
   self.credit-=count;self.counts['actual_forward_RHS']+=count
  else:self._event('DERIVATIVE_RESERVED',rhs=int(count),kind='derivative',role=role);self.counts['actual_derivative_RHS']+=count
 @contextmanager
 def span(self,label,**counts):
  st=time.perf_counter();self.counts.update(counts)
  try:yield
  finally:self.counts[label+'_seconds']+=time.perf_counter()-st
 @contextmanager
 def physical(self,label):
  if self.active is not None:raise RuntimeError('nested physics')
  self.active=uuid.uuid4().hex;st=time.perf_counter();err=None;self._event('PHYSICS_BEGIN',action=self.active,label=label)
  try:yield
  except BaseException as exc:err=repr(exc);raise
  finally:self._event('PHYSICS_END',action=self.active,label=label,error=err,wall_seconds=time.perf_counter()-st);self.active=None
 def close(self,error=None):
  if self.finished:return
  d=dict(status='FAILED' if error else 'COMPLETE',error=None if error is None else repr(error),wall_seconds=time.perf_counter()-self.started,CPU_seconds=time.process_time()-self.cpu,counters=dict(self.counts),unused_forward_credit=self.credit,peak_RSS=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,array_IO={k:IO[k]-self.io0[k] for k in set(IO)|set(self.io0)})
  self._event('JOB_END',**d);write(self.jobdir/'CLOSED.json',d);self.finished=True;self.lock.close()
 def __enter__(self):return self
 def __exit__(self,ty,err,tb):self.close(err)

class Target(hist.Target):
 def _persist(self,p,path=None):
  st=time.perf_counter();super()._persist(p,path);self.book.counts['endpoint_save_count']+=1;self.book.counts['endpoint_save_seconds']+=time.perf_counter()-st

# Metrics reproduce R7 code (no truth).
def metrics(target,p,base,y,unit,distance):
 dm=real_material(p.state.chi-base.state.chi);rms=float(np.linalg.norm(dm));phi=target.phi(p,y);phi0=target.phi(base,y);prior=.5*float(p.x@p.x);dp=prior-.5*float(base.x@base.x);pair=target.white(p.state.data-base.state.data)
 return dict(endpoint=p.key,material_RMS=rms,target_RMS=distance,distance_matched=abs(rms-distance)<=.02*distance,signed_progress=float(unit@dm),progress_pass=float(unit@dm)>=.5*distance,Phi=phi,DeltaPhi=phi-phi0,likelihood_competitive=phi-phi0<=2,noise_ball=phi<=chi2.ppf(.95,512)/2,whitened_pair_L2=float(np.linalg.norm(pair)),prior_cost=prior,DeltaPrior=dp,reasonable_prior_heuristic=dp<=10,legal=target.prior.legal(p.q),GT_used=False)

def substep(prior,base,p,v,J,residual,arm,g0,trust,distance):
 if arm not in ('B_FIXED_WEAK','C_DYNAMIC_WEAK','S2_GLOBAL'):return gn_step(prior,base.x,p.x,v,J,residual,trust=trust,material_limit=.5*distance)[0],{},None
 st=time.perf_counter()
 if arm=='C_DYNAMIC_WEAK':
  g=frame(prior,p.x,J);M=g['M'];Qa=g['Qweak'];Ualign,_,Vhalign=np.linalg.svd(Qa.T@g0['Qweak']);Qw=Qa@Ualign@Vhalign;dynamic_seconds=g['svd_seconds']
 else:
  M=real_material(material_derivative(prior,p.x));Qw=g0['Qglobal'] if arm=='S2_GLOBAL' else g0['Qweak'];dynamic_seconds=0.
 u=real_material(material_derivative(prior,base.x)@v);u/=np.linalg.norm(u)
 _,s,Vh=np.linalg.svd(M,full_matrices=False);keep=s>max(s[0]*1e-8,1e-10);E=Vh[keep].T
 radial=real_material(decode(prior,p.x)[0]-decode(prior,base.x)[0]);radial/=max(np.linalg.norm(radial),1e-30)
 constraints=np.vstack((Qw.T,u@M,radial@M));Z=E@null_space(constraints@E,rcond=1e-10)
 A=np.vstack((J@Z,Z));r=np.r_[residual,p.x];step=Z@(-np.linalg.lstsq(A,r,rcond=1e-10)[0])
 factor=min(1.,trust/max(np.linalg.norm(step),1e-300),.5*distance/max(np.linalg.norm(M@step),1e-300));step*=factor
 return step,dict(correction_rank=Z.shape[1],weak_constraint_linear_defect=float(np.linalg.norm(Qw.T@step)),Gaussian_Q_orthogonality_error=float(np.linalg.norm(Qw.T@Qw-np.eye(Qw.shape[1]))),dynamic_geometry_arithmetic_seconds=dynamic_seconds,dynamic_SVD_new_RHS=0,J_already_paid=True,substep_wall_seconds=time.perf_counter()-st),Qw

def radius_projection(prior,base,current,proposal,distance,Qw):
 if Qw is None:return match_distance(prior,base.x,lambda h:base.x+h*(proposal-base.x),distance)
 # Decoder-only Newton projection stays in Qs: exactly fixed weak coordinates.
 x=proposal.copy();m0=decode(prior,base.x)[0];Pstrong=np.eye(66)-Qw@Qw.T
 for k in range(30):
  dm=real_material(decode(prior,x)[0]-m0);norm=float(np.linalg.norm(dm));err=norm-distance
  if abs(err)<1e-10*distance:return x,dict(status='FIXED_WEAK_DECODER_RADIUS_MATCHED',iterations=k,actual_RMS=norm,finite_weak_coordinate_defect=float(np.linalg.norm(Qw.T@(x-current.x))))
  M=real_material(material_derivative(prior,x));grad=Pstrong@(M.T@(dm/max(norm,1e-30)));den=float(grad@grad)
  if den<1e-24:break
  delta=-err*grad/den
  if np.linalg.norm(delta)>.25:delta*=.25/np.linalg.norm(delta)
  improved=False
  for j in range(12):
   xn=x+(0.5**j)*delta
   if np.linalg.norm(xn-base.x)>4.:continue
   en=abs(np.linalg.norm(real_material(decode(prior,xn)[0]-m0))-distance)
   if en<abs(err):x=xn;improved=True;break
  if not improved:break
 return None,dict(status='FIXED_WEAK_RADIUS_PROJECTION_FAIL',iterations=k,finite_weak_coordinate_defect=float(np.linalg.norm(Qw.T@(x-current.x))))

def correct(target,base,y,xstart,v,distance,dest,arm,g0):
 unit=real_material(material_derivative(target.prior,base.x)@v);unit/=np.linalg.norm(unit);x=xstart.copy();best=None;bestobj=np.inf;trust=.5;previous=None;history=[];accepted=[base.x.copy()]
 for i in range(3):
  p=target.evaluate(x,arm+f'_trial{i}');rr=target.white(p.state.data-y);obj=.5*float(rr@rr+p.x@p.x);feas=feasible_progress(target.prior,base.x,p.x,unit,distance)
  rho=None
  if previous is not None:
   po,pred=previous;rho=(po-obj)/max(pred,1e-300) if pred>0 else -1.
   if rho<.25:trust*=.5
   elif rho>.75:trust=min(1.,trust*1.5)
  row=dict(iteration=i,endpoint=p.key,objective=obj,Phi=target.phi(p,y),feasible=feas,trust_radius=trust,actual_predicted_ratio=rho);history.append(row)
  if feas and obj<bestobj:best=p;bestobj=obj;accepted.append(p.x.copy())
  if i==2:break
  current=best if best is not None else p
  J=base.J if arm=='F_FROZENJ' else target.jacobian(current)
  r=target.white(current.state.data-y)
  step,extra,Qproject=substep(target.prior,base,current,v,J,r,arm,g0,trust,distance);row.update(extra)
  xn,ma=radius_projection(target.prior,base,current,current.x+step,distance,Qproject);row['projection']=ma
  if xn is None or not feasible_progress(target.prior,base.x,xn,unit,distance):row['stop']='NO_FEASIBLE_RETRACTED_STEP';break
  actual=xn-current.x;M=real_material(material_derivative(target.prior,current.x))
  for shrink in range(9):
   if np.linalg.norm(actual)<=trust*(1+1e-10) and np.linalg.norm(M@actual)<=.5*distance*(1+1e-10):break
   proposal=current.x+step*(.5**(shrink+1));xn,ma=radius_projection(target.prior,base,current,proposal,distance,Qproject)
   if xn is None:break
   actual=xn-current.x
  if xn is None or np.linalg.norm(actual)>trust*(1+1e-10) or np.linalg.norm(M@actual)>.5*distance*(1+1e-10) or not feasible_progress(target.prior,base.x,xn,unit,distance):row['stop']='POST_RETRACTION_TRUST_CONSTRAINT';break
  aug=np.r_[r,current.x];change=np.r_[J@actual,actual];pred=.5*float(aug@aug-np.linalg.norm(aug+change)**2)
  row.update(finite_weak_coordinate_defect=None if Qproject is None else float(np.linalg.norm(Qproject.T@actual)),predicted_reduction=pred,step_norm=float(np.linalg.norm(actual)),projection_induced_step_change=float(np.linalg.norm(actual-step)))
  if pred<=0 or np.linalg.norm(actual)<1e-9:row['stop']='NONPOSITIVE_PREDICTED_REDUCTION';break
  previous=(.5*float(r@r+current.x@current.x),pred);x=xn
 write(dest/'CORRECTION.json',dict(status='BUDGETED_CORRECTOR_NOT_CONVERGED_PROFILE',history=history,selected_endpoint=None if best is None else best.key))
 if best is not None:
  target.persist_cache(best,dest/'ENDPOINT.npz');save(dest/'ACTUAL_CORRECTION_PATH.npz',x_nodes=np.array(accepted))
 return best

def costprobe():
 dest=OUT/'costprobe'
 if dest.exists():raise RuntimeError('one paid costprobe only')
 dest.mkdir(parents=True);b=point(AUD/'R7_LOCKED_000/ONLINE_INPUT.npz');y=arrays(AUD/'R7_LOCKED_000/ONLINE_INPUT.npz')['observation']
 with Book('COSTPROBE',300) as book:
  target=Target(book)
  try:
   st=time.perf_counter();p=target.evaluate(b.x,'paid_costprobe_same_saved_anchor');forward=time.perf_counter()-st;st=time.perf_counter();target.jacobian(p);jwall=time.perf_counter()-st
   pred=float(np.linalg.norm(p.state.data-b.state.data));je=float(np.linalg.norm(p.J-b.J)/np.linalg.norm(b.J))
   v=np.random.default_rng(20261029).normal(size=66);v/=np.linalg.norm(v);eps=1e-5
   pp=target.evaluate(b.x+eps*v,'costprobe_FD_plus');pm=target.evaluate(b.x-eps*v,'costprobe_FD_minus');fd=(target.white(pp.state.data)-target.white(pm.state.data))/(2*eps);fde=float(np.linalg.norm(fd-b.J@v)/max(np.linalg.norm(b.J@v),1e-30))
   rec=dict(status='PASS' if max(pred,je)<1e-9 and fde<1e-3 else 'FAIL',forward_seconds=forward,J_seconds=jwall,prediction_error=pred,J_relative_error=je,finite_difference_relative_error=fde,actual=book.total,counters=dict(book.counts),projected_full_study_RHS=20280,GT_opened=False)
   write(dest/'RESULT.json',rec);print(json.dumps(rec),flush=True)
  finally:target.close()

def parent_run(parent,paths=True):
 if json.loads((OUT/'costprobe/RESULT.json').read_text())['status']!='PASS':raise RuntimeError('costprobe gate')
 dest=OUT/'parents'/parent
 if (dest/'COMPLETE.json').exists():raise RuntimeError('no automatic parent replay')
 dest.mkdir(parents=True,exist_ok=True);idx=int(parent.rsplit('_',1)[1]);inp=AUD/parent/'ONLINE_INPUT.npz';b=point(inp);y=arrays(inp)['observation'];endpoints={};rows=[]
 with Book(parent,1800) as book:
  target=Target(book)
  try:
   g=frame(target.prior,b.x,b.J);g['Qglobal']=arrays(AUD/'S2_GLOBAL_FRAME.npz')['Qweak'];gg=arrays(inp);v=g['v'];directions={a:v for a in ARMS};meta=json.loads((AUD/parent/'BASIS_METADATA.json').read_text())
   for kind in ['FOURIER','WAVELET']:directions['D_'+kind]=gg['D_'+kind]
   directions['E_SOM']=None if meta['E_SOM']['status']=='MATERIAL_LIFT_FAIL' else gg['E_SOM'];directions['R_RANDOM']=gg['random_v'];vg=g['Qglobal'][:,-1].copy();vg/=np.linalg.norm(g['M']@vg)
   if vg[np.argmax(abs(vg))]<0:vg=-vg
   directions['S2_GLOBAL']=vg
   save(dest/'BASE.npz',x=b.x,material=b.state.chi,prediction=b.state.data,J=b.J,current=b.state.current,field=b.state.field,observation=y,**g)
   for radius in [.01,.02]:
    order=ARMS[idx%len(ARMS):]+ARMS[:idx%len(ARMS)]
    for arm in order:
     dd=dest/arm/f'RMS_{radius:.4f}';dd.mkdir(parents=True,exist_ok=True);before=book.total.copy();st=time.perf_counter();vv=directions[arm]
     if vv is None:row=dict(parent=parent,arm=arm,radius=radius,**meta[arm])
     else:
      xp,match=match_distance(target.prior,b.x,lambda h:b.x+h*vv,radius)
      if xp is None:row=dict(parent=parent,arm=arm,radius=radius,status=match['status'])
      else:
       if arm=='U_UNCORRECTED':pp=target.evaluate(xp,arm);target.persist_cache(pp,dd/'ENDPOINT.npz')
       else:pp=correct(target,b,y,xp,vv,radius,dd,arm,g)
       if pp is None:row=dict(parent=parent,arm=arm,radius=radius,status='NO_FEASIBLE_ENDPOINT')
       else:
        unit=real_material(material_derivative(target.prior,b.x)@vv);unit/=np.linalg.norm(unit)
        row=metrics(target,pp,b,y,unit,radius);row.update(parent=parent,arm=arm,radius=radius,status='EXECUTED',initial_match=match,basis_metadata=meta.get(arm,{}))
        row['qualified']=all(row[k] for k in ['distance_matched','progress_pass','likelihood_competitive','noise_ball','reasonable_prior_heuristic','legal'])
        if radius==.02 and arm in ['B_FIXED_WEAK','C_DYNAMIC_WEAK']:endpoints[arm]=pp
     operational_end=book.total.copy();operational_wall=time.perf_counter()-st
     row['operational_cost']={k:operational_end[k]-before[k] for k in ['RHS','states','LU','forward_RHS','derivative_RHS']};row['operational_wall_seconds']=operational_wall
     book._event('ALGORITHM_SELECTION_CLOSED',arm=arm,radius=radius,operational_cost=row['operational_cost'],operational_wall_seconds=operational_wall)
     audit_start=time.perf_counter()
     if row['status']=='EXECUTED':
      if pp.J is None:
       if pp is target._live_point:target.jacobian(pp)
       else:
        restored=target.evaluate(pp.x,arm+'_selected_J_AUDIT_restore');target.jacobian(restored);pp.J=restored.J.copy();pp.Jnative=restored.Jnative.copy();pp.T=restored.T.copy()
      target.persist_cache(pp,dd/'ENDPOINT.npz')
     row['diagnostic_cost']={k:book.total[k]-operational_end[k] for k in ['RHS','states','LU','forward_RHS','derivative_RHS']};row['diagnostic_wall_seconds']=time.perf_counter()-audit_start
     row.update(cost={k:book.total[k]-before[k] for k in ['RHS','states','LU','forward_RHS','derivative_RHS']},wall_seconds=time.perf_counter()-st)
     write(dd/'RESULT.json',row);rows.append(row);print('ENDPOINT',parent,arm,radius,row.get('Phi'),row['status'],book.total['RHS'],flush=True)
   pathrows=[]
   if paths:
    for arm,p_end in endpoints.items():
     pd=dest/'paths'/arm;pd.mkdir(parents=True,exist_ok=True);nodes=[]
     xp=arrays(dest/arm/'RMS_0.0200/ACTUAL_CORRECTION_PATH.npz')['x_nodes'];mp=[decode(target.prior,x)[0] for x in xp];lengths=np.r_[0.,np.cumsum([np.linalg.norm(real_material(mp[i+1]-mp[i])) for i in range(len(mp)-1)])]
     def path_x(t):
      ell=t*lengths[-1];k=min(np.searchsorted(lengths,ell,side='right')-1,len(xp)-2);k=max(k,0);h=(ell-lengths[k])/max(lengths[k+1]-lengths[k],1e-30);return (1-h)*xp[k]+h*xp[k+1]
     for t in [0,.25,.5,.75,1]:
      if t==0:p=b
      elif t==1:p=p_end
      else:p=target.evaluate(path_x(t),arm+f'_path_{t}');target.jacobian(p)
      gf=frame(target.prior,p.x,p.J);dist=float(np.linalg.norm(real_material(p.state.chi-b.state.chi)));unit=real_material(material_derivative(target.prior,b.x)@v);unit/=np.linalg.norm(unit)
      row=metrics(target,p,b,y,unit,dist);row.update(t=t,material_effective_rank=gf['rank'],weak_material_projector_operator2_drift=projector_drift(g['Qw'],gf['Qw']))
      save(pd/f'NODE_{t:.2f}.npz',x=p.x,material=p.state.chi,prediction=p.state.data,current=p.state.current,field=p.state.field,J=p.J,Qw=gf['Qw'],material_singular=gf['material_singular']);nodes.append(row)
     pr=dict(parent=parent,arm=arm,nodes=nodes,sampled_path_qualified=all(r['likelihood_competitive'] and r['noise_ball'] and r['reasonable_prior_heuristic'] for r in nodes),continuous_claim=False,path_definition='actual accepted correction polyline; uniform material arclength node placement',correction_x_nodes=xp.tolist());write(pd/'PATH.json',pr);pathrows.append(pr)
   write(dest/'COMPLETE.json',dict(parent=parent,role='EXPOSED_DEV_ONLY',endpoint_results=rows,path_results=pathrows,cumulative_ledger=book.total,online_GT_used=False,fresh_locked_evidence=False))
  finally:target.close()

def e1():
 prior=prior_only();rows=[]
 for idx in range(8):
  parent=f'R7_LOCKED_{idx:03d}';b=point(AUD/parent/'ONLINE_INPUT.npz');g0=frame(prior,b.x,b.J);nodes=[]
  for t in [0,.25,.5,.75,1]:
   path=H7/'E2/parents'/parent/'paths/B_FULLJ_CHORD'/f'NODE_{t:.2f}.npz'
   if not path.exists():nodes.append(dict(t=t,status='HISTORICAL_NODE_MISSING'));continue
   a=arrays(path);g=frame(prior,a['x'],a['J']);cos=np.linalg.svd(g0['Qweak'].T@g['Qweak'],compute_uv=False)
   nodes.append(dict(t=t,status='REUSED_ARRAYS',weak_Gaussian_principal_cosines=cos.tolist(),weak_max_angle_degrees=float(np.degrees(np.arccos(np.clip(cos.min(),-1,1)))),material_rank=g['rank'],material_RMS_from_anchor=float(np.linalg.norm(real_material(a['material']-b.state.chi))),new_RHS=0))
  rows.append(dict(parent=parent,historical_path='R7 B_FULLJ_CHORD; historical geometry diagnostic only',nodes=nodes))
 write(AUD/'E1_HISTORICAL_SUBSPACE_STABILITY.json',dict(status='COMPLETE',new_RHS=0,rows=rows,scientific_scope='exposed DEV historical arrays; no prospective evidence',prior_orthonormal_frame=True))
 print('E1_COMPLETE zeroRHS',flush=True)

def summarize():
 rows=[];parents=[];paths=[]
 for p in sorted((OUT/'parents').glob('*/COMPLETE.json')):
  a=json.loads(p.read_text());rows+=a['endpoint_results'];paths+=a['path_results'];parents.append(a['parent'])
 table=[]
 for radius in [.01,.02]:
  for arm in ARMS:
   a=[r for r in rows if r['arm']==arm and r['radius']==radius];ok=[r for r in a if r['status']=='EXECUTED']
   table.append(dict(arm=arm,radius=radius,n=len(a),available=len(ok),qualified=sum(r.get('qualified',False) for r in ok),likelihood_competitive=sum(r['likelihood_competitive'] for r in ok),pair_le1=sum(r['whitened_pair_L2']<=1 for r in ok),median_DeltaPhi=None if not ok else float(np.median([r['DeltaPhi'] for r in ok])),median_pair=None if not ok else float(np.median([r['whitened_pair_L2'] for r in ok])),RHS=sum(r['cost']['RHS'] for r in a),operational_RHS=sum(r['operational_cost']['RHS'] for r in a),diagnostic_RHS=sum(r['diagnostic_cost']['RHS'] for r in a),operational_wall_seconds=sum(r['operational_wall_seconds'] for r in a),states=sum(r['cost']['states'] for r in a),wall_seconds=sum(r['wall_seconds'] for r in a)))
 write(OUT/'SUMMARY.json',dict(parents=parents,table=table,paths=paths,rows=rows,role='DEV_ONLY',GT_online=False,new_locked_evidence=False))
 lines=['# R9-S exploration (exposed development objects only)','', '| Arm | r | Available | Qualified | dPhi<=2 | pair<=1 | median dPhi | median pair | RHS |','|---|---:|---:|---:|---:|---:|---:|---:|---:|']
 for r in table:lines.append('|'+ '|'.join(str(r[k]) for k in ['arm','radius','available','qualified','likelihood_competitive','pair_le1','median_DeltaPhi','median_pair','RHS'])+'|')
 lines+=['','R7 material metric, prior, noise ball, likelihood and signed-progress gates preserved. These eight historical objects are exposed DEV. No posterior, exact fiber, or continuous-between-node certificate follows.','SOM current decomposition requires a shared-material decoder lift; MATERIAL_LIFT_FAIL remains failure.']
 (OUT/'RESULTS.md').write_text('\n'.join(lines)+'\n');print('SUMMARY',len(parents),'parents',flush=True)

def main():
 pa=argparse.ArgumentParser();pa.add_argument('command',choices=['e0','e1','costprobe','parent','all','summarize']);pa.add_argument('--parent');pa.add_argument('--no-paths',action='store_true');a=pa.parse_args()
 if a.command=='e0':e0()
 elif a.command=='e1':e1()
 elif a.command=='costprobe':costprobe()
 elif a.command=='parent':parent_run(a.parent,not a.no_paths)
 elif a.command=='all':
  for i in range(8):parent_run(f'R7_LOCKED_{i:03d}',not a.no_paths);summarize()
 else:summarize()
if __name__=='__main__':main()
