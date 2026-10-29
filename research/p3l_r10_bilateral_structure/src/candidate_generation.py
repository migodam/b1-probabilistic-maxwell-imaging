"""Frozen matched randomized inverse clouds. Exact historical ROI physics, no posterior assertion."""
from __future__ import annotations
import os
for n in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'):os.environ[n]='1'
import argparse,fcntl,hashlib,json,sys,time,traceback
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
P=Path(__file__).resolve().parents[1]
H=P.parents[2]/'b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace'
sys.path.insert(0,str(H/'src'))
import possibility_volume as pv
CONFIG=P/'config/candidate_protocol.json'
ENABLE=P/'config/candidate_enable.json'
OUT=P/'results/candidates'
ARMS=['FULL66_RML','SHAPE_TARGET_RML_BLOCK','BILATERAL_TARGET_RML_BLOCK']
COUNTS=('reduced_RHS256','reduced_LU256','full_RHS576','full_LU576','actual_forward_RHS256','actual_J_RHS256','actual_LU256','completed_forward_RHS256','completed_J_RHS256','failed_actions','array_saves','array_bytes','physics_wall_seconds','save_wall_seconds')
class Stop(RuntimeError):pass

def write(p,a):pv.write(p,a)
def ident(p):
 p=Path(p);s=p.stat();return dict(path=str(p),bytes=s.st_size,mtime_ns=s.st_mtime_ns)
def source_id():return hashlib.blake2b(Path(__file__).read_bytes(),digest_size=32).hexdigest()
def freeze():
 if CONFIG.exists():raise RuntimeError('protocol already exists')
 c=dict(id='R10_CANDIDATE_DEV_V1',approved_by_root=False,created_utc=time.time(),source_blake2b=source_id(),historical_root=str(H),output=str(OUT),arms=ARMS,K=16,parents=[f'R9V_FRESH_{i:03d}' for i in range(8)],exposure='DEV first8 previously exposed R9fresh32; neither fresh confirmation nor posterior validation',seed_namespace=[20261030,1000701],candidate=dict(max_forward=12,max_J=8,initial_trust=2.,trust_min=.03125,trust_max=4.,damping_initial=.01,line_search=[1.,.5,.25],gradient_stop=1e-6,step_stop=1e-7),caps=dict(reduced_RHS256=250000,reduced_LU256=7500,full_RHS576=0,full_LU576=0,controlled_wall_seconds=900),shared_setup=dict(path=str(H/'results/stable_roi/SETUP.npz'),new_setup_charge=0,already_paid_setup_full_RHS576=264),noise_ball=282.8737921622988,epsilon_primary=.075,qualified_endpoint='legal and true observed Phi<=noise_ball; no displacement gate',selection='minimum own randomized objective among all visited states, including anchor; no GT',objective='0.5||white(F-y)-eta512||²+0.5||x-xi66||²; shape arms xi[64:]=0; bilateral arm xi[:64]=L_B@xi[:64]',profile='NOT_RUN coupled nonlinear profile; 2x2 amplitude Schur is ordinary joint66 GN block solve, every arm identical; differences vs full66 arise from randomized target distribution',bilateral='fixed host cache ER,O; x=0 decoderD; H=2/sigma² Re(D*[(O*O) hadamard(conj(ER)ER.T)]D); inverse(I+H/trace(H)*64) eigenvalues clip[.25,4], renormalize trace64; retain all64 directions',fairness='same16 common random draws, sameall66 GN freedom and same12F8J caps for all3; identical Schur solver in allarms; interleaved proposal order; all rejected/failed trials charged',prior_coordinate='x64:66 Gaussian via q64:66=normalCDF(x); originalGaussian64xUniform2 prior pushforward',resources=dict(BLAS_threads=1,single_worker=True),truth='scorer separate after PREDICTIONS_SEALED; online never reads truth',inputs=[ident(H/'results/possibility_volume/inputs/online'/f'R9V_FRESH_{i:03d}.npz') for i in range(8)])
 write(CONFIG,c)
 (P/'METHOD_CANDIDATES_AND_SELECTION.md').write_text('''# Candidate mechanisms frozen before physics\n\nAll three mechanisms address direction coverage rather than enlarging the old weak8 plane. The strong baseline is ordinary full66 randomized maximum likelihood with independent real standard measurement noise and Gaussian prior centers. Shape64 randomization isolates the material field from two global inclusion amplitudes: amplitudes are optimized conditionally by exact block elimination of the local Gauss–Newton system. That block elimination is algebraically joint66 Gauss–Newton; it is not a new algorithm or exact nonlinear amplitude profiling. The Schur implementation is used identically in all arms, so differences from FULL66 are caused by target distributions. A direct joint solve of the same damped normal system is the zero-physics code-review reference. Coupled nonlinear amplitude profiling is NOT_RUN. The third mechanism changes the full64 shape-center covariance using a fixed host source-and-receiver Gram with correct real decoder pullback; covariance retains all64 support and trace64.\n\nSelection: implement all three because they provide a strong ordinary baseline, an interpretable shape/amplitude randomization ablation, and a bilateral full-support geometry ablation. Randomized minimizers are finite heuristic possibilities, not posterior samples. No Jacobian determinant density or RTO claim is made.\n\nEvery proposal has the same full66 corrector freedom, 12 paid forward states and8 paid fullJ states maximum; exact cached ROI256 backend, identical common random draws, K16, shared old anchor, true-observation noise gate. No truth guides endpoint choice. Cap hit and unexecuted work are right-censored. DEV data are first8 of previously exposed R9fresh32. Source and protocol freeze precede a separate root run token.\n''')
 print('FROZEN',CONFIG,flush=True)

class Book:
 def __init__(self,c):
  self.c=c;self.job='CANDIDATE_DEV';self.context={};self.jobdir=OUT/'physics';self.jobdir.mkdir(parents=True,exist_ok=True);self.started=time.perf_counter();self.cpu0=time.process_time();self.total={k:0 for k in COUNTS};self.seq=0
  self.f=(OUT/'RESOURCE_LEDGER.jsonl').open('x');self.lock=(P/'candidate_physics_worker.lock').open('a');fcntl.flock(self.lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  self.event('JOB_BEGIN',pid=os.getpid(),prior_backend='exact cached ROI256',setup_new_RHS576=0,setup_historical_RHS576=264)
 def event(self,event,**kw):
  for k in COUNTS:self.total[k]+=kw.get(k,0)
  self.seq+=1;row=dict(sequence=self.seq,utc=time.time(),event=event,context=self.context.copy(),**kw,after=self.total.copy());self.f.write(json.dumps(pv.plain(row),allow_nan=False)+'\n');self.f.flush()
 def check(self):
  if time.perf_counter()-self.started>=self.c['caps']['controlled_wall_seconds']:raise Stop('controlled900s wallcap')
 def reserve_forward(self):
  self.check()
  if self.total['reduced_RHS256']+8>self.c['caps']['reduced_RHS256'] or self.total['reduced_LU256']+1>self.c['caps']['reduced_LU256']:raise Stop('physical hardcap')
  self.event('FORWARD_RESERVED',reduced_RHS256=8,reduced_LU256=1)
 def reserve_J(self):
  self.check()
  if self.total['reduced_RHS256']+32>self.c['caps']['reduced_RHS256']:raise Stop('physical hardcap')
  self.event('J_RESERVED',reduced_RHS256=32)
 def save(self,path,**data):
  t=time.perf_counter();pv.save(path,**data);self.event('ARRAY_SAVE',array_saves=1,array_bytes=Path(path).stat().st_size,save_wall_seconds=time.perf_counter()-t)
 def close(self):
  self.event('JOB_END',wall_seconds=time.perf_counter()-self.started,CPU_seconds=time.process_time()-self.cpu0);write(OUT/'RESOURCE_SNAPSHOT.json',dict(**self.total,wall_seconds=time.perf_counter()-self.started,CPU_seconds=time.process_time()-self.cpu0));self.f.close();self.lock.close()

def bilateral(ph):
 t=time.perf_counter();D=pv.derivative(ph.prior,np.zeros(66))[ph.idx,:64];O=ph.a['O'];E=ph.a['ER'];Gr=O.conj().T@O;Gs=E.conj()@E.T;Hb=(2/pv.SIGMA**2)*np.real(D.conj().T@(Gr*Gs)@D);Hb=(Hb+Hb.T)/2
 ev,V=eigh(Hb);scale=max(float(np.mean(np.maximum(ev,0))),1e-20);lam=np.clip(1/(1+np.maximum(ev,0)/scale),.25,4.);lam*=64/lam.sum();L=(V*np.sqrt(lam))@V.T
 Hs=(2/pv.SIGMA**2)*np.real(D.conj().T@Gs@D);Hr=(2/pv.SIGMA**2)*np.real(D.conj().T@Gr@D)
 pv.save(OUT/'BILATERAL_FIXED_PRIOR.npz',source_only_H=Hs,receiver_only_H=Hr,H=Hb,eigenvalues=ev,covariance_eigenvalues=lam,center_transform=L,decoder_real=D.real,decoder_imag=D.imag)
 return L,dict(wall_seconds=time.perf_counter()-t,live_physics_solves=0,description='cached host-field algebra only; no nonlinear forward/J; all64 support retained',covariance_trace=float(lam.sum()),covariance_min=float(lam.min()),covariance_max=float(lam.max()),controls='source_only_H and receiver_only_H stored; normalized trace controls and cached nonlinearJ alignment written after parent anchor load')

def schur(H,g):
 # Same full66 solve in every arm; exact block elimination, no freedom difference.
 A=H[:64,:64];B=H[:64,64:];C=H[64:,64:];CiBt=np.linalg.solve(C,B.T);Cig=np.linalg.solve(C,g[64:]);s=np.linalg.solve(A-B@CiBt,-g[:64]+B@Cig);a=np.linalg.solve(C,-g[64:]-B.T@s);return np.r_[s,a]

def candidate(ph,base,y,eta,xi,c,parent,arm,k):
 cfg=c['candidate'];t=time.perf_counter();cpu=time.process_time();begin=ph.book.total.copy();p=base;best=base;nf=nj=0;trust=cfg['initial_trust'];damp=cfg['damping_initial'];history=[]
 def objective(pt):return .5*float(np.sum((pv.white(pt.prediction-y)-eta)**2)+np.sum((pt.x-xi)**2))
 bestobj=obj=objective(p);status='MAX_ITERATIONS';reason=None
 try:
  while nf<cfg['max_forward'] and nj<cfg['max_J']:
   ph.book.check()
   if p.J is None and p is not ph.live:
    if nf>=cfg['max_forward']:break
    p=ph.evaluate(p.x,f'{parent}_{arm}_{k}_restore');nf+=1
   J=ph.jacobian(p);nj+=1;res=pv.white(p.prediction-y)-eta;g=J.T@res+(p.x-xi)
   if np.max(abs(g))<cfg['gradient_stop']:status='GRADIENT_STOP';break
   H=J.T@J+np.eye(66);step=schur(H+damp*np.eye(66),g);sn=np.linalg.norm(step)
   if sn>trust:step*=trust/sn
   if np.linalg.norm(step)<cfg['step_stop']:status='STEP_STOP';break
   accepted=False
   for alpha in cfg['line_search']:
    if nf>=cfg['max_forward']:break
    q=ph.evaluate(p.x+alpha*step,f'{parent}_{arm}_{k}_trial');nf+=1;qo=objective(q);history.append(dict(forward=nf,J=nj,objective=qo,Phi=ph.phi(q,y),alpha=alpha,accepted=bool(qo<obj),xnorm=float(np.linalg.norm(q.x))))
    if qo<bestobj:best=q;bestobj=qo
    if qo<obj:p=q;obj=qo;accepted=True;break
   if accepted:trust=min(cfg['trust_max'],trust*1.5);damp=max(1e-6,damp*.5)
   else:trust*=.5;damp*=10
   if trust<cfg['trust_min']:status='TRUST_STOP';break
  if nf>=cfg['max_forward']:status='FORWARD_CAP'
  elif nj>=cfg['max_J']:status='J_CAP'
 except Stop as e:status='RIGHT_CENSORED_GLOBAL_CAP';reason=str(e)
 except Exception as e:status='FAILED';reason=repr(e);ph.book.event('CANDIDATE_FAILURE',failed_actions=1,error=reason,traceback=traceback.format_exc())
 dest=OUT/'parents'/parent/arm/f'{k:02d}';dest.mkdir(parents=True,exist_ok=True);ph.persist(best,dest/'ENDPOINT.npz');pv.save(dest/'RANDOMIZATION.npz',eta=eta,xi=xi)
 phi=ph.phi(best,y);row=dict(parent=parent,arm=arm,k=k,status=status,error=reason,forward_trials=nf,J_states=nj,Phi=phi,legal=ph.prior.legal(best.q),qualified=bool(ph.prior.legal(best.q) and phi<=c['noise_ball']),random_objective=bestobj,anchor_random_objective=objective(base),anchor_selected=bool(np.array_equal(best.x,base.x)),material_RMS_from_anchor=ph.prior.material_rms(best.material,base.material),GT_used=False,posterior_sample=False,path=str(dest/'ENDPOINT.npz'),wall_seconds=time.perf_counter()-t,CPU_seconds=time.process_time()-cpu,cost={key:ph.book.total[key]-begin[key] for key in COUNTS},history=history)
 write(dest/'RECORD.json',row);print('CANDIDATE',parent,arm,k,status,round(phi,3),row['qualified'],nf,nj,flush=True);return row

def geometry_controls(base,parent):
 """Only cached old J and fixed host arrays; no new physics, no numerical solve test."""
 z=pv.load_npz(OUT/'BILATERAL_FIXED_PRIOR.npz');matrices={'source_only':z['source_only_H'],'receiver_only':z['receiver_only_H'],'bilateral':z['H'],'cached_nonlinear_J':base.J[:,:64].T@base.J[:,:64]};N={};V={};out={}
 for name,h in matrices.items():
  h=(h+h.T)/2;ev,v=eigh(h);N[name]=h/max(float(np.trace(h)),1e-300);V[name]=v[:,-16:];out[name]=dict(trace=float(np.trace(h)),eigenvalues=ev.tolist(),effective_rank=float(np.trace(h)**2/max(float(np.sum(h*h)),1e-300)))
 for name in ('source_only','receiver_only','bilateral'):
  a=N[name];b=N['cached_nonlinear_J'];out[name].update(cached_J_normalized_Frobenius_cosine=float(np.sum(a*b)/max(np.linalg.norm(a)*np.linalg.norm(b),1e-300)),cached_J_top16_overlap=float(np.linalg.norm(V[name].T@V['cached_nonlinear_J'])**2/16))
 write(OUT/'parents'/parent/'GEOMETRY_CONTROLS.json',dict(parent=parent,new_physics_solves=0,cache='old R9 anchor J; historical charge retained, not a newly computed J',metrics=out,block_solver_review='schur(H,g) returns solve(H,-g) by exact2x2 elimination when invertible; no additional physical or microscopic identity experiment',nonlinear_profile='NOT_RUN'))


def run():
 c=json.loads(CONFIG.read_text());en=json.loads(ENABLE.read_text())
 if not c['approved_by_root'] or en.get('protocol_id')!=c['id'] or en.get('enabled') is not True:raise RuntimeError('root run token missing')
 if c['source_blake2b']!=source_id():raise RuntimeError('source changed after freeze')
 if (OUT/'PREDICTIONS_SEALED.json').exists():raise RuntimeError('no replay')
 book=Book(c);ph=pv.ExactROI(book);rows=[];stopped=False
 try:
  L,bmeta=bilateral(ph);write(OUT/'BILATERAL_SETUP.json',bmeta)
  for pi,parent in enumerate(c['parents']):
   y=pv.online_observation(dict(output=str(H/'results/possibility_volume')),parent);base=pv.load_anchor(H/'results/possibility_volume/parents'/parent/'anchor/ANCHOR.npz');geometry_controls(base,parent);rng=np.random.default_rng(np.random.SeedSequence(c['seed_namespace']+[pi]));draws=[(rng.normal(size=512),rng.normal(size=66)) for _ in range(c['K'])]
   for k,(eta,rawxi) in enumerate(draws):
    for arm in c['arms']:
     book.context=dict(parent=parent,arm=arm,k=k)
     if stopped:
      row=dict(parent=parent,arm=arm,k=k,status='NOT_RUN_RIGHT_CENSORED',qualified=False,GT_used=False,path=None,cost={key:0 for key in COUNTS});rows.append(row);continue
     xi=rawxi.copy()
     if arm!='FULL66_RML':xi[64:]=0
     if arm=='BILATERAL_TARGET_RML_BLOCK':xi[:64]=L@xi[:64]
     row=candidate(ph,base,y,eta,xi,c,parent,arm,k);rows.append(row)
     if row['status']=='RIGHT_CENSORED_GLOBAL_CAP':stopped=True
   write(OUT/'parents'/parent/'COMPLETE.json',dict(parent=parent,rows=[r for r in rows if r['parent']==parent],GT_used=False))
  write(OUT/'PREDICTIONS_SEALED.json',dict(protocol_id=c['id'],utc=time.time(),GT_used=False,rows=rows,stopped=stopped,source_blake2b=source_id(),dataset_exposure=c['exposure']))
 finally:ph.close();book.close()

if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('action',choices=['freeze','run']);args=a.parse_args();globals()[args.action]()
