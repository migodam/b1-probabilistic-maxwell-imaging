"""Actual nonlinear shape-path reuse adjudication; no GT, no material coverage scoring."""
from __future__ import annotations
import os
for n in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'):os.environ[n]='1'
import json,time,csv,sys,hashlib,fcntl,datetime,traceback
from pathlib import Path
import numpy as np
from scipy.linalg import eigh,schur,solve_triangular,lu_factor,lu_solve
D=Path(__file__).resolve().parent;P=D.parent;sys.path.insert(0,str(P/'src'));import candidate_generation as cg;import reduced_response as dev
pv=cg.pv;SIGMA=pv.SIGMA;ARMS=['DIRECT_LU','REBUILD_R64','REUSE_ANCHOR_R64']
class Stop(RuntimeError):pass

def write(p,x):pv.write(p,x)
def arr(p):return pv.load_npz(p)
class Book:
 def __init__(self,c):
  self.c=c;self.f=(D/'COST_LEDGER.jsonl').open('x');self.start=time.perf_counter();self.cpu=time.process_time();self.rhs={};self.lu=0;self.matvec=0;self.seq=0;self.context={};self.lock=(P/'candidate_physics_worker.lock').open('a');fcntl.flock(self.lock,fcntl.LOCK_EX|fcntl.LOCK_NB);self.event('BEGIN',pid=os.getpid())
 def check(self):
  if time.time()>=self.c['stop_new_compute_unix'] or time.perf_counter()-self.start>=self.c['wall_cap_seconds']:raise Stop('absolute19:06orcontrolledwallSTOP')
 def event(self,event,**kw):
  self.seq+=1;self.f.write(json.dumps(pv.plain(dict(sequence=self.seq,event=event,utc=time.time(),context=self.context.copy(),**kw)),allow_nan=False)+'\n');self.f.flush()
 def reserve(self,order,label,lu=0):
  self.check()
  if sum(self.rhs.values())+8>10000 or self.lu+lu>600:raise Stop('RHSorLUcap')
  self.rhs[order]=self.rhs.get(order,0)+8;self.lu+=lu;self.event('SOLVE_RESERVED',order=order,RHS=8,LU256=lu,label=label,after_RHS=sum(self.rhs.values()),after_LU256=self.lu)
 def close(self,status,error=None):
  x=dict(status=status,error=error,RHS_by_order=self.rhs,RHS_conservative_total=sum(self.rhs.values()),LU256=self.lu,RHS576=0,Green_columns_applied=self.matvec,wall_seconds=time.perf_counter()-self.start,CPU_seconds=time.process_time()-self.cpu);self.event('END',**x);write(D/'RECEIPT.json',x);self.f.close();self.lock.close()

def formation(prior,setup,x):
 t=time.perf_counter();q=pv.q_from_x(x);material,meta=prior.decode(q);a=meta['fraction'];p=np.sqrt(a[setup['ROI']]);B=p[:,None]*setup['C']*p[None,:];L=setup['O']*p[None,:];R=p[:,None]*setup['ER'];return B,L,R,a,q,material,time.perf_counter()-t

def reduced(B,L,R,V,b):
 t=time.perf_counter()
 if V.shape[1]<64:return None,time.perf_counter()-t
 b.matvec+=64;BV=B@V;H=V.conj().T@BV;T,Q=schur(H,output='complex',check_finite=False);VQ=V@Q;BVQ=BV@Q;LQ=L@VQ;QR=Q.conj().T@(V.conj().T@R);return (T,VQ,BVQ,LQ,QR),time.perf_counter()-t

def run():
 c=json.loads((D/'CONTRACT.json').read_text());assert c['source_blake2b']==hashlib.blake2b(Path(__file__).read_bytes(),digest_size=32).hexdigest();b=Book(c);rows=[];setups=[];anchors=[];status='FAILED';error=None
 try:
  setup=arr(cg.H/'results/stable_roi/SETUP.npz');cfg=json.loads((pv.R4/'research/p3l_r4_conditional_transport/configs/TARGET.json').read_text())
  with np.load(pv.R4/'research/p3l_r4_conditional_transport/inputs/background/weak.npz') as z:points=z['points'].copy()
  prior=pv.MaterialPrior(cfg,points);I256=np.eye(256,dtype=complex)
  for pi,parent in enumerate(c['parents']):
   b.check();b.context=dict(parent=parent,role='ANCHOR_INITIAL_BASIS');base=arr(P/'results/candidates_fresh10/parents'/parent/'anchor/ANCHOR.npz');r0=json.loads((P/'results/candidates_fresh10/parents'/parent/'FULL66_RML/00/RECORD.json').read_text());x0=base['x'];x1=arr(r0['path'])['x'];y=pv.online_observation(dict(output=str(P/'results/candidates_fresh10')),parent);deltas=arr(P/'results/reduced_response_fresh10'/f'{parent}.npz')['delta'][:8]
   B0,L0,R0,a0,q0,m0,anchorform=formation(prior,setup,x0);ts=time.perf_counter();V0,widths0=dev.basis(B0,R0,64,1e-12,b);anchorbasis=time.perf_counter()-ts;anchors.append(dict(parent=parent,formation_seconds=anchorform,basis64_seconds=anchorbasis,total_reuse_initial_seconds=anchorform+anchorbasis,rank=int(V0.shape[1])));pv.save(D/f'{parent}_ANCHOR.npz',x0=x0,x1=x1,V0=V0,delta=deltas,fraction=a0)
   for ti,tval in enumerate(c['t']):
    b.check();b.context=dict(parent=parent,t=tval,role='SHAPE_SETUP');x=(1-tval)*x0+tval*x1;B,L,R,a,q,material,form=formation(prior,setup,x);assert prior.legal(q)
    ts=time.perf_counter();beta=float(eigh((B-B.conj().T)/(2j),eigvals_only=True,subset_by_index=(0,0),check_finite=False)[0]);Lnorm=float(np.linalg.norm(L,'fro'));energy=time.perf_counter()-ts
    ts=time.perf_counter();Vnew,widths=dev.basis(B,R,64,1e-12,b);newbasis=time.perf_counter()-ts;rebuild,rbmodel=reduced(B,L,R,Vnew,b);reuse,rumodel=reduced(B,L,R,V0,b);ss=dict(parent=parent,t=tval,formation_seconds=form,beta_Lnorm_seconds=energy,beta=beta,Lnorm=Lnorm,rebuild_basis_seconds=newbasis,rebuild_model_seconds=rbmodel,reuse_model_seconds=rumodel,new_basis_rank=int(Vnew.shape[1]),reused_basis_rank=int(V0.shape[1]),shape_legal=True,material_RMS_from_anchor=float(np.linalg.norm(material-m0)/24));setups.append(ss);b.event('SHAPE_SETUP_COMPLETE',**ss);raw={arm:[] for arm in ARMS};approx={arm:[] for arm in ARMS[1:]};node_material=[];node_q=[]
    def direct(d,label):
     b.reserve(256,label,1);ts=time.perf_counter();fac=lu_factor(I256-d*B,check_finite=False);u=lu_solve(fac,R,check_finite=False);F=setup['Fh']+d*(L@u);dt=time.perf_counter()-ts;b.event('SOLVE_COMPLETE',label=label,order=256,wall_seconds=dt);return F,dt
    for ni,d in enumerate(deltas):
     b.context=dict(parent=parent,t=tval,node=ni);order=ARMS[(pi+ti+ni)%3:]+ARMS[:(pi+ti+ni)%3];out={};meta={}
     for arm in order:
      if arm=='DIRECT_LU':out[arm],dt=direct(d,arm);meta[arm]=dict(pure_seconds=dt,total_seconds=dt,gate=True,fallback=False,bound=0.,potential=0.,residual=0.,coercivity=None)
      else:
       model=rebuild if arm=='REBUILD_R64' else reuse;gate=False;bound=pot=rn=bn=None;solve=gt=fbtime=0.;fa=setup['Fh'].copy()
       if model is not None:
        T,VQ,BVQ,LQ,QR=model;b.reserve(64,arm);ts=time.perf_counter();u=solve_triangular(np.eye(64)/d-T,QR,check_finite=False);fa=setup['Fh']+LQ@u;solve=time.perf_counter()-ts;b.event('SOLVE_COMPLETE',label=arm,order=64,wall_seconds=solve);ts=time.perf_counter();res=R-(VQ@u)/d+BVQ@u;rn=float(np.linalg.norm(res));bn=float(d.imag/abs(d)**2+beta)
        if bn>0:bound=float(np.sqrt(2)/SIGMA*Lnorm*rn/bn);pot=float(np.sqrt(2)/SIGMA*np.linalg.norm(fa-y)*bound+.5*bound**2);gate=bool(bound<=.01 and pot<=.1)
        gt=time.perf_counter()-ts
       if gate:out[arm]=fa
       else:out[arm],fbtime=direct(d,arm+'_EXACT_FULL_LU_FALLBACK')
       approx[arm].append(fa);meta[arm]=dict(pure_seconds=solve+gt,total_seconds=solve+gt+fbtime,gate=gate,fallback=not gate,bound=bound,potential=pot,residual=rn,coercivity=bn,solve_seconds=solve,gate_seconds=gt,fallback_seconds=fbtime)
     full=out['DIRECT_LU'];uamp=np.array([(d.real-.8), (d.imag-.03)/.15]);qn=np.r_[x[:64],uamp];mn=prior.decode(qn)[0];assert prior.legal(qn);node_q.append(qn);node_material.append(mn)
     for arm in ARMS:
      raw[arm].append(out[arm]);mm=meta[arm];ae=0. if arm=='DIRECT_LU' else float(np.sqrt(2)/SIGMA*np.linalg.norm(approx[arm][-1]-full));row=dict(parent=parent,t=tval,node=ni,arm=arm,delta_real=float(d.real),delta_imag=float(d.imag),actual_approx_white_error=ae,returned_white_error=float(np.sqrt(2)/SIGMA*np.linalg.norm(out[arm]-full)),**mm);rows.append(row);b.event('NODE_RESULT',**row)
    pv.save(D/f'{parent}_t{ti}.npz',t=np.array(tval),latent_x=x,q_interpolated=q,interpolated_material=material,node_q=np.stack(node_q),node_material=np.stack(node_material),fraction=a,delta=deltas,B=B,L=L,R=R,Vnew=Vnew,V0=V0,**{arm:np.stack(raw[arm]) for arm in ARMS},**{arm+'_approx':np.stack(approx[arm]) for arm in ARMS[1:]})
   write(D/'PARTIAL.json',dict(completed_parents=[s['parent'] for s in anchors],completed_node_arm_rows=len(rows)))
  totals={}
  for arm in ARMS:
   rr=[r for r in rows if r['arm']==arm];common=sum(s['formation_seconds'] for s in setups);initial=sum(s['total_reuse_initial_seconds'] for s in anchors) if arm=='REUSE_ANCHOR_R64' else 0.;modelsetup=0. if arm=='DIRECT_LU' else sum(s['beta_Lnorm_seconds']+(s['rebuild_basis_seconds']+s['rebuild_model_seconds'] if arm=='REBUILD_R64' else s['reuse_model_seconds']) for s in setups);pure=sum(r['pure_seconds'] for r in rr);task=sum(r['total_seconds'] for r in rr);totals[arm]=dict(nodes=len(rr),gate_pass=sum(r['gate'] for r in rr),fallback_nodes=sum(r['fallback'] for r in rr),common_shape_setup_seconds=common,anchor_initial_basis_seconds=initial,own_update_model_setup_seconds=modelsetup,pure_online_seconds=pure,online_including_fallback_seconds=task,pure_approx_cold_task_seconds=common+initial+modelsetup+pure,full_cold_task_including_fallback_seconds=common+initial+modelsetup+task,max_white_bound=max((r['bound'] or 0.) for r in rr),max_actual_approx_white_error=max(r['actual_approx_white_error'] for r in rr))
  directtime=totals['DIRECT_LU']['full_cold_task_including_fallback_seconds']
  for arm in ARMS[1:]:totals[arm]['saving_full_task_vs_DIRECT']=1-totals[arm]['full_cold_task_including_fallback_seconds']/directtime
  write(D/'SUMMARY.json',dict(status='COMPLETE',scope=c['scope'],parents=8,shapes=24,nodes=192,totals=totals,shape_setups=setups,anchor_initial_setups=anchors,allocation='Each standalone arm pays full common newshape formation; both reduced arms pay full shared necessary beta/Lnorm once per shape in their standalone view. REUSE additionally pays actual coldanchor formation+basis onceperparent; REBUILD pays newbasis eachshape and needs no anchorbasis. Physically shared setup computedonce; studyallocatedformation/3 and beta/2 may be read from rawbreakdown but not usedto invent advantage.',GT_access=False,coverage_scoring=False,statistical_speedup_claim=False));status='COMPLETE'
 except BaseException as e:error=repr(e);write(D/'SUMMARY.json',dict(status='PARTIAL_STOP',error=error,completed_node_arm_rows=len(rows),new_compute_stopped=True,GT_access=False));write(D/'ERROR.json',dict(error=error,traceback=traceback.format_exc()));raise
 finally:
  if rows:
   keys=list(dict.fromkeys(k for row in rows for k in row))
   with (D/'NODE_RESULTS.csv').open('w',newline='') as h:w=csv.DictWriter(h,fieldnames=keys);w.writeheader();w.writerows(rows)
  b.close(status,error)
if __name__=='__main__':run()
