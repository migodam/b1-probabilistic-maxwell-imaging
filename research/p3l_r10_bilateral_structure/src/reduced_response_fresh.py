"""Frozen heldout-shape fixed64 Krylov kernel study; no material truth read."""
from __future__ import annotations
import os
for n in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'):os.environ[n]='1'
import argparse,csv,fcntl,hashlib,json,sys,time,resource,traceback
from pathlib import Path
import numpy as np
from scipy.linalg import schur,lu_factor,lu_solve,solve_triangular,eigh
P=Path(__file__).resolve().parents[1];CFG=P/'config/reduced_response_fresh10_frozen.json';OUT=P/'results/reduced_response_fresh10';SIGMA=.00119211798514519
sys.path.insert(0,str(P/'src'));import candidate_generation as cg;import reduced_response as dev
pv=cg.pv
class Stop(RuntimeError):pass

def write(p,x):pv.write(p,x)
def csvwrite(name,rows):
 if not rows:return
 with (OUT/name).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for row in rows for k in row)));w.writeheader();w.writerows(rows)
class Book:
 def __init__(self,c):
  self.c=c;OUT.mkdir(parents=True,exist_ok=True);self.f=(OUT/'COST_LEDGER.jsonl').open('x');self.lock=(P/'candidate_physics_worker.lock').open('a');fcntl.flock(self.lock,fcntl.LOCK_EX|fcntl.LOCK_NB);self.start=time.perf_counter();self.cpu=time.process_time();self.seq=0;self.rhs={};self.lu=0;self.matvec=0;self.schurs={};self.context={};self.failures=0;self.event('BEGIN',pid=os.getpid(),protocol=c['id'])
 def event(self,event,**kw):
  self.seq+=1;self.f.write(json.dumps(pv.plain(dict(sequence=self.seq,utc=time.time(),event=event,context=self.context.copy(),**kw)),allow_nan=False)+'\n');self.f.flush()
 def check(self):
  if time.perf_counter()-self.start>self.c['caps']['wall_seconds']:raise Stop('controlledwallcap')
 def reserve(self,order,method,lu=0):
  self.check()
  if sum(self.rhs.values())+8>self.c['caps']['complex_RHS_all_orders'] or self.lu+lu>self.c['caps']['LU256']:raise Stop('physicalhardcap')
  self.rhs[order]=self.rhs.get(order,0)+8;self.lu+=lu;self.event('SOLVE_RESERVED',order=order,method=method,RHS=8,LU256=lu,after_RHS_all_orders=sum(self.rhs.values()),after_LU256=self.lu)
 def schur(self,order):self.schurs[order]=self.schurs.get(order,0)+1;self.event('SCHUR_ATTEMPT',order=order)
 def close(self,status,error=None):
  d=dict(status=status,error=error,RHS_by_order=self.rhs,RHS_conservative_total=sum(self.rhs.values()),RHS256=self.rhs.get(256,0),LU256=self.lu,RHS576=0,Schur_decompositions_by_order=self.schurs,Green_matrix_columns_applied=self.matvec,failed_actions=self.failures,wall_seconds=time.perf_counter()-self.start,CPU_seconds=time.process_time()-self.cpu,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
  self.event('END',**d);write(OUT/'RECEIPT.json',d);self.f.close();self.lock.close()

def run():
 c=json.loads(CFG.read_text());expected=c['source_blake2b'];actual=hashlib.blake2b(Path(__file__).read_bytes(),digest_size=32).hexdigest()
 if expected!=actual:raise RuntimeError('sourcechangedafterfreeze')
 b=Book(c);rows=[];setups=[];status='FAILED';error=None
 try:
  t=time.perf_counter();a=pv.load_npz(cg.H/'results/stable_roi/SETUP.npz');cfg=json.loads((pv.R4/'research/p3l_r4_conditional_transport/configs/TARGET.json').read_text())
  with np.load(pv.R4/'research/p3l_r4_conditional_transport/inputs/background/weak.npz') as z:points=z['points'].copy()
  prior=pv.MaterialPrior(cfg,points);idx=a['ROI'].astype(int);loadtime=time.perf_counter()-t;b.event('COMMON_HOST_CACHE_LOAD',wall_seconds=loadtime,historical_RHS576=264,new_RHS576=0)
  I256=np.eye(256,dtype=complex)
  for pi,parent in enumerate(c['parents']):
   b.check();b.context=dict(parent=parent,stage='COLD_SETUP');t=time.perf_counter()
   with np.load(P/'results/candidates_fresh10/parents'/parent/'anchor/ANCHOR.npz') as z:x=z['x'].copy()
   # Decode only fixed shape from a legitimately acquired observation-only anchor.
   _,meta=prior.decode(pv.q_from_x(x));fraction=meta['fraction'][idx];p=np.sqrt(fraction);B=p[:,None]*a['C']*p[None,:];L=a['O']*p[None,:];R=p[:,None]*a['ER'];formation=time.perf_counter()-t;y=pv.online_observation(dict(output=str(P/'results/candidates_fresh10')),parent);Fh=a['Fh'];allhost=bool(not np.any(fraction))
   rng=np.random.default_rng(np.random.SeedSequence(c['amplitude_seed_namespace']+[pi]));U=rng.uniform(size=(32,2));deltas=.8+U[:,0]+1j*(.03+.15*U[:,1]);rank=0;beta=0.;Lnorm=0.;widths=[];basis_time=energy_time=model_time=fullschur_time=0.;V=np.zeros((256,0),complex);model=None;T=Q=LQ=QR=None
   if not allhost:
    t=time.perf_counter();b.schur(256);T,Q=schur(B,output='complex',check_finite=False);LQ=L@Q;QR=Q.conj().T@R;fullschur_time=time.perf_counter()-t
    t=time.perf_counter();beta=float(eigh((B-B.conj().T)/(2j),eigvals_only=True,subset_by_index=(0,0),check_finite=False)[0]);Lnorm=float(np.linalg.norm(L,'fro'));energy_time=time.perf_counter()-t
    t=time.perf_counter();V,widths=dev.basis(B,R,64,c['rank_relative_cutoff'],b);basis_time=time.perf_counter()-t;rank=V.shape[1]
    if rank>=64:
     t=time.perf_counter();Vr=V[:,:64];b.matvec+=64;BV=B@Vr;Hr=Vr.conj().T@BV;b.schur(64);Tr,Qr=schur(Hr,output='complex',check_finite=False);VQ=Vr@Qr;BVQ=BV@Qr;LVQ=L@VQ;RQ=Qr.conj().T@(Vr.conj().T@R);model=(Tr,VQ,BVQ,LVQ,RQ);model_time=time.perf_counter()-t
   setup=dict(parent=parent,formation_seconds=formation,full_Schur_setup_seconds=fullschur_time,coercivity_and_Lnorm_seconds=energy_time,basis_to64_seconds=basis_time,reduced_Schur_and_port_transforms_seconds=model_time,beta=beta,bmin_prior_rectangle=.03/(1.8**2+.03**2)+beta,L_Frobenius=Lnorm,basis_rank=int(rank),basis_block_widths=widths,allhost_analytic=allhost,reduced_available=bool(model is not None),min_fraction=float(fraction.min()),max_fraction=float(fraction.max()));setups.append(setup);b.event('COLD_SETUP_COMPLETE',**setup)
   exacts=[];approxs=[];delivered=[];records=[]
   def exact(delta,method):
    if allhost:return Fh.copy(),0.
    b.reserve(256,method,int(method=='EXACT_LU256'));ts=time.perf_counter()
    try:
     if method=='EXACT_LU256':fac=lu_factor(I256-delta*B,check_finite=False);u=lu_solve(fac,R,check_finite=False);f=Fh+delta*(L@u)
     else:u=solve_triangular(I256-delta*T,QR,check_finite=False);f=Fh+delta*(LQ@u)
     elapsed=time.perf_counter()-ts;b.event('SOLVE_COMPLETE',order=256,method=method,wall_seconds=elapsed);return f,elapsed
    except BaseException as e:b.failures+=1;b.event('SOLVE_FAILURE',order=256,method=method,error=repr(e));raise
   for ni,delta in enumerate(deltas):
    b.context=dict(parent=parent,node=ni);order=c['methods'][(pi+ni)%3:]+c['methods'][:(pi+ni)%3];out={};times={};gate=False;bound=None;potential=None;bn=None;resnorm=None;rsolve=gatetime=fallbacktime=0.;fallback=False;approx=Fh.copy();returned=Fh.copy()
    try:
     for method in order:
      if method.startswith('EXACT_'):out[method],times[method]=exact(delta,method)
      else:
       if allhost:gate=True;bound=potential=0.;bn=float(delta.imag/abs(delta)**2);resnorm=0.
       elif model is not None:
        Tr,VQ,BVQ,LVQ,RQ=model;b.reserve(64,'KRYLOV64');ts=time.perf_counter();u=solve_triangular(np.eye(64)/delta-Tr,RQ,check_finite=False);approx=Fh+LVQ@u;rsolve=time.perf_counter()-ts;b.event('SOLVE_COMPLETE',order=64,method='KRYLOV64',wall_seconds=rsolve)
        ts=time.perf_counter();res=R-(VQ@u)/delta+BVQ@u;resnorm=float(np.linalg.norm(res));bn=float(delta.imag/abs(delta)**2+beta)
        if bn>0:
         bound=float(np.sqrt(2)/SIGMA*Lnorm*resnorm/bn);rw=float(np.sqrt(2)/SIGMA*np.linalg.norm(approx-y));potential=rw*bound+.5*bound*bound;gate=bool(np.isfinite(bound) and bound<=c['white_bound_gate'] and potential<=c['potential_bound_gate'])
        gatetime=time.perf_counter()-ts
       if gate:returned=approx
       else:fallback=True;returned,fallbacktime=exact(delta,'EXACT_FALLBACK_SCHUR256')
    except Stop as e:
     rows.append(dict(parent=parent,node=ni,status='RIGHT_CENSORED_CAP',gate=False,error=str(e)));raise
    fe=out['EXACT_COMPLEX_SCHUR256'];fd=out['EXACT_LU256'];err=float(np.sqrt(2)/SIGMA*np.linalg.norm(approx-fe));ret_err=float(np.sqrt(2)/SIGMA*np.linalg.norm(returned-fe));directerr=float(np.sqrt(2)/SIGMA*np.linalg.norm(fd-fe));row=dict(parent=parent,node=ni,status='COMPLETE',gate=gate,error=None,rank=int(rank),allhost_analytic=allhost,delta_real=float(delta.real),delta_imag=float(delta.imag),coercivity=bn,residual_Fnorm=resnorm,white_bound=bound,potential_bound=potential,actual_approx_white_error=err,actual_returned_white_error=ret_err,bound_dominates_sampled_error=bool(bound is not None and err<=bound*(1+1e-8)+1e-10),direct_vs_Schur_white_error=directerr,exact_fallback=fallback,LU256_seconds=times['EXACT_LU256'],full_Schur_seconds=times['EXACT_COMPLEX_SCHUR256'],rank64_solve_seconds=rsolve,full_residual_and_gate_seconds=gatetime,fallback_seconds=fallbacktime,reduced_delivered_seconds=rsolve+gatetime+fallbacktime,execution_order=','.join(order));rows.append(row);records.append(row);exacts.append(fe);approxs.append(approx);delivered.append(returned)
   pv.save(OUT/f'{parent}.npz',x_shape=x[:64],fraction=fraction,delta=deltas,full_response=np.stack(exacts),approx_response=np.stack(approxs),delivered_response=np.stack(delivered),V64=V,B=B,L=L,R=R,observation=y)
   write(OUT/f'{parent}_SETUP.json',setup);print('HELDOUT_SHAPE',parent,'gate',sum(r['gate'] for r in records),'/32','maxbound',max(r['white_bound'] or 0 for r in records),flush=True)
  csvwrite('NODE_RESULTS.csv',rows);csvwrite('SETUP_COSTS.csv',setups);cold=[]
  for setup in setups:
   rs=[r for r in rows if r['parent']==setup['parent']];fallback=any(r['exact_fallback'] for r in rs);redsetup=setup['formation_seconds']+setup['coercivity_and_Lnorm_seconds']+setup['basis_to64_seconds']+setup['reduced_Schur_and_port_transforms_seconds']+(setup['full_Schur_setup_seconds'] if fallback else 0);fullsetup=setup['formation_seconds']+setup['full_Schur_setup_seconds'];directsetup=setup['formation_seconds']
   for M in c['M']:
    rr=redsetup+M*np.mean([r['reduced_delivered_seconds'] for r in rs]);ff=fullsetup+M*np.mean([r['full_Schur_seconds'] for r in rs]);dd=directsetup+M*np.mean([r['LU256_seconds'] for r in rs]);cold.append(dict(parent=setup['parent'],M=M,reduced_total_seconds=float(rr),full_Schur_total_seconds=float(ff),direct_LU_total_seconds=float(dd),saving_vs_full_Schur=float(1-rr/ff),saving_vs_direct_LU=float(1-rr/dd),actual_measured_query_count=32,scope='M32measuredsum; otherMmeasuredmeanprojection; commonhostcacheIOseparate'))
  csvwrite('COLD_AMORTIZATION.csv',cold);pooled=[]
  for M in c['M']:
   rr=[r for r in cold if r['M']==M];red=sum(r['reduced_total_seconds'] for r in rr);full=sum(r['full_Schur_total_seconds'] for r in rr);direct=sum(r['direct_LU_total_seconds'] for r in rr);pooled.append(dict(M=M,reduced_seconds=red,full_Schur_seconds=full,direct_LU_seconds=direct,saving_vs_full_Schur=1-red/full,saving_vs_direct_LU=1-red/direct))
  write(OUT/'RESULTS.json',dict(protocol_id=c['id'],status='COMPLETE',nodes=len(rows),parents=len(setups),gate_pass=sum(r['gate'] for r in rows),exact_fallback_nodes=sum(r['exact_fallback'] for r in rows),max_actual_approx_white_error=max(r['actual_approx_white_error'] for r in rows),max_actual_returned_white_error=max(r['actual_returned_white_error'] for r in rows),max_direct_vs_Schur_white_error=max(r['direct_vs_Schur_white_error'] for r in rows),max_white_bound=max(r['white_bound'] or 0 for r in rows),max_potential_bound=max(r['potential_bound'] or 0 for r in rows),sampled_bounds_all_hold=all(r['bound_dominates_sampled_error'] for r in rows),cold_pooled=pooled,common_host_cache_load_seconds=loadtime,setup_details=setups,GT_access=False,scope=c['source_scope'],qualification=c['qualification'],method='standardblockKrylov/Galerkin+scalaramplitudereuse; no newtheorem; onlyr64selectedfromDEVthenheldoutshapetested'))
  status='COMPLETE'
 except BaseException as e:error=repr(e);write(OUT/'FAILURE.json',dict(error=error,traceback=traceback.format_exc(),completed_rows=len(rows)));csvwrite('PARTIAL_NODE_RESULTS.csv',rows);raise
 finally:b.close(status,error)
if __name__=='__main__':run()
