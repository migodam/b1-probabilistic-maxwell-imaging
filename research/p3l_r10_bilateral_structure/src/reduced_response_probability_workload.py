"""Frozen DEV conditional quadrature workload replay, no fitting, no material truth."""
from __future__ import annotations
import os
for n in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'):os.environ[n]='1'
import csv,hashlib,json,sys,time
from pathlib import Path
import numpy as np
from scipy.linalg import schur,eigh,solve_triangular
from scipy.special import logsumexp
P=Path(__file__).resolve().parents[1];sys.path.insert(0,str(P/'src'));import reduced_response_fresh as f;import reduced_response as dev
pv=f.pv;cg=f.cg;OUT=P/'results/reduced_response_probability_workload';CFG=P/'config/reduced_response_probability_workload_frozen.json';f.OUT=OUT

def run():
 c=json.loads(CFG.read_text())
 if c['source_blake2b']!=hashlib.blake2b(Path(__file__).read_bytes(),digest_size=32).hexdigest():raise RuntimeError('frozen source changed')
 b=f.Book(c);status='FAILED';error=None;rows=[]
 try:
  ts=time.perf_counter();a=pv.load_npz(cg.H/'results/stable_roi/SETUP.npz');cfg=json.loads((pv.R4/'research/p3l_r4_conditional_transport/configs/TARGET.json').read_text())
  with np.load(pv.R4/'research/p3l_r4_conditional_transport/inputs/background/weak.npz') as z:points=z['points'].copy()
  prior=pv.MaterialPrior(cfg,points);parent=c['parent'];ap=cg.H/'results/possibility_volume/parents'/parent/'anchor/ANCHOR.npz';x=pv.load_npz(ap)['x'];_,meta=prior.decode(pv.q_from_x(x));p=np.sqrt(meta['fraction'][a['ROI']]);B=p[:,None]*a['C']*p[None,:];L=a['O']*p[None,:];R=p[:,None]*a['ER'];Fh=a['Fh'];formation=time.perf_counter()-ts;y=pv.online_observation(dict(output=str(cg.H/'results/possibility_volume')),parent)
  ts=time.perf_counter();b.schur(256);T,Q=schur(B,output='complex',check_finite=False);LQ=L@Q;QR=Q.conj().T@R;fullsetup=time.perf_counter()-ts
  ts=time.perf_counter();beta=float(eigh((B-B.conj().T)/(2j),eigvals_only=True,subset_by_index=(0,0),check_finite=False)[0]);Lnorm=float(np.linalg.norm(L,'fro'));energy=time.perf_counter()-ts
  ts=time.perf_counter();V,widths=dev.basis(B,R,64,c['rank_relative_cutoff'],b);basis_time=time.perf_counter()-ts;rank=V.shape[1];model=None;model_time=0.
  if rank>=64:
   ts=time.perf_counter();b.matvec+=64;BV=B@V;Hr=V.conj().T@BV;b.schur(64);Tr,Qr=schur(Hr,output='complex',check_finite=False);VQ=V@Qr;BVQ=BV@Qr;LVQ=L@VQ;RQ=Qr.conj().T@(V.conj().T@R);model=(Tr,VQ,BVQ,LVQ,RQ);model_time=time.perf_counter()-ts
  raw=list(csv.DictReader((P/c['node_source']).open()));ns=[r for r in raw if r['method']==c['node_method'] and int(r['level'])==c['node_level']];assert len(ns)==836;v=np.array([[float(r['delta_real']),float(r['delta_imag'])] for r in ns]);weights=np.array([float(r['integration_weight']) for r in ns]);ds=v[:,0]+1j*v[:,1];Fs=[];As=[];Ys=[];phi_full=[];phi_red=[];I=np.eye(256,dtype=complex)
  def exact(d,method):
   b.reserve(256,method);t=time.perf_counter();u=solve_triangular(I-d*T,QR,check_finite=False);F=Fh+d*(LQ@u);dt=time.perf_counter()-t;b.event('SOLVE_COMPLETE',method=method,order=256,wall_seconds=dt);return F,dt
  for i,d in enumerate(ds):
   b.context=dict(parent=parent,node=i);out={};tim={};gate=False;bound=potential=None;rs=gt=ft=0.;approx=Fh.copy();returned=Fh.copy();fallback=False
   order=['FULL','REDUCED'] if i%2==0 else ['REDUCED','FULL']
   for method in order:
    if method=='FULL':out['FULL'],tim['FULL']=exact(d,'EXACT_COMPLEX_SCHUR256')
    else:
     if model is not None:
      Tr,VQ,BVQ,LVQ,RQ=model;b.reserve(64,'KRYLOV64');t=time.perf_counter();u=solve_triangular(np.eye(64)/d-Tr,RQ,check_finite=False);approx=Fh+LVQ@u;rs=time.perf_counter()-t;b.event('SOLVE_COMPLETE',method='KRYLOV64',order=64,wall_seconds=rs)
      t=time.perf_counter();res=R-(VQ@u)/d+BVQ@u;bn=float(d.imag/abs(d)**2+beta)
      if bn>0:bound=float(np.sqrt(2)/f.SIGMA*Lnorm*np.linalg.norm(res)/bn);potential=float(np.sqrt(2)/f.SIGMA*np.linalg.norm(approx-y)*bound+.5*bound**2);gate=bound<=c['white_bound_gate'] and potential<=c['potential_bound_gate']
      gt=time.perf_counter()-t
     if gate:returned=approx
     else:fallback=True;returned,ft=exact(d,'EXACT_FALLBACK_SCHUR256')
   full=out['FULL'];pf=float(np.linalg.norm(full-y)**2/f.SIGMA**2);pr=float(np.linalg.norm(returned-y)**2/f.SIGMA**2);err=float(np.sqrt(2)/f.SIGMA*np.linalg.norm(approx-full));phi_full.append(pf);phi_red.append(pr);Fs.append(full);As.append(approx);Ys.append(returned);rows.append(dict(node=i,delta_real=float(d.real),delta_imag=float(d.imag),integration_weight=float(weights[i]),gate=bool(gate),fallback=fallback,white_bound=bound,potential_bound=potential,actual_white_error=err,bound_dominates_sampled_error=bool(bound is not None and err<=bound*(1+1e-8)+1e-10),full_Phi=pf,delivered_Phi=pr,full_seconds=tim['FULL'],reduced_solve_seconds=rs,gate_seconds=gt,fallback_seconds=ft,reduced_seconds=rs+gt+ft))
  def moments(phis):
   t=time.perf_counter();lw=np.log(weights)-np.asarray(phis)-np.log(.15);logz=float(logsumexp(lw));w=np.exp(lw-logz);mean=w@v;var=w@((v-mean)**2);events=w@np.column_stack((v[:,0]>1.3,v[:,1]>.105));result=dict(log_evidence=logz,amplitude_mean=mean.tolist(),amplitude_variance=var.tolist(),event_probability=events.tolist());return result,time.perf_counter()-t
  mf,qf=moments(phi_full);mr,qr=moments(phi_red);diff=dict(relative_evidence=float(abs(np.expm1(mr['log_evidence']-mf['log_evidence']))),amplitude_mean=float(np.max(abs(np.array(mr['amplitude_mean'])-mf['amplitude_mean']))),event_probability=float(np.max(abs(np.array(mr['event_probability'])-mf['event_probability']))));tol=c['tolerances'];passed=all(diff[k]<=tol[k] for k in tol);fallback_setup=fullsetup if any(r['fallback'] for r in rows) else 0.;fulltotal=formation+fullsetup+sum(r['full_seconds'] for r in rows)+qf;redtotal=formation+energy+basis_time+model_time+fallback_setup+sum(r['reduced_seconds'] for r in rows)+qr
  result=dict(protocol_id=c['id'],status='NUMERICAL_FIXED_SHAPE_TASK_REPLAY' if passed else 'NUMERICAL_TASK_REPLAY_UNQUALIFIED',parent=parent,nodes=836,rank=int(rank),gate_pass=sum(r['gate'] for r in rows),fallback_nodes=sum(r['fallback'] for r in rows),max_white_bound=max(r['white_bound'] or 0 for r in rows),max_actual_white_error=max(r['actual_white_error'] for r in rows),max_potential_bound=max(r['potential_bound'] or 0 for r in rows),sampled_bounds_all_hold=all(r['bound_dominates_sampled_error'] for r in rows),exact_moments=mf,reduced_moments=mr,differences=diff,tolerances=tol,numerical_task_tolerances_pass=passed,cost=dict(shared_formation_and_dependency_load_seconds=formation,full_complexSchur_setup_seconds=fullsetup,beta_and_Lnorm_seconds=energy,basis_to64_seconds=basis_time,reduced_Schur_and_ports_seconds=model_time,exact_quadrature_seconds=qf,reduced_quadrature_seconds=qr,exact_cold_workload_seconds=fulltotal,reduced_cold_workload_seconds=redtotal,saving=1-redtotal/fulltotal),historical_profile_and_grid='sharedexistingREV2mode/Hessian/eventsplit836nodes/weights; no refitting; historicalcomplete transfer69784RHS/270LU includes this construction but no separatelymeasured puregridtime; not invented or subtracted',GT_access=False,independent_shapes=False,fullsupport_certificate=False,scope='DEVparent000 fixedshape conditional numericaltask replay at preexistingrule; nodewiseresidualbound, notfullsupport or Full66posterior')
  pv.write(OUT/'RESULTS.json',result);f.csvwrite('NODE_RESULTS.csv',rows);pv.save(OUT/'RESPONSES.npz',delta=ds,weights=weights,full_response=np.stack(Fs),approx_response=np.stack(As),delivered_response=np.stack(Ys),full_Phi=np.array(phi_full),reduced_Phi=np.array(phi_red),B=B,L=L,R=R,V64=V);status='COMPLETE';print(json.dumps(result,indent=2))
 except BaseException as e:error=repr(e);pv.write(OUT/'FAILURE.json',dict(error=error,completed_nodes=len(rows)));raise
 finally:b.close(status,error)
if __name__=='__main__':run()
