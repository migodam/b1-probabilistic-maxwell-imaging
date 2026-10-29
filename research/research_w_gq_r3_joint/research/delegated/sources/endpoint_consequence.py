"""Offline dev-only same-state paid Jacobian consequence analysis; no physics imports."""
from pathlib import Path
import json,time,os
import numpy as np
ROOT=Path(__file__).resolve().parents[3]
DEV=ROOT/'results/online/dev_r0'
OUT=Path(__file__).resolve().parent/'ENDPOINT_CONSEQUENCE_RUN'
def dump(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def rel(a,b):
 n=np.linalg.norm(b)
 return float(np.linalg.norm(a)/n) if n else None

def main():
 OUT.mkdir(exist_ok=False);t=time.perf_counter();cpu=time.process_time();rows=[]
 for method in ['FIXED_CURRENT_CHART','ADAPTIVE_CURRENT_CHART']:
  directory=DEV/method;contract=json.loads((directory/'ONLINE_INPUT_CONTRACT.json').read_text());spec=contract['spec'];sigma=float(contract['config']['sigma']);scale=np.asarray(spec['scale']);lo=np.asarray(spec['lower']);hi=np.asarray(spec['upper']);endpoint=np.load(directory/'endpoint_exact_jacobian.npz');y=np.load(ROOT/'inputs/online/dev_r0/observation.npz')['y']
  sources=[directory/'endpoint_exact_jacobian.npz']+sorted(directory.glob('joint_start_*_FINAL_JACOBIAN.npz'))+sorted((directory/'jacobian_audits').glob('CURRENT_AUDIT_*.npz'))
  for path in sources:
   d=np.load(path);z=d['z'];jc=d['exact_jacobian'];J=np.vstack((jc.real,jc.imag))*scale/sigma;data=d['data'];rc=data-y;r=np.r_[rc.real.ravel(),rc.imag.ravel()]/sigma
   rp=np.r_[(z[:6]-3)/1.5,(z[6:12]-.06)/.15];jp=np.zeros((12,19));jp[np.arange(12),np.arange(12)]=np.r_[np.full(6,1/1.5),np.full(6,1/.15)]*scale[:12];precision=jp.T@jp;g=J.T@r+jp.T@rp;diag=np.maximum(np.sqrt(np.sum(J*J,axis=0)+np.sum(jp*jp,axis=0)),1.)
   def pg(gradient):
    v=gradient.copy();v[(z<=lo+1e-9)&(v>0)]=0;v[(z>=hi-1e-9)&(v<0)]=0
    return float(np.max(np.abs(v)/diag))
   row=dict(method=method,source=str(path.relative_to(ROOT)),shared_state_within_file=True,exact_projected_stationarity=pg(g),same_full_response_objective=float(.5*(r@r+rp@rp)),sigma=sigma,scaled_coordinates=True)
   approx=d['approximate_jacobian'] if 'approximate_jacobian' in d else None
   if approx is None and path.name=='endpoint_exact_jacobian.npz':
    for candidate in sorted(directory.glob('joint_start_*_FINAL_JACOBIAN.npz')):
     other=np.load(candidate)
     if 'approximate_jacobian' in other and np.array_equal(other['z'],z) and np.array_equal(other['data'],data):approx=other['approximate_jacobian'];row['matched_approximate_source']=str(candidate.relative_to(ROOT));break
   if approx is None:row['chart_pair_status']='NOT_SAVED';rows.append(row);continue
   row['chart_pair_status']='SAME_STATE_PAIR';Jh=np.vstack((approx.real,approx.imag))*scale/sigma;E=Jh-J;gh=Jh.T@r+jp.T@rp;row['chart_projected_stationarity_same_full_diagonal']=pg(gh)
   row['tangent_errors']={k:dict(relative_frobenius=rel(E[:,s],J[:,s]),absolute_frobenius=float(np.linalg.norm(E[:,s]))) for k,s in [('whole',slice(None)),('wall',slice(12)),('target',slice(12,None))]}
   row['gradient_miss_norm']={k:float(np.linalg.norm((E.T@r)[s])) for k,s in [('whole',slice(None)),('wall',slice(12)),('target',slice(12,None))]};lamw=precision[:12,:12]
   def schur(j,damping):
    a,b=j[:,:12],j[:,12:];hw=a.T@a+lamw+damping*np.eye(12);c=a.T@b
    return b.T@b+damping*np.eye(7)-c.T@np.linalg.solve(hw,c),hw,c
   sd,hw0,c0=schur(J,0);sdh,_,_=schur(Jh,0);ss,hw,c=schur(J,1e-3);ssh,_,_=schur(Jh,1e-3)
   vals=np.linalg.eigvalsh((sd+sd.T)/2);row['dataSchur']=dict(relative_frobenius_error=rel(sdh-sd,sd),absolute_operator_error=float(np.linalg.norm(sdh-sd,2)),full_eigenvalues=vals.tolist(),target_prior_precision=0,algorithm_damping=0)
   N=np.linalg.solve(hw0,c0);R=J[:,12:]-J[:,:12]@N;target_effective=E[:,12:]-E[:,:12]@N;leading=target_effective.T@R+R.T@target_effective;actual=sdh-sd
   row['first_order_target_effective']=dict(convention='DeltaJ=chart-full; derivative evaluated at full J',effective_tangent_norm=float(np.linalg.norm(target_effective)),relative_to_profiled_response=rel(target_effective,R),leading_Schur_norm=float(np.linalg.norm(leading)),actual_Schur_norm=float(np.linalg.norm(actual)),relative_linearization_remainder=rel(actual-leading,actual),remainder_relative_to_fullSchur=rel(actual-leading,sd),finite_certificate=False)
   H=J.T@J+precision+1e-3*np.eye(19);Hh=Jh.T@Jh+precision+1e-3*np.eye(19);delta=-np.linalg.solve(H,g);dh=-np.linalg.solve(Hh,gh);e=H@dh+g;ep=e[12:]-c.T@np.linalg.solve(hw,e[:12]);target_pred=np.linalg.solve(ss,ep)
   row['damped_unconstrained_step']=dict(damping=1e-3,damping_scope='ALL19 scaled coordinates; not posterior',relative_whole_difference=rel(dh-delta,delta),relative_target_difference=rel(dh[12:]-delta[12:],delta[12:]),target_difference_norm=float(np.linalg.norm(dh[12:]-delta[12:])),profiled_residual_norm=float(np.linalg.norm(ep)),F29_target_prediction_norm=float(np.linalg.norm(target_pred)),full_stepSchur_min_eigenvalue=float(np.linalg.eigvalsh((ss+ss.T)/2)[0]),stepSchur_relative_frobenius_error=rel(ssh-ss,ss))
   arrays=dict(z=z,parameter_scales=scale,sigma=np.asarray(sigma),J_full_scaled_white=J,J_chart_scaled_white=Jh,residual_white=r,gradient_full=g,gradient_chart=gh,prior_precision_scaled=precision,dataSchur_full=sd,dataSchur_chart=sdh,stepSchur_full=ss,stepSchur_chart=ssh,step_full=delta,step_chart=dh,normal_residual=e,profiled_residual=ep,F29_target_error=target_pred,wall_compensation_N=N,profiled_response_R=R,target_effective_tangent=target_effective,Schur_leading_term=leading)
   if 'proposed_direction' in d:
    v=d['proposed_direction']/scale;row['paid_proposed_direction']=dict(scaled_norm=float(np.linalg.norm(v)),full_white_change_norm=float(np.linalg.norm(J@v)),white_defect_norm=float(np.linalg.norm(E@v)),relative_defect=rel(E@v,J@v));arrays['proposed_direction_scaled']=v
   target=OUT/(method+'__'+path.stem+'.npz');np.savez_compressed(target,**arrays);row['derived_arrays']=str(target.relative_to(ROOT));rows.append(row)
 dump(OUT/'SUMMARY.json',dict(scope='Already-paid dev same-state arrays only; no PDE, truth or tuning',rows=rows,limits='Local derivative/QP consequences, not truth accuracy or global identifiability'))
 receipt=dict(physical_cost=dict(full_fe=0,rhs=0),process_cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-t,rows=len(rows),same_state_pairs=sum(r.get('chart_pair_status')=='SAME_STATE_PAIR' for r in rows),pid=os.getpid());dump(OUT/'RECEIPT.json',receipt);print(json.dumps(receipt))
if __name__=='__main__':main()
