#!/usr/bin/env python3
"""Paid exact Hessian curvature experiment; no truth loaded in this program."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
import json,time,importlib.util
from pathlib import Path
import numpy as np
from scipy.linalg import lu_solve
from scipy.sparse.linalg import LinearOperator,eigsh,ArpackNoConvergence
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parents[1];OLD=ROOT.parent/'research_z_alg_r3_digits'
def module(name,path):
 spec=importlib.util.spec_from_file_location(name,path);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod
physics=module('old_physics_hessian',OLD/'src/physics.py');multiscale=module('old_regularizer_hessian',OLD/'src/multiscale.py')
def clean(a):
 if isinstance(a,np.ndarray):return clean(a.tolist())
 if isinstance(a,np.generic):return clean(a.item())
 if isinstance(a,dict):return {k:clean(v) for k,v in a.items()}
 if isinstance(a,(list,tuple)):return [clean(v) for v in a]
 if isinstance(a,float) and not np.isfinite(a):return None
 return a
def save(results):
 (ROOT/'diagnostics/hessian_geometry.json').write_text(json.dumps(clean(results),indent=2,allow_nan=False)+'\n')
class MatvecCap(RuntimeError):pass

def run_case(name,source,oid):
 begun,cpu=time.perf_counter(),time.process_time();frozen=json.loads(source.read_text());m=np.asarray(frozen['top1']['m']);cfg=frozen['config']
 public=dict(np.load(OLD/'public_data/operators.npz'));target=np.load(OLD/f'public_data/holdout/{oid}.npz')['Y']
 model=physics.Model(public);reg=multiscale.Regularizer(public['grid_shape'],cfg['smooth_weight'],cfg['tv_weight'],cfg['tv_epsilon'])
 state=model.forward(m,label=name+'_exact_H_base');Q=state['E'];R=state['y']-target;scale=float(np.vdot(target,target).real);basefactor=model.factor
 with model.operation('HESSIAN_BASE_ADJOINT',rhs=model.ns,adjoint_RHS=model.ns):
  adj=lu_solve(basefactor,model.SH@R,trans=2,check_finite=False)
 with model.operation('EXACT_REGULARIZER_HESSIAN_SETUP',exact_reg_setup=1):
  regvalue,reggrad,_=reg.evaluate(m);d=reg.D@m;q=np.sqrt(d*d+reg.epsilon**2)
  w=(2*reg.smooth+reg.tv*reg.epsilon**2/q**3)/reg.edges
 g=np.real(np.sum(adj.conj()*Q,axis=1))/scale+reggrad
 lower,upper=float(cfg['lower']),float(cfg['upper']);lo=m<=lower+1e-6;hi=m>=upper-1e-6;free=~(lo|hi)
 pg=g.copy();pg[lo & (g>0)]=0;pg[hi & (g<0)]=0
 mapping=m-np.clip(m-g,lower,upper)
 indices=np.where(free)[0];calls=0
 def hvp(h):
  with model.operation('EXACT_HESSIAN_VECTOR',rhs=2*model.ns,HVP=1,tangent_RHS=model.ns,incremental_adjoint_RHS=model.ns):
   U=lu_solve(basefactor,h[:,None]*Q,check_finite=False)
   dQ=model.G@U
   rhs=model.GH@(h[:,None]*adj)+model.SH@(model.S@U)
   dadj=lu_solve(basefactor,rhs,trans=2,check_finite=False)
   out=np.real(np.sum(dadj.conj()*Q+adj.conj()*dQ,axis=1))/scale
  with model.operation('EXACT_REGULARIZER_HESSIAN_VECTOR',reg_HVP=1):out+=reg.D.T@(w*(reg.D@h))
  return out
 def mv(v):
  nonlocal calls
  if calls>=250:raise MatvecCap('MAX_250_MATVEC')
  calls+=1;h=np.zeros(model.n);h[indices]=v;return hvp(h)[indices]
 result=dict(case=name,status='OFFLINE_PUBLIC_Y_CURVATURE_EXPERIMENT_NO_GT',source=str(source.relative_to(ROOT.parent)),truth_access='NONE',
  relative_residual=float(np.linalg.norm(R)/np.sqrt(scale)),objective=float(.5*np.vdot(R,R).real/scale+regvalue),material=m,
  KKT=dict(projected_gradient_inf=float(np.linalg.norm(pg,np.inf)),projected_gradient_l2=float(np.linalg.norm(pg)),unit_projected_mapping_inf=float(np.linalg.norm(mapping,np.inf)),raw_gradient_inf=float(np.linalg.norm(g,np.inf)),free_count=int(free.sum()),active_lower=int(lo.sum()),active_upper=int(hi.sum()),lower_KKT_violations=int(np.count_nonzero(lo&(g < -1e-8))),upper_KKT_violations=int(np.count_nonzero(hi&(g>1e-8))),tolerance_for_active=1e-6,stationarity_acceptance_threshold='NOT_INVENTED; raw metrics reported'),
  eigensolver=dict(method='scipy eigsh smallest algebraic',k=min(3,max(len(indices)-1,0)),tol=1e-3,maxmatvec=250,maxiter=250,free_indices=indices,symmetry='Exact analytic symmetric Hessian action; no numerical formula/symmetry accuracy test added'))
 eigbeg,eigcpu=time.perf_counter(),time.process_time();values=np.empty(0);vectors=np.empty((len(indices),0))
 try:
  if len(indices)>=2:
   values,vectors=eigsh(LinearOperator((len(indices),len(indices)),matvec=mv,dtype=float),k=result['eigensolver']['k'],which='SA',tol=1e-3,maxiter=250,ncv=min(len(indices),20),v0=np.cos(np.arange(len(indices))+.314159))
   result['eigensolver']['status']='CONVERGED_TOL_1E_MINUS3'
  else:result['eigensolver']['status']='INSUFFICIENT_FREE_DIMENSION'
 except ArpackNoConvergence as exc:
  values,vectors=exc.eigenvalues,exc.eigenvectors;result['eigensolver']['status']='ARPACK_NO_CONVERGENCE_PARTIAL_PAID';result['eigensolver']['error']=str(exc)
 except MatvecCap as exc:result['eigensolver']['status']='MAX_250_MATVEC_PAID';result['eigensolver']['error']=str(exc)
 result['eigensolver'].update(matvec_count=calls,wall_seconds=time.perf_counter()-eigbeg,cpu_seconds=time.process_time()-eigcpu,eigenvalues=values)
 curves=[];directions=[]
 for i,val in enumerate(values):
  h=np.zeros(model.n);h[indices]=vectors[:,i]
  with model.operation('GN_RAYLEIGH_TANGENT',rhs=model.ns,GN_Rayleigh=1,tangent_RHS=model.ns):U=lu_solve(basefactor,h[:,None]*Q,check_finite=False);SU=model.S@U
  data_gn=float(np.vdot(SU,SU).real/scale);reg_r=float(np.dot(reg.D@h,w*(reg.D@h)))
  directions.append(dict(index=i,exact_eigenvalue=float(val),GN_data_Rayleigh=data_gn,exact_reg_Rayleigh=reg_r,GN_plus_exact_reg_Rayleigh=data_gn+reg_r,residual_Hessian_Rayleigh=float(val-data_gn-reg_r),direction_unit_l2=h,gradient_dot_direction=float(g@h)))
 result['directions']=directions
 print(json.dumps(clean(dict(case=name,KKT=result['KKT'],eigenvalues=values,eigensolver_status=result['eigensolver']['status'],matvecs=calls))),flush=True)
 if len(values) and values[0]<0:
  h=np.asarray(directions[0]['direction_unit_l2']);h=h/np.sqrt(np.mean(h*h))
  interval=[-np.inf,np.inf]
  for mi,di in zip(m,h):
   if abs(di)>1e-16:
    a,b=(lower-mi)/di,(upper-mi)/di;interval[0]=max(interval[0],min(a,b));interval[1]=min(interval[1],max(a,b))
  result['negative_curvature_ray']=dict(unit_RMS_direction=h,feasible_alpha_interval=interval,scope='Box feasible scalar interval; never coordinate clipping; 0.99 boundary safety on truncated requested alpha',original_requested_RMS=[.1,.5,1,2])
  # Preserve the decisive eigenvector before spending on ray forwards.
  for sign in [-1,1]:
   for rms in [.1,.5,1,2]:
    requested=sign*rms;alpha=max(.99*interval[0],min(.99*interval[1],requested))
    if abs(alpha)<1e-12:
     curves.append(dict(sign=sign,requested_RMS=rms,status='ZERO_FEASIBLE_INTERVAL',actual_alpha=alpha));continue
    trial=model.forward(m+alpha*h,label=f'{name}_NC_{alpha:.6g}')
    rr=float(np.linalg.norm(trial['y']-target)/np.sqrt(scale));rv=reg.evaluate(trial['m'])[0];loss=.5*rr**2+rv
    curves.append(dict(sign=sign,requested_RMS=rms,actual_alpha=float(alpha),actual_RMS=float(abs(alpha)),relative_residual=rr,objective=float(loss),objective_decrease=float(result['objective']-loss),quadratic_predicted_decrease=float(-alpha*(g@h)-.5*alpha*alpha*values[0]*np.dot(h,h)),material=trial['m'],status='PHYSICS_EVALUATED_FEASIBLE_NO_CLIP'))
 result['ray_probes']=curves
 result['costs']=dict(counts=dict(model.counts),operation_wall_seconds=dict(model.timings),operation_cpu_seconds=dict(model.cpu_timings),wall_seconds=time.perf_counter()-begun,cpu_seconds=time.process_time()-cpu,LU_reuse='Base LU retained for all HVP and GN Rayleigh; ray probes each pay one new LU',HVP_RHS_per_call=2*model.ns)
 return result

def report(results):
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 lines=['# Exact Hessian and box KKT mechanism experiment','','Status: offline public-Y-only diagnostic. No GT is opened by this program. Exact data Hessian includes the residual-Hessian term; the regularizer Hessian differentiates smoothed TV exactly. The eigenproblem is restricted to coordinates more than 1e-6 from both box bounds. No exact-rank, second-order global optimality, or formula accuracy claim is made.','','| Case | Free/lower/upper | Projected gradient infinity norm | Smallest algebraic eigenvalues | Eigensolver status | Best ray objective decrease |','|---|---|---:|---|---|---:|']
 fig,axes=plt.subplots(1,2,figsize=(12,4.5),constrained_layout=True)
 for i,r in enumerate(results['cases']):
  k=r['KKT'];ds=r['directions'];ev=r['eigensolver']['eigenvalues'];probes=r['ray_probes'];improve=max([q.get('objective_decrease',0) for q in probes]+[0])
  lines.append(f"| {r['case']} | {k['free_count']}/{k['active_lower']}/{k['active_upper']} | {k['projected_gradient_inf']:.4g} | {', '.join(f'{z:.5g}' for z in ev)} | {r['eigensolver']['status']} | {improve:.5g} |")
  for d in ds:axes[0].scatter(i-.08,d['exact_eigenvalue'],marker='v',color=f'C{i}');axes[0].scatter(i+.08,d['GN_plus_exact_reg_Rayleigh'],marker='o',color=f'C{i}')
  for q in probes:
   if 'objective' in q:axes[1].scatter(q['actual_alpha'],q['objective_decrease'],label=r['case'] if q is probes[0] else None,color=f'C{i}')
  lines += ['',f"## {r['case']}",'',f"Original objective={r['objective']:.9g}, relative residual={r['relative_residual']:.9g}. KKT lower/upper violations at gradient threshold1e-8: {k['lower_KKT_violations']}/{k['upper_KKT_violations']}. The reported KKT residual is a raw metric rather than a newly imposed acceptance gate.",'','| Eigenpair | Exact curvature | GN data + exact regularizer | Residual-Hessian contribution |','|---|---:|---:|---:|']
  for d in ds:lines.append(f"| {d['index']} | {d['exact_eigenvalue']:.7g} | {d['GN_plus_exact_reg_Rayleigh']:.7g} | {d['residual_Hessian_Rayleigh']:.7g} |")
  if len(ev) and ev[0]<0:lines += ['', 'A negative free-coordinate eigenvalue is a feasible local negative-curvature direction, up to the eigensolver tolerance. A strict saddle interpretation additionally requires sufficient first-order stationarity; the KKT metrics are retained so that this is not assumed. Finite rays evaluate the actual original regularized objective and use the scalar box-feasible interval, never coordinate clipping.','','| requested signed RMS | actual signed alpha | actual objective decrease | relative residual |','|---|---:|---:|---:|']
  for q in probes:
   if 'objective' in q:lines.append(f"| {q['sign']*q['requested_RMS']:+.2f} | {q['actual_alpha']:+.5g} | {q['objective_decrease']:+.7g} | {q['relative_residual']:.7g} |")
  lines += ['',f"Paid cost: {r['costs']['counts'].get('LU',0)} LU, {r['costs']['counts'].get('actual_complex_RHS',0)} complex RHS; {r['costs']['counts'].get('HVP',0)} HVP. Wall={r['costs']['wall_seconds']:.4f}s; CPU={r['costs']['cpu_seconds']:.4f}s. Eigensolver matvec cap=250, tol=1e-3. Failed/partial work remains paid."]
 axes[0].axhline(0,color='black',lw=.8);axes[0].set_xticks(range(len(results['cases'])),[r['case'] for r in results['cases']]);axes[0].set_title('Exact Hessian (v) versus GN + exact reg (o)');axes[0].set_ylabel('Rayleigh curvature, unit L2 directions');axes[0].grid(alpha=.2)
 axes[1].axhline(0,color='black',lw=.8);axes[1].set_title('Physical negative-curvature ray probes');axes[1].set_xlabel('signed alpha, direction unit RMS');axes[1].set_ylabel('original regularized objective decrease');axes[1].grid(alpha=.2)
 handles,labels=axes[1].get_legend_handles_labels()
 if handles:axes[1].legend()
 for ext in ['png','pdf','svg']:fig.savefig(ROOT/f'figures/hessian_geometry.{ext}',dpi=150,bbox_inches='tight')
 plt.close(fig);(ROOT/'diagnostics/hessian_geometry.md').write_text('\n'.join(lines)+'\n')

def main():
 begun,cpu=time.perf_counter(),time.process_time();results=dict(status='RUNNING',truth_access='NONE',formula='Exact normalized data Hessian plus exact generic regularizer Hessian; no formula checks or solver modifications',cases=[])
 cases=[('V1_holdout7_data_lbfgs',ROOT/'results/development_holdout7_v1/holdout7_data_lbfgs.json','obs_c3604cf3555d'),('old_holdout9_FullJ',OLD/'results/holdout/fullj_native/obs_42c48e795766.json','obs_42c48e795766'),('old_holdout8_FullJ',OLD/'results/holdout/fullj_native/obs_a63df7f20566.json','obs_a63df7f20566')]
 with threadpool_limits(limits=1):
  for name,source,oid in cases:
   results['cases'].append(run_case(name,source,oid));save(results)
 results.update(status='COMPLETED_PAID_CURVATURE_EXPERIMENT',wall_seconds=time.perf_counter()-begun,cpu_seconds=time.process_time()-cpu);save(results);report(results)
 print(json.dumps(dict(status=results['status'],wall_seconds=results['wall_seconds'],cpu_seconds=results['cpu_seconds'])),flush=True)
if __name__=='__main__':main()


def dense_freeface_followup():
 """One authorized representation change; preserve capped ARPACK cost and stop."""
 results=json.loads((ROOT/'diagnostics/hessian_geometry.json').read_text())
 with threadpool_limits(limits=1):
  for row in results['cases']:
   if row['case']=='old_holdout8_FullJ':continue
   nfree=row['KKT']['free_count']
   if nfree>80:
    row['dense_freeface_followup']=dict(status='NOT_RUN_FREE_DIMENSION_ABOVE_80',free_count=nfree,extra_LU=0,extra_complex_RHS=0)
    continue
   begun,cpu=time.perf_counter(),time.process_time();m=np.asarray(row['material']);source=ROOT.parent/row['source'];frozen=json.loads(source.read_text());cfg=frozen['config']
   op=dict(np.load(OLD/'public_data/operators.npz'));target=np.load(OLD/'public_data/holdout/obs_c3604cf3555d.npz')['Y']
   model=physics.Model(op);reg=multiscale.Regularizer(op['grid_shape'],cfg['smooth_weight'],cfg['tv_weight'],cfg['tv_epsilon'])
   state=model.forward(m,label='authorized_dense_freeface_followup');basefactor=model.factor;Q=state['E'];R=state['y']-target;scale=float(np.vdot(target,target).real)
   with model.operation('DENSE_FACE_BASE_ADJOINT',rhs=model.ns,adjoint_RHS=model.ns):adj=lu_solve(basefactor,model.SH@R,trans=2,check_finite=False)
   with model.operation('DENSE_FACE_REG_HESSIAN_SETUP',exact_reg_setup=1):
    d=reg.D@m;q=np.sqrt(d*d+reg.epsilon**2);w=(2*reg.smooth+reg.tv*reg.epsilon**2/q**3)/reg.edges
   indices=np.asarray(row['eigensolver']['free_indices']);H=np.empty((nfree,nfree))
   for c,index in enumerate(indices):
    h=np.zeros(model.n);h[index]=1
    with model.operation('DENSE_FACE_EXACT_HVP',rhs=2*model.ns,HVP=1,tangent_RHS=model.ns,incremental_adjoint_RHS=model.ns):
     U=lu_solve(basefactor,h[:,None]*Q,check_finite=False);dQ=model.G@U
     dadj=lu_solve(basefactor,model.GH@(h[:,None]*adj)+model.SH@(model.S@U),trans=2,check_finite=False)
     v=np.real(np.sum(dadj.conj()*Q+adj.conj()*dQ,axis=1))/scale
    with model.operation('DENSE_FACE_EXACT_REG_HVP',reg_HVP=1):v+=reg.D.T@(w*(reg.D@h))
    H[:,c]=v[indices]
   with model.operation('DENSE_FACE_EIGH',dense_eigh=1):values,vectors=np.linalg.eigh(.5*(H+H.T))
   fullvectors=np.zeros((model.n,min(3,nfree)));fullvectors[indices]=vectors[:,:min(3,nfree)]
   row['dense_freeface_followup']=dict(status='COMPLETED_ONE_REPRESENTATION_CHANGE',free_count=nfree,free_indices=indices,eigenvalues=values,first_three_direction_unit_L2=fullvectors,
    interpretation='Exact analytic fixed-free-face Hessian assembled by 31 coordinate HVP; symmetric representation averaged before dense eigh. No formula/precision comparison, no ARPACK restart, no tolerance sweep. Nonnegative fixed-face spectrum cannot certify full constrained local/global optimum.',
    Hessian_free_symmetric=.5*(H+H.T),truth_access='NONE',costs=dict(counts=dict(model.counts),operation_wall_seconds=dict(model.timings),operation_cpu_seconds=dict(model.cpu_timings),wall_seconds=time.perf_counter()-begun,cpu_seconds=time.process_time()-cpu,base_setup='Prior process had exited, so one new paid base forward+adjoint (1 LU/24 RHS) was required; retained followup LU serves all 31 HVP.'))
   print(json.dumps(clean(dict(case=row['case'],dense_min_eigenvalues=values[:3],counts=dict(model.counts),wall_seconds=time.perf_counter()-begun))),flush=True)
 results['followup_scope']='One authorized low-dimensional dense representation change, h7 only; h9 free200 NOT_RUN; h8 not repeated; original 250-matvec failures retained.'
 save(results)
 lines=['\n## Authorized dense fixed-face representation followup\n']
 for row in results['cases']:
  d=row.get('dense_freeface_followup')
  if not d:continue
  if d['status'].startswith('NOT_RUN'):lines.append(f"{row['case']}: NOT_RUN, free dimension {d['free_count']} exceeds80. No added physical work.\n");continue
  ev=d['eigenvalues'];cost=d['costs'];lines += [f"{row['case']}: explicit {d['free_count']}×{d['free_count']} fixed-face Hessian from coordinate HVP, followed by one symmetric dense eigh. Smallest three eigenvalues: "+', '.join(f'{x:.10g}' for x in ev[:3])+'.\n',f"Added paid work: {cost['counts']['HVP']} HVP, {cost['counts']['LU']} base LU, {cost['counts']['actual_complex_RHS']} complex RHS; wall {cost['wall_seconds']:.5f}s, CPU {cost['cpu_seconds']:.5f}s. The previous process had exited; this includes its necessary new base forward/adjoint, with LU reused by every HVP.\n",'The original 250-matvec ARPACK failure remains in the ledger. This is one representation change, not a tolerance sweep, derivative accuracy comparison, or formula test. A nonnegative fixed-free-face Hessian does not certify full box-constrained local or global optimality. No further ray or GT test is performed.\n']
 with (ROOT/'diagnostics/hessian_geometry.md').open('a') as f:f.write('\n'.join(lines))
