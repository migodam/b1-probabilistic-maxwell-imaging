"""Prospective E3_REV1 minimal weak-ray repair. Original E0-E3 immutable."""
from __future__ import annotations
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'):os.environ[k]='1'
import json,time,sys,argparse
from pathlib import Path
import numpy as np
from scipy.linalg import null_space
import exploration as old
P=old.P;OUT=P/'results/exploration_rev1';AUD=P/'audit/exploration_rev1';CONF=P/'config/exploration_rev1_freeze.json';CAP=4500
ARMS=['B_RAY','C_RAY']
# Runtime output override affects this process only; frozen original file is unchanged.
old.OUT=OUT

def freeze():
 AUD.mkdir(parents=True,exist_ok=True);OUT.mkdir(parents=True,exist_ok=True)
 d=dict(id='R9_S_E3_REV1',frozen_utc=time.time(),reason='Exposed-DEV mechanism repair after E3 fixed-coordinate decoder projection caused nonpositive predicted decrease; no truth read',parents=[f'R7_LOCKED_{i:03d}' for i in range(8)],radii=[.01,.02],arms=ARMS,weak_rank=4,initial_direction='same original R7 anchor prior-aware material-normalized v0 for B/C',B='Fixed anchor Gaussian QR weak frame; GN correction in weak orthogonal complement plus original signed-progress tangent; original R7 radial retraction scales all weak coordinates preserving weak ray orientation',C='Recomputed current Gaussian QR weak frame each required J; same initial direction and signed-progress/radial retraction rules; no fixed anchor weak-ray claim',extra_radial_tangent_constraint=False,weak_amplitude='free to change under radial retraction',max_states=3,trust_initial=.5,material_trust_fraction=.5,max_latent_radius=4.,prior_and_gates='original R7 unchanged',algorithm_J='only before required next correction; endpoint J after operational clock separate paid audit',all8_complete_required=True,additional_RHS_cap=CAP,original_E3_immutable=True,fresh_confirmation=False,noGT=True,extra_path_physics=False,root_aggregate_before=13976,aggregate_bound=18476,root_reserve_atleast=2500)
 if CONF.exists():
  d['frozen_utc']=json.loads(CONF.read_text())['frozen_utc']
  if json.loads(CONF.read_text())!=d:raise RuntimeError('revision freeze mismatch')
 else:old.write(CONF,d)
 src=AUD/'FROZEN_REV1_SOURCE.py'
 if src.exists() and src.read_bytes()!=Path(__file__).read_bytes():raise RuntimeError('revision source mismatch')
 if not src.exists():src.write_bytes(Path(__file__).read_bytes())
 print('REV1_FROZEN zeroRHS',flush=True)

class Book(old.Book):
 def __init__(self,job,max_wall=600):
  if not CONF.exists() or (AUD/'FROZEN_REV1_SOURCE.py').read_bytes()!=Path(__file__).read_bytes():raise RuntimeError('REV1 not frozen')
  super().__init__('REV1_'+job,max_wall);self.cap=CAP
  if self.total['RHS']>CAP:raise RuntimeError('revision cap')

def substep(prior,base,p,v,J,residual,arm,g0,trust,distance):
 st=time.perf_counter()
 if arm=='C_RAY':
  g=old.frame(prior,p.x,J);M=g['M'];Qa=g['Qweak'];U,_,Vh=np.linalg.svd(Qa.T@g0['Qweak']);Qw=Qa@U@Vh;dynamic=g['svd_seconds']
 else:M=old.real_material(old.material_derivative(prior,p.x));Qw=g0['Qweak'];dynamic=0.
 u=old.real_material(old.material_derivative(prior,base.x)@v);u/=np.linalg.norm(u)
 _,sv,Vh=np.linalg.svd(M,full_matrices=False);keep=sv>max(sv[0]*1e-8,1e-10);E=Vh[keep].T
 Z=E@null_space(np.vstack((Qw.T,u@M))@E,rcond=1e-10)
 A=np.vstack((J@Z,Z));r=np.r_[residual,p.x];step=Z@(-np.linalg.lstsq(A,r,rcond=1e-10)[0]);factor=min(1.,trust/max(np.linalg.norm(step),1e-300),.5*distance/max(np.linalg.norm(M@step),1e-300));step*=factor
 return step,dict(correction_rank=Z.shape[1],weak_linear_defect=float(np.linalg.norm(Qw.T@step)),dynamic_geometry_arithmetic_seconds=dynamic,dynamic_SVD_new_RHS=0,J_already_paid=True,substep_wall_seconds=time.perf_counter()-st),Qw

def radial(prior,base,current,proposal,distance,Qweak):
 return old.match_distance(prior,base.x,lambda h:base.x+h*(proposal-base.x),distance)
old.substep=substep;old.radius_projection=radial

def ray_metrics(g0,x,base,v):
 w0=g0['Qweak'].T@v;w=g0['Qweak'].T@(x-base.x);coef=float(w@w0/(w0@w0));normal=w-coef*w0
 return dict(anchor_weak_ray_absolute_defect=float(np.linalg.norm(normal)),anchor_weak_ray_relative_defect=float(np.linalg.norm(normal)/max(np.linalg.norm(w),1e-30)),anchor_weak_ray_signed_amplitude=coef,weak_coordinate_norm=float(np.linalg.norm(w)))

def run(parent):
 dest=OUT/'parents'/parent
 if (dest/'COMPLETE.json').exists():raise RuntimeError('no revision replay')
 dest.mkdir(parents=True,exist_ok=True);inp=old.AUD/parent/'ONLINE_INPUT.npz';b=old.point(inp);y=old.arrays(inp)['observation'];idx=int(parent.rsplit('_',1)[1]);rows=[]
 with Book(parent) as book:
  target=old.Target(book)
  try:
   g0=old.frame(target.prior,b.x,b.J);v=g0['v'];unit=old.real_material(old.material_derivative(target.prior,b.x)@v);unit/=np.linalg.norm(unit)
   for radius in [.01,.02]:
    order=ARMS if idx%2==0 else ARMS[::-1]
    for arm in order:
     dd=dest/arm/f'RMS_{radius:.4f}';dd.mkdir(parents=True,exist_ok=True);before=book.total.copy();st=time.perf_counter();xp,match=old.match_distance(target.prior,b.x,lambda h:b.x+h*v,radius)
     if xp is None:row=dict(parent=parent,arm=arm,radius=radius,status=match['status']);pp=None
     else:
      pp=old.correct(target,b,y,xp,v,radius,dd,arm,g0)
      if pp is None:row=dict(parent=parent,arm=arm,radius=radius,status='NO_FEASIBLE_ENDPOINT')
      else:
       row=old.metrics(target,pp,b,y,unit,radius);row.update(parent=parent,arm=arm,radius=radius,status='EXECUTED',initial_match=match,**ray_metrics(g0,pp.x,b,v))
       row['qualified']=all(row[k] for k in ['distance_matched','progress_pass','likelihood_competitive','noise_ball','reasonable_prior_heuristic','legal'])
       trace=old.arrays(dd/'ACTUAL_CORRECTION_PATH.npz')['x_nodes'];old.write(dd/'WEAK_RAY_TRACE.json',dict(arm=arm,fixed_anchor_ray_required=arm=='B_RAY',nodes=[dict(node=i,**ray_metrics(g0,x,b,v)) for i,x in enumerate(trace)]))
     opend=book.total.copy();row['operational_cost']={k:opend[k]-before[k] for k in ['RHS','states','LU','forward_RHS','derivative_RHS']};row['operational_wall_seconds']=time.perf_counter()-st;book._event('ALGORITHM_SELECTION_CLOSED',arm=arm,radius=radius,operational_cost=row['operational_cost'],operational_wall_seconds=row['operational_wall_seconds']);ast=time.perf_counter()
     if pp is not None:
      if pp.J is None:
       if pp is target._live_point:target.jacobian(pp)
       else:
        restored=target.evaluate(pp.x,arm+'_selected_J_AUDIT_restore');target.jacobian(restored);pp.J=restored.J.copy();pp.Jnative=restored.Jnative.copy();pp.T=restored.T.copy()
      target.persist_cache(pp,dd/'ENDPOINT.npz')
     row['diagnostic_cost']={k:book.total[k]-opend[k] for k in ['RHS','states','LU','forward_RHS','derivative_RHS']};row['diagnostic_wall_seconds']=time.perf_counter()-ast;row['cost']={k:book.total[k]-before[k] for k in ['RHS','states','LU','forward_RHS','derivative_RHS']};row['wall_seconds']=time.perf_counter()-st;old.write(dd/'RESULT.json',row);rows.append(row);print('REV1_ENDPOINT',parent,arm,radius,row.get('Phi'),row['status'],book.total['RHS'],flush=True)
   old.write(dest/'COMPLETE.json',dict(parent=parent,role='EXPOSED_DEV_REV1',rows=rows,ledger=book.total,online_GT_used=False))
  finally:target.close()

def summarize():
 rows=[]
 for p in sorted((OUT/'parents').glob('*/COMPLETE.json')):rows+=json.loads(p.read_text())['rows']
 table=[];paired=[];stops=[]
 for radius in [.01,.02]:
  for arm in ARMS:
   a=[r for r in rows if r['arm']==arm and r['radius']==radius];ok=[r for r in a if r['status']=='EXECUTED'];chs=[]
   for r in ok:
    ch=json.loads((OUT/'parents'/r['parent']/arm/f'RMS_{radius:.4f}'/'CORRECTION.json').read_text())['history'];chs.append(ch)
    for h in ch:
     if 'stop' in h:stops.append(dict(parent=r['parent'],arm=arm,radius=radius,stop=h['stop']))
   table.append(dict(arm=arm,radius=radius,n=len(a),available=len(ok),qualified=sum(r.get('qualified',False) for r in ok),likelihood_competitive=sum(r['likelihood_competitive'] for r in ok),pair_le1=sum(r['whitened_pair_L2']<=1 for r in ok),median_pair=None if not ok else float(np.median([r['whitened_pair_L2'] for r in ok])),median_DeltaPhi=None if not ok else float(np.median([r['DeltaPhi'] for r in ok])),operational_RHS=sum(r['operational_cost']['RHS'] for r in a),audit_RHS=sum(r['diagnostic_cost']['RHS'] for r in a),total_RHS=sum(r['cost']['RHS'] for r in a),operational_wall=sum(r['operational_wall_seconds'] for r in a),one_state_only=sum(len(ch)==1 for ch in chs),ray_defect_max=max((r['anchor_weak_ray_relative_defect'] for r in ok),default=0.),dynamic_geometry_seconds=sum(h.get('dynamic_geometry_arithmetic_seconds',0.) for ch in chs for h in ch)))
  for i in range(8):
   p=f'R7_LOCKED_{i:03d}';a=[r for r in rows if r['parent']==p and r['radius']==radius]
   if len(a)!=2 or any(r['status']!='EXECUTED' for r in a):continue
   q={r['arm']:r for r in a};b,c=q['B_RAY'],q['C_RAY'];ba=old.arrays(OUT/'parents'/p/'B_RAY'/f'RMS_{radius:.4f}'/'ENDPOINT.npz');ca=old.arrays(OUT/'parents'/p/'C_RAY'/f'RMS_{radius:.4f}'/'ENDPOINT.npz');paired.append(dict(parent=p,radius=radius,B_qualified=b['qualified'],C_qualified=c['qualified'],latent_endpoint_difference=float(np.linalg.norm(ba['x']-ca['x'])),material_endpoint_difference_RMS=float(np.linalg.norm(ba['material']-ca['material'])/24),B_minus_C_Phi=b['Phi']-c['Phi'],B_pair=b['whitened_pair_L2'],C_pair=c['whitened_pair_L2'],B_operational_RHS=b['operational_cost']['RHS'],C_operational_RHS=c['operational_cost']['RHS'],B_operational_wall=b['operational_wall_seconds'],C_operational_wall=c['operational_wall_seconds']))
 led=[json.loads(x) for x in (OUT/'RESOURCE_LEDGER.jsonl').read_text().splitlines()];jobs=[r for r in led if r['event']=='JOB_END'];resource=dict(total=led[-1]['after'],controlled_wall=sum(r['wall_seconds'] for r in jobs),CPU_seconds=sum(r['CPU_seconds'] for r in jobs),actual_RHS=sum(r['counters'].get('actual_forward_RHS',0)+r['counters'].get('actual_derivative_RHS',0) for r in jobs),array_IO={k:sum(r['array_IO'].get(k,0) for r in jobs) for k in ['array_save_attempts','array_save_completed','array_save_bytes','array_save_seconds']})
 old.write(OUT/'SUMMARY.json',dict(status='COMPLETE' if len(rows)==32 else 'PARTIAL',rows=rows,table=table,paired=paired,stops=stops,resource=resource,original_negative_preserved=True,development_revision=True,noGT=True,no_posterior_claim=True))
 lines=['# E3_REV1 weak-ray repair (exposed DEV only)','', 'This separate preregistered development revision changes the failed fixed-coordinate radius projection to original R7 radial retraction, preserving the fixed weak ray in B while allowing its amplitude to change. Original E3 remains immutable.','', '| arm | r | available | qualified | pair<=1 | median pair | op RHS | audit RHS | op wall(s) | one-state-only | max ray defect |','|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
 for r in table:lines.append('|'+ '|'.join(str(r[k]) for k in ['arm','radius','available','qualified','pair_le1','median_pair','operational_RHS','audit_RHS','operational_wall','one_state_only','ray_defect_max'])+'|')
 lines+=['',f'Additional actual physics: {resource["total"]["RHS"]} RHS, {resource["total"]["states"]} states/LUs, {resource["total"]["failures"]} solver failures. Aggregate exploration RHS={13976+resource["total"]["RHS"]}.', '','All operational and audit costs include failed/rejected attempts. Current selected endpoint J diagnostics occur after operational selection timing. A and Frozen-J remain original-run references; their wall times are from a separate timing window. No fresh confirmation, posterior, or continuous-fiber claim.']
 (OUT/'REPORT.md').write_text('\n'.join(lines)+'\n');print('REV1_SUMMARY',len(rows),resource['total']['RHS'],flush=True)

def main():
 p=argparse.ArgumentParser();p.add_argument('command',choices=['freeze','all','summarize']);a=p.parse_args()
 if a.command=='freeze':freeze()
 elif a.command=='all':
  for i in range(8):run(f'R7_LOCKED_{i:03d}');summarize()
 else:summarize()
if __name__=='__main__':main()
