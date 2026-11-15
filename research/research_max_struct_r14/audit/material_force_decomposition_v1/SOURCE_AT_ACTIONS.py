"""Finite no-GT diagnostic of the already-paid first material updates.

Reads exactly 7 existing intervention checkpoints x 8 branches. It never calls
step(), restores a physics model, evaluates a forward solve, or optimizes.
The metric derivative is the inherited complete gradient minus the frozen-P
constitutive gradient and inherited prior; no derivative is reimplemented.
"""
from __future__ import annotations
import os
for _k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[_k]='1'
import argparse,csv,inspect,json,math,sys,time,traceback
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path


def dump(path,value):
 with Path(path).open('x') as f:json.dump(value,f,indent=2,allow_nan=False);f.write('\n')
def table(out,name,rows):
 dump(out/(name+'.json'),rows)
 keys=list(dict.fromkeys(k for r in rows for k in r))
 with (out/(name+'.csv')).open('x',newline='') as f:
  writer=csv.DictWriter(f,fieldnames=keys);writer.writeheader();writer.writerows({k:json.dumps(v) if isinstance(v,(dict,list)) else v for k,v in r.items()} for r in rows)
def stamp(path):
 p=Path(path).resolve();s=p.stat();return dict(path=str(p),bytes=s.st_size,mtime_ns=s.st_mtime_ns)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--legacy-repo',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
 wall=time.perf_counter();start=time.process_time();R=a.root.resolve();W=a.legacy_repo.resolve();O=a.out.resolve()
 if O.exists():raise FileExistsError(O)
 O.mkdir(parents=True)
 rows=json.loads((R/'audit/dynamics_all192_v1/INTERVENTION_BRANCHES.json').read_text())
 expected_ids={'r13_fresh_0002_k00128','r13_fresh_0002_k00512','r13_fresh_0017_k00128','r13_fresh_0017_k00512','r13_fresh_0017_k02500','r13_fresh_0037_k00128','r13_fresh_0037_k00512'}
 expected_branches={'pair'+p+'_'+policy for p in ('00','01','10','11') for policy in ('retainPR','resetPR')}
 if len(rows)!=56 or {r['intervention'] for r in rows}!=expected_ids:raise ValueError('Only the original seven intervention sets are permitted')
 for ident in expected_ids:
  if {r['branch'] for r in rows if r['intervention']==ident}!=expected_branches:raise ValueError('Expected exactly four pairs x retain/reset PR')
 sourcepaths=[]
 for r in rows:
  folder=R/'audit/interventions'/r['intervention']/'branches'/r['branch']
  sourcepaths.extend([folder/'step_00000_state.npz',folder/'step_00001_state.npz',folder/'step_00001_physical.json'])
 dump(O/'PROTOCOL_BEFORE_ACTIONS.json',dict(created_utc=datetime.now(timezone.utc).isoformat(),diagnostic_only=True,truth_access='NONE',samples=56,interventions=sorted(expected_ids),branches=sorted(expected_branches),evaluation='m0 from step0; J1/Q1 and actual material beta/source Os from first completed step1',residual='e=J1-diag(m0)Q1',g_const='-rho Re sum_s conj(Q1_s) * Pm(m0)^*[(Pm(m0)e) Os]_s',g_metric='inherited full material gradient minus g_const minus inherited g_prior',work='dot(g_component,m1-m0); negative is alignment with negative component gradient',force_convention='Names g_const/g_metric/g_prior denote gradients. A descent force is their negative.',no_tests='No finite differences, amplitude search, new identity check, optimizer, forward, success label, or recovery predictor.',sources=[stamp(p) for p in sourcepaths],source_script=stamp(__file__)))
 sys.path[:0]=[str(W),str(R/'src')]
 import numpy as np
 from dualmetric_core import StatefulDualMetric,_load_state_npz,OneHopAbsorption
 from research.research_z_alg_r4_diversity.src.physics import Regularizer
 def norm(x):return float(np.linalg.norm(x))
 def angle(x,y):
  nx,ny=norm(x),norm(y)
  if nx==0 or ny==0:return None,'zero norm'
  return float(np.degrees(np.arccos(np.clip(float(np.dot(x,y))/(nx*ny),-1.,1.)))),None
 def ratio(x,y):return float(x/y) if y!=0 else None
 results=[];failures=[];events=[];total=Counter();arrays={};setup=Counter();reg_cache={}
 for index,record in enumerate(rows):
  before_cpu=time.process_time();before_wall=time.perf_counter();cost=Counter();folder=R/'audit/interventions'/record['intervention']/'branches'/record['branch'];label=record['intervention']+'__'+record['branch']
  try:
   s0=_load_state_npz(folder/'step_00000_state.npz');s1=_load_state_npz(folder/'step_00001_state.npz')
   if s1['iterations']!=s0['iterations']+1 or s1['pair']!=record['pair'] or s0['pair']!=record['pair']:raise ValueError('Saved first completed update identity mismatch')
   if not s1['history'] or s1['history'][-1]['global_iteration']!=s1['iterations']:raise ValueError('Missing actual first material update row')
   native=s1['history'][-1];m0=np.asarray(s0['m'],float);m1=np.asarray(s1['m'],float);J1=s1['J'];Q1=s1['Q'];dm=m1-m0;shape=tuple(s1['shape'])
   if shape not in reg_cache:
    reg_cache[shape]=Regularizer(shape);setup['original_regularizer_sparse_difference_setup']+=1
   regularizer=reg_cache[shape];prior,gprior,_=regularizer.evaluate(m0);gprior=np.asarray(gprior,float);cost['validation_regularizer_gradient_evaluations']+=1
   def reg(z):
    if not np.array_equal(z,m0):raise ValueError('Diagnostic gradient is evaluated only at saved m0')
    return prior,gprior
   # Read-only facade: only inherited _P/_PH/material_value_gradient are called.
   # No constructor, Fourier basis, Green action, current update or physics model.
   core=StatefulDualMetric.__new__(StatefulDualMetric);core.n=s1['n'];core.pair=s1['pair'];core.beta0=s1['beta0'];core.beta1=s1['beta1'];core.rho=float(s1['frozen_config'].get('rho',1.));core.Os=s1['Os'];core.cost=cost;core.reg=reg
   core.blocks=OneHopAbsorption.__new__(OneHopAbsorption);core.blocks.__dict__.update(s1['onehop_cache']);core.blocks.cost=cost
   e=J1-m0[:,None]*Q1;family=core.pair[1]
   value,gfull=core.material_value_gradient(m0,J1,Q1)
   Pe=core._P(m0,e,family);Wres=Pe@core.Os;cost['source_metric_products']+=1
   PhW=core._PH(m0,Wres,family)
   cost['validation_P_apply_calls']+=1;cost['validation_P_adjoint_apply_calls']+=1
   gconst=-core.rho*np.real(np.sum(Q1.conj()*PhW,axis=1))
   gmetric=np.asarray(gfull,float)-gconst-gprior
   gs=dict(const=gconst,metric=gmetric,prior=gprior,full=np.asarray(gfull,float))
   beta=core.beta0 if family=='0' else core.beta1
   actual=float(native['beta_material'])
   if beta!=actual:raise ValueError('Material beta differs from actual native first-step metadata')
   physical=json.loads((folder/'step_00001_physical.json').read_text())
   r=dict(intervention=record['intervention'],case_id=record['case_id'],branch=record['branch'],pair=record['pair'],reset_PR=record['reset_PR'],source_global_iteration=record['source_global_iteration'],completed_global_iteration=s1['iterations'],stage_before=s0['stage'],stage_material=s1['stage'],stage_iteration_material=s1['iteration'],material_family='P'+family,material_beta=beta,beta0=core.beta0,beta1=core.beta1,rho=core.rho,cells=len(m0),sources=J1.shape[1],source_metric='actual saved step1.Os; full right-source metric',source_Os_trace=float(np.trace(core.Os).real),source_Os_offdiagonal_Frobenius=norm(core.Os-np.diag(np.diag(core.Os))),material_value_at_start_with_J1_Q1=float(value),prior_value_at_start=float(prior),constitutive_defect_Frobenius=norm(e),transformed_defect_Frobenius=norm(Pe),material_step_L2=norm(dm),material_step_Linf=float(np.max(np.abs(dm))),material_inner_nit=native.get('material_inner_nit'),material_inner_success=native.get('material_inner_success'),material_fallback=native.get('material_fallback'),direction_restarted=native.get('direction_restarted'),PR_plus=native.get('PR_plus'),line_alpha=native.get('line_alpha'),existing_step_angle_to_negative_physical_Phi_gradient_degrees=physical.get('material_step_angle_to_negative_physical_Phi_gradient_degrees'),existing_step_angle_scope='Reuse already-paid step1 physical audit; no new gradient/forward')
   for name,g in gs.items():
    r['g_'+name+'_L2']=norm(g);r['g_'+name+'_Linf']=float(np.max(np.abs(g)));r['g_'+name+'_work_on_actual_dm']=float(np.dot(g,dm));ang,reason=angle(dm,-g);r['step_angle_to_negative_g_'+name+'_degrees']=ang;r['step_angle_to_negative_g_'+name+'_null_reason']=reason
   r['metric_to_const_L2_ratio']=ratio(norm(gmetric),norm(gconst));r['metric_to_full_L2_ratio']=ratio(norm(gmetric),norm(gfull));r['prior_to_full_L2_ratio']=ratio(norm(gprior),norm(gfull))
   for x,y in [('const','metric'),('const','prior'),('metric','prior'),('full','const'),('full','metric'),('full','prior')]:
    ang,reason=angle(gs[x],gs[y]);r['angle_g_'+x+'_g_'+y+'_degrees']=ang;r['angle_g_'+x+'_g_'+y+'_null_reason']=reason
   for name,data in dict(m0=m0,m1=m1,J1=J1,Q1=Q1,e=e,Pe=Pe,Os=core.Os,delta_m=dm,g_const=gconst,g_metric=gmetric,g_prior=gprior,g_full=gfull).items():arrays[label+'__'+name]=data
   r['new_validation_operation_counts']=dict(cost);results.append(r)
  except Exception as exc:failures.append(dict(intervention=record['intervention'],branch=record['branch'],error=repr(exc),traceback=traceback.format_exc()))
  finally:
   total.update(cost);events.append(dict(intervention=record['intervention'],branch=record['branch'],CPU_seconds=time.process_time()-before_cpu,wall_seconds=time.perf_counter()-before_wall,new_validation_operation_counts=dict(cost),physical_LU=0,physical_complex_RHS=0,optimizer_calls=0))
 table(O,'FIRST_STEP_MATERIAL_FORCES',results);table(O,'VALIDATION_COST_EVENTS',events);table(O,'FAILURES',failures);np.savez_compressed(O/'FIRST_STEP_FORCE_ARRAYS.npz',**arrays)
 # Mechanical summaries only, with no truth, success labels or predictor fitting.
 grouped=[]
 for p in ('00','01','10','11'):
  for reset in (False,True):
   group=[r for r in results if r['pair']==p and r['reset_PR']==reset]
   item=dict(pair=p,reset_PR=reset,samples=len(group),unit='seven saved checkpoint states; exposed posthoc dependent branches, not independent recovery samples')
   for k in ['g_const_L2','g_metric_L2','g_prior_L2','g_full_L2','metric_to_const_L2_ratio','metric_to_full_L2_ratio','angle_g_const_g_metric_degrees','step_angle_to_negative_g_const_degrees','step_angle_to_negative_g_metric_degrees','step_angle_to_negative_g_full_degrees','g_const_work_on_actual_dm','g_metric_work_on_actual_dm','g_prior_work_on_actual_dm','g_full_work_on_actual_dm']:
    vals=[r[k] for r in group if r[k] is not None]
    item[k+'_median']=float(np.median(vals)) if vals else None;item[k+'_min']=min(vals) if vals else None;item[k+'_max']=max(vals) if vals else None
   grouped.append(item)
 table(O,'MECHANICAL_GROUP_SUMMARY',grouped)
 import matplotlib;matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'savefig.bbox':'tight'})
 fig,axes=plt.subplots(2,2,figsize=(13,8));colors={'00':'#0072B2','01':'#D55E00','10':'#009E73','11':'#CC79A7'}
 for p in ('00','01','10','11'):
  for reset in (False,True):
   group=[r for r in results if r['pair']==p and r['reset_PR']==reset];marker='x' if reset else 'o';label=p+(' reset PR' if reset else ' retain PR')
   axes[0,0].scatter([r['g_const_L2'] for r in group],[r['g_metric_L2'] for r in group],c=colors[p],marker=marker,label=label,s=32,alpha=.8)
   axes[0,1].scatter([r['step_angle_to_negative_g_const_degrees'] for r in group],[r['step_angle_to_negative_g_full_degrees'] for r in group],c=colors[p],marker=marker,s=32,alpha=.8)
   axes[1,0].scatter([r['g_const_work_on_actual_dm'] for r in group],[r['g_metric_work_on_actual_dm'] for r in group],c=colors[p],marker=marker,s=32,alpha=.8)
   axes[1,1].scatter([r['angle_g_const_g_metric_degrees'] for r in group],[r['metric_to_const_L2_ratio'] for r in group],c=colors[p],marker=marker,s=32,alpha=.8)
 axes[0,0].set(xlabel='‖g_const‖₂',ylabel='‖g_metric‖₂',xscale='log',yscale='log');axes[0,1].set(xlabel='Actual step angle to −g_const (degrees)',ylabel='Actual step angle to −g_full (degrees)',xlim=(0,180),ylim=(0,180));axes[0,1].plot([0,180],[0,180],c='gray',ls=':',lw=.8)
 axes[1,0].set(xlabel='g_const · actual Δm',ylabel='g_metric · actual Δm');axes[1,0].set_xscale('symlog',linthresh=1e-10);axes[1,0].set_yscale('symlog',linthresh=1e-10);axes[1,0].axhline(0,c='gray',lw=.6);axes[1,0].axvline(0,c='gray',lw=.6)
 axes[1,1].set(xlabel='Angle(g_const, g_metric) (degrees)',ylabel='‖g_metric‖₂ / ‖g_const‖₂',xlim=(0,180),yscale='log')
 for ax in axes.ravel():ax.grid(alpha=.2)
 handles,labels=axes[0,0].get_legend_handles_labels();fig.legend(handles,labels,ncol=4,loc='upper center',bbox_to_anchor=(.5,.94));fig.suptitle('56 saved first material steps: frozen-P constitutive term and metric derivative remainder');fig.tight_layout(rect=(0,.02,1,.85));fig.text(.5,.005,'Source: existing seven intervention full-state checkpoints × 8 branches; actual m0 and J1/Q1/Os. No truth, new forward, or optimizer.',ha='center',fontsize=8)
 for ext in ('png','pdf'):fig.savefig(O/('MATERIAL_FORCE_COMPONENTS.'+ext),dpi=175)
 plt.close(fig)
 receipt=dict(status='COMPLETE' if len(results)==56 and not failures else 'INCOMPLETE_RETAINED_FAILURES',diagnostic_only=True,truth_access='NONE',expected_branches=56,completed_branches=len(results),failures=len(failures),checkpoint_states=7,gradient_location='m0 from branch step0, fixed updated auxiliary J1/Q1 from branch step1',prior=dict(source=inspect.getfile(Regularizer),class_name='Regularizer',shape=list(next(iter(reg_cache))) if reg_cache else None,smooth=1e-4,tv=1e-4,epsilon=.01,qualification='Exact original MaxwellModel.reg constructor defaults, inherited intervention build_config regularizer=model.reg.evaluate(m)[:2]. Sparse setup only, no MaxwellModel constructed.'),gradient_implementation=dict(source=inspect.getfile(StatefulDualMetric),full_gradient='StatefulDualMetric.material_value_gradient; P0 scalar full derivative and P1 inherited onehop_material_value_gradient',frozen_P='StatefulDualMetric._P/_PH at m0 using actual saved material family/beta and onehop cache',metric_derivative='g_full-g_const-g_prior; intentionally no redundant identity check'),setup_operation_counts=dict(setup),new_validation_operation_counts=dict(total),historical_counts_recharged=False,new_physical_LU=0,new_complex_RHS=0,new_Green_or_FFT_calls=0,optimizer_calls=0,CPU_seconds=time.process_time()-start,whole_process_CPU_seconds_at_receipt=time.process_time(),wall_seconds=time.perf_counter()-wall,CPU_scope='Imports + checkpoint I/O/copies + cached block actions + inherited prior/full gradients + tables/arrays/plots; receipt serialization follows timestamp. Launcher separately records whole OS child CPU.',source_script=stamp(__file__),finished_utc=datetime.now(timezone.utc).isoformat(),interpretation_scope='Finite seven-state mechanism diagnostic only. Not a new chain rule, Jacobian/stability theorem, recovery predictor, or success fit.')
 dump(O/'RECEIPT.json',receipt)
 report=['# 首个材料步的非根梯度分解','',f"固定7个已有checkpoint × 四metric pair × PR保留/重置，共{len(results)}/56枝完成；失败{len(failures)}。只读step0的m0和首个完成step1的J1/Q1、实际beta/Os、缓存onehop块。没有GT读取、forward、Green/FFT传播、current更新或optimizer。",'', '材料梯度在实际readout开始位置m0、current更新后的固定J1/Q1处计算。g_const冻结P，只对e=J1−mQ1求导；g_metric定义为原完整gradient−g_const−g_prior；P0/P1都计算。先验直接复用原Regularizer的默认构造，与原intervention的MaxwellModel.reg一致。','', 'g均表示梯度分量，下降力是−g。CSV列work_on_actual_dm=dot(g,m1−m0)，负号表示与该负梯度分量同向。记录三项和完整梯度的L2/Linf范数、实际材料步与−g的夹角、各梯度分量间夹角，以及已有paid physical-step gradient angle；零向量角保持null。','', '样本是7个已曝光完整状态的56个依赖分枝，不是56个独立恢复任务；不使用truth/recovery标签，不拟合成功，不扫尺度、不重写完整导数、不重复机器精度恒等式验证。新增block动作/源权乘积/先验计算全部作为validation收费；历史checkpoint计数不重收。','', '第一步分解用于描述ideal根线性化之外的非根metric项。是否支持进一步数学解释由主线程判断，不能由本图推出全局恢复、渐近稳定或新链式法则。','', f"本次正文CPU {receipt['CPU_seconds']:.6f}s；receipt前整个进程CPU {receipt['whole_process_CPU_seconds_at_receipt']:.6f}s；新LU/RHS/optimizer均0。精确完成后OS子进程CPU另见LAUNCH_RECEIPT.json。"]
 (O/'REPORT_ZH.md').write_text('\n'.join(report)+'\n');print(json.dumps(receipt,indent=2))
 return 0 if receipt['status']=='COMPLETE' else 1
if __name__=='__main__':raise SystemExit(main())
