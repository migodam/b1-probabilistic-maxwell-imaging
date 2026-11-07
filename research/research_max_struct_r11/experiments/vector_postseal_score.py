"""Postseal true3D scores/plots; no new physics. Fresh GT needs root flag."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
 os.environ[name]='1'
import time
ENTRY=time.perf_counter()
from pathlib import Path
import sys,json,csv,argparse,shutil
import numpy as np
from scipy.stats import chi2
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
from research.research_max_struct_r10.experiments.score_vector_comparison import material_metrics

def read(p):return json.loads(Path(p).read_text())
def dump(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def table(p,rows):
 keys=list(dict.fromkeys(k for row in rows for k in row))
 with Path(p).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
def sq(x):return float(np.vdot(x,x).real)
def label(t):return t.get('method_label',t['arm'] if t['material_mode']=='ordinary' else t['material_mode'])
def update_scope(method):
 if method in ('sourceCIE','sourceDual','GlobalDual','globalDual'):
  return 'ordinary CIE current; ordinary scalar CIE material readout'
 return 'ordinary CIE current; '+method+' Green-denominator material readout only, not both blocks'
def stopping(r):
 s=r['physical_corrector']['stop']
 if 'ITERATIONS REACHED LIMIT' in s:return 'ACCEPTED_ITERATION_CAP_NOT_CONVERGENCE'
 if 'REQUEST_CAP' in s or 'EVALUATIONS' in s:return 'OBJECTIVE_REQUEST_CAP_NOT_CONVERGENCE'
 if 'ABNORMAL' in s:return 'LINE_SEARCH_STOP'
 if r['PG_inf']<=1e-8:return 'PROJECTED_GRADIENT_TARGET_NOT_GT_CERTIFICATE'
 return 'OTHER_STOP_PG_ABOVE_TARGET'

def figures(out,groups,materials,methods,rows,curves):
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 from matplotlib.colors import Normalize
 from matplotlib.cm import ScalarMappable
 out.mkdir(parents=True,exist_ok=False);norm=Normalize(0,8);cmap=plt.get_cmap('viridis')
 xyz=np.meshgrid(*[np.linspace(-.3,.3,13)]*3,indexing='ij')
 for group in groups:
  names=['GT']+methods;cols=len(names)
  fig=plt.figure(figsize=(3.5*cols,4.8),layout='constrained')
  for i,name in enumerate(names):
   ax=fig.add_subplot(1,cols,i+1,projection='3d');m=materials[group].get(name)
   if m is not None:
    v=m.reshape(12,12,12);colors=cmap(norm(v));colors[...,3]=.86
    ax.voxels(*xyz,v>.5,facecolors=colors,edgecolor=None)
   else:ax.text(0,0,0,'MISSING/FAILED',ha='center')
   ax.set(xlim=(-.3,.3),ylim=(-.3,.3),zlim=(-.3,.3),title=name,xlabel='x / lambda',ylabel='y / lambda',zlabel='z / lambda')
   ax.set_box_aspect((1,1,1));ax.view_init(24,38)
  fig.suptitle(group+' | true3D effective voxels chi>0.5 | fixed chi0..8')
  fig.colorbar(ScalarMappable(norm=norm,cmap=cmap),ax=fig.axes,shrink=.55,label='contrast chi')
  fig.savefig(out/(group+'_volume.png'),dpi=180);plt.close(fig)
  fig,axes=plt.subplots(3,cols,figsize=(3.1*cols,9),constrained_layout=True)
  for col,name in enumerate(names):
   m=materials[group].get(name)
   for axis in range(3):
    ax=axes[axis,col]
    if m is None:ax.text(.5,.5,'MISSING/FAILED',ha='center');ax.set_axis_off();continue
    ax.imshow(np.take(m.reshape(12,12,12),6,axis=axis).T,origin='lower',extent=(-.3,.3,-.3,.3),cmap=cmap,norm=norm)
    if axis==0:ax.set_title(name)
    if col==0:ax.set_ylabel(('x','y','z')[axis]+'=+0.025 lambda')
  fig.suptitle(group+' | fixed orthogonal cell-center slices | all1728 materials')
  fig.colorbar(ScalarMappable(norm=norm,cmap=cmap),ax=axes,shrink=.65,label='contrast chi')
  fig.savefig(out/(group+'_orthogonal.png'),dpi=180);plt.close(fig)
 for condition in ('matched12','cross18to12'):
  selected=[g for g in groups if g.endswith(condition)]
  if not selected:continue
  fig,axes=plt.subplots(1,2,figsize=(13,4.5),constrained_layout=True)
  x=np.arange(len(selected));width=.8/max(len(methods),1)
  for i,method in enumerate(methods):
   data={r['public_id']:r for r in rows if r['method']==method}
   for ax,key in zip(axes,['amplitude_normalized_RMSE','fixed_IoU']):
    ax.bar(x+(i-(len(methods)-1)/2)*width,[data[g].get(key,np.nan) for g in selected],width,label=method)
    ax.set_xticks(x,selected,rotation=30);ax.set_ylabel(key);ax.grid(axis='y',alpha=.2)
  axes[0].legend(fontsize=8);axes[1].set_ylim(0,1)
  fig.suptitle(condition+' | prescribed designs; matched/cross never pooled')
  fig.savefig(out/(condition+'_all_material_metrics.png'),dpi=180);plt.close(fig)
  fig,axes=plt.subplots(1,2,figsize=(13,4.5),constrained_layout=True)
  for i,method in enumerate(methods):
   data={r['public_id']:r for r in rows if r['method']==method}
   for ax,key in zip(axes,['Phi3','sigma_CN99_ratio']):
    ax.bar(x+(i-(len(methods)-1)/2)*width,[data[g].get(key,np.nan) for g in selected],width,label=method)
    ax.set_xticks(x,selected,rotation=30);ax.set_ylabel(key);ax.grid(axis='y',alpha=.2)
  axes[0].legend(fontsize=8);axes[1].axhline(1,color='red',linestyle='--',label='CN99 radius')
  fig.suptitle(condition+' | source-normalized Phi3 and sigma noise-ball ratio are separate metrics')
  fig.savefig(out/(condition+'_Phi3_and_sigma_CN99_separate.png'),dpi=180);plt.close(fig)
  for axis_name in ('physical_LU','CPU_seconds','wall_seconds'):
   fig,axes=plt.subplots(len(selected),1,figsize=(10,3.1*len(selected)),squeeze=False,constrained_layout=True)
   for ax,group in zip(axes[:,0],selected):
    for method in methods:
     points=[p for p in curves if p['public_id']==group and p['method']==method and p.get(axis_name) is not None]
     if points:ax.semilogy([p[axis_name] for p in points],[p['Phi3'] for p in points],label=method)
    host=next((r.get('host_dispatch_role') or r.get('host_platform') for r in rows if r['public_id']==group),'unknown host')
    ax.set(title=group+' | '+str(host),xlabel=axis_name,ylabel='source-normalized Phi3');ax.grid(alpha=.2)
   axes[0,0].legend(fontsize=8)
   fig.suptitle('Already-paid Phi3; old ordinary DEV intermediate CPU/wall unavailable')
   fig.savefig(out/(condition+'_Phi3_'+axis_name+'.png'),dpi=180);plt.close(fig)

def main():
 p=argparse.ArgumentParser();p.add_argument('--cohort',choices=['dev','confirmation'],required=True);p.add_argument('--authorization',required=True);p.add_argument('--authorization-note',required=True);p.add_argument('--truth-root',type=Path);args=p.parse_args()
 if args.cohort=='dev':
  if args.authorization!='OLD_EXPOSED_GT_ONLY_AFTER12SEALS':raise ValueError('Wrong old-GT authorization')
  tasks=read(ROOT/'configs/vector_dev_ordinary_queue.json')['tasks']+read(ROOT/'configs/vector_dev_denominator_queue.json')['tasks']
  truth_root=REPO.parent/'r10_worktree/research/research_max_struct_r10/sealed_truth/vector';public_root=ROOT/'delegated/vector/development_observations'
 else:
  if args.authorization!='ROOT_EXPLICIT_FRESH50_POSTSEAL':raise ValueError('Fresh GT needs later explicit root authorization')
  tasks=read(ROOT/'configs/vector_confirmation_queue.json')['tasks'];truth_root=ROOT/'delegated/vector/private_confirmation_truth';public_root=ROOT/'delegated/vector/confirmation_observations'
 if args.truth_root is not None:truth_root=args.truth_root
 required=12 if args.cohort=='dev' else 50
 if len(tasks)!=required:raise ValueError('Scorer must include entire frozen cohort')
 root=ROOT/'delegated/vector/inversions';seals=[];receipts={}
 for task in tasks:
  dest=root/task['task_id'];receipt=read(dest/'RECEIPT.json');seal=read(dest/'SEALED.json')
  if receipt['status']!='SEALED_PUBLIC_PREDICTION_NO_SCORE' or seal['status']!='PREDICTION_SEALED_BEFORE_ROOT_OFFLINE_SCORING':raise ValueError('Nonterminal/unsealed')
  if not (dest/'prediction.npz').exists():raise ValueError('Missing sealed prediction')
  if args.cohort=='confirmation' and not (dest/'COST_PREFIX_SEAL.json').exists():raise ValueError('Missing pre-GT prefix seal')
  receipts[task['task_id']]=receipt;seals.append(dict(task=task,seal=seal,receipt_bytes=(dest/'RECEIPT.json').stat().st_size,prediction_bytes=(dest/'prediction.npz').stat().st_size))
 out=ROOT/'delegated/vector'/(args.cohort+'_postseal_scores');out.mkdir(parents=True,exist_ok=False)
 seal_dir=out/'prediction_seals';seal_dir.mkdir()
 for task in tasks:shutil.copy2(root/task['task_id']/'SEALED.json',seal_dir/(task['task_id']+'.json'))
 dump(out/'ALL_PUBLIC_SEALS_BEFORE_GT.json',dict(status='FULL_COHORT_SEAL_BARRIER_BEFORE_FIRST_GT_READ',authorization=args.authorization_note,tasks=seals,GT_access_so_far=False,new_forward_RHS=0))
 # First private array read strictly follows the full-cohort public seal barrier.
 truth={};materials={};rows=[];states=[];curves=[];failure=[];model_errors=[];prefix_scores=[]
 affected_cases={t['public_id'].split('_')[0] for t in tasks if receipts[t['task_id']].get('infrastructure_attempt_metadata')}
 with np.load(ROOT/'delegated/vector/source_basis/BASIS.npz',allow_pickle=False) as z:U=z['right_vectors']
 with np.load(ROOT/'delegated/vector/PUBLIC_V04_FIXED_STATE.npz',allow_pickle=False) as z:ws=z['weights']
 for task in tasks:
  ident=task['task_id'];group=task['public_id'];case=group.split('_')[0];method=label(task);dest=root/ident;r=receipts[ident]
  if case not in truth:
   with np.load(truth_root/case/'grid12.npz',allow_pickle=False) as z:truth[case]={k:z[k].copy() for k in ('material','J','Y_clean')}
  gt=truth[case];materials.setdefault(group,{'GT':gt['material']})
  with np.load(public_root/group/'observations.npz',allow_pickle=False) as z:Y=z['Y']
  sigma=read(public_root/group/'PUBLIC.json')['sigma'];radius=float(chi2.ppf(.99,2*Y.size)/2)
  with np.load(dest/'prediction.npz',allow_pickle=False) as z:m=z['material'];pred=z['prediction'];J=z['J']
  residual=pred-Y;white=sq(residual)/sigma**2;phi=float(np.mean(np.sum(abs(residual)**2,axis=0)/np.sum(abs(Y)**2,axis=0)))
  error=J-gt['J'];total=float(np.sum(abs(error)**2*ws));captured=float(np.sum(abs(U.conj().T@error)**2*ws))
  linear_max=max(e.get('maximum_relative_linear_residual',0.) for e in r['paid_events'])
  row=dict(task_id=ident,case_id=case,public_id=group,condition='cross18to12' if 'cross' in group else 'matched12',method=method,
   **material_metrics(m,gt['material']),Phi3=r['physical_Phi'],Phi3_recomputed_saved_prediction=phi,
   sigma_white_squared_norm=white,sigma_CN99_radius_squared=radius,sigma_CN99_ratio=white/radius,sigma_CN99_qualified=white<=radius,
   current3N_source_weighted_relative_L2=float(np.sqrt(total/np.sum(abs(gt['J'])**2*ws))),
   current3N_unweighted_error_energy_h3=.05**3*sq(error),current3N_source_weighted_error_energy_h3=.05**3*total,
   current3N_source_basis_captured_energy=.05**3*captured,
   current3N_source_basis_uncaptured_energy=.05**3*max(0,total-captured),current_domain='complex5184x60 current, not scalar material projection',
   material_domain='1728 real shared isotropic chi',nested_CPU_seconds=r['CPU_seconds'],nested_wall_seconds=r['wall_seconds'],physical_LU=r['counts']['physical_LU'],
   primal_RHS=r['counts'].get('primal_RHS',r['counts'].get('physical_source_RHS')),adjoint_RHS=r['counts'].get('adjoint_RHS',r['counts'].get('physical_adjoint_RHS')),
   corrector_accepted=r['physical_corrector']['accepted'],corrector_requests=r['physical_corrector']['requests'],corrector_stop=r['physical_corrector']['stop'],PG_inf=r['PG_inf'],
   maximum_original_linear_residual=linear_max,physical_linear_qualified=bool(np.isfinite(linear_max) and linear_max<=2e-6),
   auxiliary_capped_stages=sum(s['stop']=='CAP' for s in r['stages']),
   optimizer_class=stopping(r),linalg_library=r.get('torch_linalg_library','Inherited default, ordinary DEV6'),update_scope=update_scope(method))
  row.update(infrastructure_retry=bool(r.get('infrastructure_attempt_metadata')),infrastructure_affected_object=case in affected_cases,
   nested_cost_scope='RETRY_ATTEMPT_ONLY_NOT_END_TO_END' if r.get('infrastructure_attempt_metadata') else 'NORMAL_FIRST_ATTEMPT_ARM',
   end_to_end_equal_cost_qualified=not bool(r.get('infrastructure_attempt_metadata')))
  row.update(host_platform=r.get('host_platform','win32_DEV'),host_dispatch_role=r.get('host_dispatch_role') or 'WINDOWS_CUDA_DEV')
  rows.append(row);materials[group][method]=m.copy()
  classes=[row['optimizer_class']]
  if r.get('infrastructure_attempt_metadata'):classes.append('INFRASTRUCTURE_QUERY_TIMEOUT_RETRY_PRIOR_PAID_COST_SEPARATE')
  if not row['material_legal']:classes.append('ILLEGAL_MATERIAL')
  if not row['physical_linear_qualified']:classes.append('PHYSICAL_LINEAR_SYSTEM_UNQUALIFIED')
  if not row['sigma_CN99_qualified']:classes.append('CROSSGRID_MISMATCH_OR_OPTIMIZATION_UNRESOLVED' if 'cross' in group else 'MATCHED_FINITE_MODEL_DATA_MISFIT_ABOVE_NOISE_BALL')
  failure.append(dict(task_id=ident,classes=';'.join(classes),continuous_material_metrics=True,no_new_scalar_success_gate=True))
  cumulative_LU=0
  for event in r['paid_events']:
   cumulative_LU+=event.get('counts_delta',{}).get('physical_LU',0)
   curves.append(dict(task_id=ident,public_id=group,method=method,state_id=event['state_id'],physical_LU=cumulative_LU,
    CPU_seconds=event.get('completion_arm_CPU_seconds'),wall_seconds=event.get('completion_arm_wall_seconds'),Phi3=event['Phi'],
    sigma_white_squared_norm=event.get('sigma_whitened_squared_norm'),legacy_intermediate_clock_unavailable='completion_arm_CPU_seconds' not in event,
    cost_scope=row['nested_cost_scope'],end_to_end_equal_cost_qualified=row['end_to_end_equal_cost_qualified']))
  for name in ['physical_initial','stage1_physical','stage2_physical','stage3_physical','physical_accepted80','physical_accepted160','physical_accepted240','physical_accepted320','prediction']:
   f=dest/(name+'.npz')
   if not f.exists():continue
   with np.load(f,allow_pickle=False) as z:sm=z['material']
   states.append(dict(task_id=ident,public_id=group,method=method,state=name,**material_metrics(sm,gt['material'])))
  if (dest/'COST_PREFIX_SEAL.json').exists():
   with np.load(dest/'COST_PREFIX_MATERIALS.npz',allow_pickle=False) as z:pm=z['material']
   for item in read(dest/'COST_PREFIX_SEAL.json')['prefixes']:
    score=dict(task_id=ident,public_id=group,method=method,**item,cost_scope=row['nested_cost_scope'],
     end_to_end_equal_cost_qualified=row['end_to_end_equal_cost_qualified'],infrastructure_affected_object=case in affected_cases)
    if item['material_row'] is not None:score.update(material_metrics(pm[item['material_row']],gt['material']))
    prefix_scores.append(score)
 for row in rows:
  control=next((x for x in rows if x['public_id']==row['public_id'] and x['method'] in ('GlobalDual','globalDual')),None)
  if control and row['method'] not in ('GlobalDual','globalDual'):
   row['NRMSE_delta_vs_matched_Global']=row['amplitude_normalized_RMSE']-control['amplitude_normalized_RMSE']
   if row['NRMSE_delta_vs_matched_Global']>0:
    entry=next(x for x in failure if x['task_id']==row['task_id']);entry['classes']+=';MATERIAL_NRMSE_WORSE_THAN_MATCHED_GLOBAL'
    if row['Phi3']<control['Phi3']:entry['classes']+=';LOWER_PHI_WITH_WORSE_MATERIAL_NRMSE'
 for group in materials:
  if 'cross' not in group:continue
  case=group.split('_')[0]
  with np.load(truth_root/case/'grid18.npz',allow_pickle=False) as z:fine=z['Y_clean']
  model_errors.append(dict(public_id=group,clean_grid12_grid18_relative=float(np.sqrt(sq(truth[case]['Y_clean']-fine)/sq(fine))),measurement_sigma_not_inflated=True,new_independent_materials=0))
 for name,data in [('MATERIAL_AND_PHYSICS',rows),('SAVED_STATE_MATERIAL_SCORES',states),('PAID_PHYSICAL_CURVES',curves),('COST_PREFIX_MATERIAL_SCORES',prefix_scores),('FAILURE_CLASSES',failure)]:table(out/(name+'.csv'),data)
 dump(out/'MATERIAL_AND_PHYSICS.json',rows);dump(out/'CROSSGRID_MODEL_ERROR.json',model_errors)
 if affected_cases:
  table(out/'MATCHED_SENSITIVITY_EXCLUDING_INFRASTRUCTURE_AFFECTED_OBJECT.csv',[z for z in rows if z['condition']=='matched12' and z['case_id'] not in affected_cases])
  table(out/'END_TO_END_QUALIFIED_FIXED_COST_PREFIX_SCORES.csv',[z for z in prefix_scores if z['end_to_end_equal_cost_qualified']])
  dump(out/'INFRASTRUCTURE_RETRY_SENSITIVITY_SCOPE.json',dict(affected_cases=sorted(affected_cases),all8matched_objects_retained_in_main_table=True,
   sensitivity_excludes_whole_affected_object=True,prior315step_costs_must_be_added_from_original_failed_process_ledger=True,
   failed_arm_retry_prefixes_are_attempt_only_not_end_to_end=True,no_GT_based_exclusion=True))
 dump(out/'FIXED_PREFIX_SELECTOR_STATUS.json',dict(
  status='NOT_REQUESTED_FOR_EXPOSED_DEV' if args.cohort=='dev' and not prefix_scores else 'REQUESTED_PUBLIC_FIXED_PREFIXES_SCORED',
  fixed_prefix_material_rows=len(prefix_scores),paid_trajectory_rows=len(curves),
  explanation='DEV did not request frozen fixedCPU/wall/LU prefix selectors; empty prefix table is not a result. Paid trajectories and80/160/240/320states are distinct.' if args.cohort=='dev' and not prefix_scores else
   'Only already-sealed public minPhi cost-prefix choices are scored; no oracleGT iterate selection'))
 methods=list(dict.fromkeys(label(t) for t in tasks));figures(out/'figures',list(materials),materials,methods,rows,curves)
 dump(out/'SCORING_COST.json',dict(status='COMPLETE',cohort=args.cohort,CPU_seconds=time.process_time(),wall_seconds=time.perf_counter()-ENTRY,
  new_Maxwell_RHS=0,new_inverse_calls=0,GT_access='Old exposed R10 DEV only' if args.cohort=='dev' else 'Fresh root explicitly authorized after50seals',
  independent_materials=2 if args.cohort=='dev' else 8,crossgrid_not_pooled=True,
  limitations=['No universal3D scalar success gate added','Budget/ftol stop and true PG qualification separated','Phi3 source-normalized, not sigma-ML','Source-basis remainder uncaptured, not weak/kernel','Old DEV default versus later selected library not pure algorithm speed comparison','All failed/worse cases retained'],
  scope='Whole local scorer imports/metrics/plots/exports, zero physics solve'))
 print(json.dumps(dict(status='COMPLETE',cohort=args.cohort,rows=len(rows),output=str(out))),flush=True)

if __name__=='__main__':main()
