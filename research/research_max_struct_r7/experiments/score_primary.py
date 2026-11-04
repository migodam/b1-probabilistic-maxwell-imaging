"""Independent offline R7 scorer: no solver imports, no new physical calls.

`seal` reads public predictions/metadata only. `score` requires the parent's
explicit ALL_PRIMARY_SEALED message and matching literal CLI authorization.
Never run score before that message. Truth/mapping opens occur only in score.
"""
from pathlib import Path
from datetime import datetime,timezone
import argparse,csv,json,time
import numpy as np
from scipy.stats import chi2
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
ARMS=('cie_single','cie_dual','global_dual','a_material','a_plus','physical_lbfgs')
TERMINAL=('SEALED_PUBLIC_PREDICTION','FAILED','WORKER_FAILED','CANCELLED','HARD_TIMEOUT')

def load(p):return json.loads(Path(p).read_text())
def resolve(p):return Path(p) if Path(p).is_absolute() else REPO/p
def stamp(p):
 p=Path(p);s=p.stat();return dict(path=str(p),size=s.st_size,mtime_ns=s.st_mtime_ns)
def unchanged(s):
 p=Path(s['path']);return p.is_file() and stamp(p)==s
def write(p,obj):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
 with p.open('x') as f:json.dump(obj,f,indent=2,allow_nan=False)
def csvwrite(p,rows):
 if not rows:return
 with Path(p).open('x',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in rows for k in r)));w.writeheader();w.writerows(rows)
def authorize(value):
 if value!='ALL_PRIMARY_SEALED':raise RuntimeError('Parent must explicitly send ALL_PRIMARY_SEALED before scoring; no GT access authorized')
def otsu(m):
 m=np.asarray(m).ravel()
 if np.ptp(m)<=1e-12:return np.zeros(m.size,bool),float(m[0])
 h,e=np.histogram(m,bins=256);c=(e[:-1]+e[1:])/2;w=np.cumsum(h);s=np.cumsum(h*c);den=w*(w[-1]-w)
 v=np.full(256,-np.inf);good=den>0;v[good]=(s[-1]*w[good]-s[good]*w[-1])**2/den[good]
 threshold=float(c[np.argmax(v)]);return m>threshold,threshold

def material_metrics(m,truth,support,amplitude):
 m=np.asarray(m).ravel();truth=np.asarray(truth).ravel();support=np.asarray(support).ravel().astype(bool)
 if m.shape!=(576,) or truth.shape!=m.shape or support.shape!=m.shape:raise ValueError('Full576-variable material array required')
 if not np.isfinite(m).all():raise ValueError('Nonfinite material')
 rmse=float(np.linalg.norm(m-truth)/np.sqrt(m.size));nrmse=rmse/float(amplitude);seg,t=otsu(m);union=np.count_nonzero(seg|support)
 iou=float(np.count_nonzero(seg&support)/union) if union else 1.
 quantized=m.dtype==np.float32;ltol=float(np.spacing(np.float32(.02))) if quantized else 1e-12;utol=float(np.spacing(np.float32(20.))) if quantized else 1e-12
 legal=bool(np.min(m)>=.02-ltol and np.max(m)<=20+utol)
 return dict(material_nrmse=nrmse,amplitude_normalized_nrmse=nrmse,RMSE=rmse,MAE=float(np.mean(abs(m-truth))),relative_l2_material_error=float(np.linalg.norm(m-truth)/max(np.linalg.norm(truth),1e-30)),otsu_iou=iou,otsu_threshold=t,material_legal=legal,material_nrmse_pass=bool(nrmse<=.15),iou_pass=bool(iou>=.7),material_image_success=bool(nrmse<=.15 and iou>=.7 and legal),material_min=float(np.min(m)),material_max=float(np.max(m)),legal_box_quantization_tolerance=dict(lower=ltol,upper=utol,reason='one float32ULP for saved budget materials only' if quantized else 'float64 numerical tolerance'))

def public_observation(task):
 result={}
 for p in task['public_files']:
  with np.load(resolve(p),allow_pickle=False) as z:
   for k in ('Y','sigma','noise_radius'):
    if k in z.files:result[k]=z[k].copy()
 if 'Y' not in result:raise ValueError('Missing public Y')
 return result

def noise_metrics(yp,public,inverse_crime):
 Y=public['Y'];yp=np.asarray(yp)
 if yp.shape!=Y.shape or not np.isfinite(yp).all():raise ValueError('Prediction/observation shape or finiteness mismatch')
 residual=float(np.linalg.norm(yp-Y));relative=residual/max(float(np.linalg.norm(Y)),1e-300);sigma=float(public.get('sigma',0));radius=float(np.sqrt(chi2.ppf(.99,2*Y.size)/2))
 return dict(residual_norm=residual,relative_residual=relative,complex_observations=int(Y.size),sigma=sigma,CN99_radius=radius,whitened_residual=residual/sigma if sigma>0 else None,noise_CN99_pass=bool(residual/sigma<=radius) if sigma>0 else None,noise_gate_scope='same-grid independent circular complex Gaussian noise discrepancy' if inverse_crime else 'RAW NOISE-ONLY CRITERION: cross-grid physical model mismatch is additional and unbudgeted',model_mismatch_budget='NONE' if inverse_crime else 'NOT_SPECIFIED: sigma is not a full error budget')

def fixed_prior(m):
 # Exact frozen MaxwellModel->R4 Regularizer defaults, inspected source; no solver import.
 image=np.asarray(m,float).reshape(24,24);d=np.r_[np.diff(image,axis=0).ravel(),np.diff(image,axis=1).ravel()]/19.5
 return float(1e-4*np.mean(d*d)+1e-4*np.mean(np.sqrt(d*d+.01**2)-.01))

def feasible_witness(rows,cases,mapping,truths):
 records=[]
 for ident,case in cases.items():
  if not case['inverse_crime']:
   records.append(dict(case_id=ident,status='NOT_AVAILABLE_DIFFERENT_OPERATOR',reason='Stored fine-gridY_clean is not a physical prediction of coarse-gridm_ref; no new forward authorized for this scorer',new_physical_calls=0));continue
  meta=mapping[ident];truth=truths[ident][0]
  with np.load(resolve(meta['truth_path']),allow_pickle=False) as z:clean=z['Y_clean'].copy();noise=z['noise'].copy()
  public=public_observation(case);Y=public['Y'];rn=float(np.linalg.norm(clean-Y));rr=rn/float(np.linalg.norm(Y));reg=fixed_prior(truth);phi=.5*rr**2+reg
  entries=[]
  for row in [r for r in rows if r['case_id']==ident]:
   value=row.get('original_Phi');gap=None if value is None else value-phi;better=bool(gap is not None and gap>1e-12+1e-8*abs(phi));stop=row.get('corrector_stop') or '';converged=row.get('corrector_success') is True;cap='LIMIT' in stop.upper() or 'MAXIMUM' in stop.upper()
   label='KNOWN_BETTER_FEASIBLE_PHI_WITH_BUDGET_STOP' if better and cap else 'KNOWN_BETTER_FEASIBLE_PHI_WITH_REPORTED_NUMERICAL_CONVERGENCE' if better and converged else 'KNOWN_BETTER_FEASIBLE_PHI_OTHER_STOP' if better else 'NO_STRICTLY_BETTER_GT_WITNESS_PHI_IDENTIFIED'
   entries.append(dict(arm=row['arm'],status=row['status'],material_success=row['material_image_success'],prediction_Phi=value,prediction_minus_feasible_GT_Phi=gap,known_better_feasible_point=better,prediction_relative_residual=row.get('relative_residual'),prediction_to_GT_residual_ratio=row['relative_residual']/rr if row.get('relative_residual') is not None and rr>0 else None,corrector_stop=stop,reported_optimizer_success=converged,projected_gradient_inf=row.get('projected_gradient_inf'),projected_gradient_availability='PAID_FINAL_INCUMBENT_GRADIENT_SAVED' if 'projected_gradient_inf' in row else 'NOT_SAVED_AT_THIS_INCUMBENT',diagnostic=label))
  records.append(dict(case_id=ident,status='SAME_GRID_OFFLINE_FEASIBLE_WITNESS',GT_Phi=phi,GT_prior=reg,GT_residual_norm=rn,GT_relative_residual=rr,GT_whitened_residual=rn/float(public['sigma']),GT_legal_box=bool(truth.min()>=.02 and truth.max()<=20),saved_noise_consistency_error=float(np.linalg.norm((Y-clean)-noise)),predictions=entries))
 return dict(status='POST_SEAL_OFFLINE_ONLY',same_grid_witnesses=24,cross_grid_not_available=6,prior_source='MaxwellModel.reg=R4Regularizer((24,24)): smooth=tv=1e-4,epsilon=.01,D normalized by19.5',interpretation='A lower known feasible objective diagnoses failure to attain that objective; it does not prove material uniqueness, a global optimum, or a full-space basin barrier. Reported numerical convergence and projected gradient are separate from recovery. No witness is used to select predictions or alter a solver.',new_physical_calls=0,new_physical_LU=0,new_physical_RHS=0,records=records)

def seal(manifest,out):
 manifest=resolve(manifest);tasks=load(manifest)['tasks']
 if len(tasks)!=180 or len({t['case_id'] for t in tasks})!=30 or {t['arm'] for t in tasks}!=set(ARMS):raise RuntimeError('Require frozen30×6 primary task set')
 pairs=[(t['case_id'],t['arm']) for t in tasks]
 if len(set(pairs))!=180:raise RuntimeError('Duplicate case/arm task in primary manifest')
 if any(sum(t['case_id']==ident and t['arm']==arm for t in tasks)!=1 for ident in {t['case_id'] for t in tasks} for arm in ARMS):raise RuntimeError('Unbalanced primary manifest')
 entries=[];evidence=[stamp(manifest),stamp(ROOT/'data_public/cases.json')]
 cfg=load(manifest).get('config')
 if cfg:evidence.append(stamp(resolve(cfg)))
 for config in sorted({t['config'] for t in tasks if t.get('config')}):evidence.append(stamp(resolve(config)))
 for task in tasks:
  directory=resolve(task['output_dir']);rec=load(directory/'RECEIPT.json');status=rec['status']
  if status not in TERMINAL:raise RuntimeError('Nonterminal primary task '+task['task_id'])
  if status=='SEALED_PUBLIC_PREDICTION':
   for name in ['prediction.npz','prediction.json','BUDGET_STATES.npz']:
    if not(directory/name).is_file():raise RuntimeError('Missing final prediction/budget artifact '+task['task_id']+'/'+name)
  entries.append(dict(task=task,status=status,files=[stamp(p) for p in sorted(directory.iterdir()) if p.is_file()],public_files=[stamp(resolve(p)) for p in task['public_files']]))
 write(out,dict(status='ALL_PRIMARY_PUBLIC_OUTPUTS_SEALED',created_utc=datetime.now(timezone.utc).isoformat(),manifest=str(manifest),evidence=evidence,tasks=entries,truth_access='NONE',hash_checks='NONE',objects=30,tasks_count=180))
 print(json.dumps(dict(status='PUBLIC_SEAL_COMPLETE',objects=30,tasks=180,output=str(out))))

def verify_seal(path):
 seal=load(path)
 if seal['status']!='ALL_PRIMARY_PUBLIC_OUTPUTS_SEALED':raise RuntimeError('Not a complete primary public seal')
 for s in seal['evidence']+[s for t in seal['tasks'] for s in t['files']+t['public_files']]:
  if not unchanged(s):raise RuntimeError('Public sealed artifact changed '+s['path'])
 return seal

def mapping_after_authorization():
 # This function must only be called after authorize()+verify_seal().
 data=load(ROOT/'sealed_truth/mapping_PRIVATE.json')
 return {r['case_id']:r for r in data['records']}

def truth_after_authorization(meta):
 with np.load(resolve(meta['truth_path']),allow_pickle=False) as z:return z['m'].ravel().copy(),z['mask'].ravel().astype(bool)

def score(seal_path,out,authorization):
 authorize(authorization);started,cpu=time.perf_counter(),time.process_time();sealed=verify_seal(seal_path)
 out=Path(out)
 if out.exists():raise FileExistsError('Write-once score output already exists '+str(out))
 cases={r['case_id']:r for r in load(ROOT/'data_public/cases.json')['records']};mapping=mapping_after_authorization()
 rows=[];pictures={};truths={};labels={};checkpoint_rows=[];budget_available=[]
 for entry in sealed['tasks']:
  task=entry['task'];ident=task['case_id'];case=cases[ident];meta=mapping[ident];directory=resolve(task['output_dir']);rec=load(directory/'RECEIPT.json');cost=rec.get('physical_model_cost') or {};corr=rec.get('corrector') or {}
  if ident not in truths:truths[ident]=truth_after_authorization(meta)
  truth,support=truths[ident];family=meta.get('parameters',{}).get('family',meta.get('family','UNKNOWN'));labels[ident]=family
  row=dict(task_id=task['task_id'],case_id=ident,object_id=case['object_id'],arm=task['arm'],status=entry['status'],condition=case['condition'],family=family,cohort='R7_NEW30_CONFIRMATION',generation_grid=case['generation_grid'],inversion_grid=case['inversion_grid'],inverse_crime=case['inverse_crime'],amplitude=case['amplitude'],sources=case['sources'],aperture=case['aperture'],noise_relative_complex_RMS=case['noise_relative_complex_RMS'],denominator_included=True,material_image_success=False,same_grid_material_noise_success=False if case['inverse_crime'] else None,material_nrmse=None,otsu_iou=None,relative_residual=None,original_Phi=rec.get('physical_objective'),cpu_seconds=rec.get('cpu_seconds'),wall_seconds=rec.get('wall_seconds'),LU=cost.get('LU'),source_RHS=cost.get('source_RHS'),adjoint_RHS=cost.get('adjoint_RHS'),actual_complex_RHS=cost.get('actual_complex_RHS'),corrector_stop=corr.get('stop'),corrector_success=corr.get('success'),corrector_iterations=corr.get('iterations'),corrector_evaluations=corr.get('evaluations'),auxiliary_cpu_seconds=rec.get('auxiliary_cpu_seconds'),auxiliary_iterations=rec.get('auxiliary_iterations'),auxiliary_execution_receipt=rec.get('auxiliary_execution_receipt'))
  if entry['status']=='SEALED_PUBLIC_PREDICTION':
   public=public_observation(task)
   with np.load(directory/'prediction.npz',allow_pickle=False) as z:pred={k:z[k].copy() for k in z.files}
   m=pred['m'].ravel();row['prediction_material_dtype']=str(m.dtype);row['prediction_material_cells']=int(m.size);metrics=material_metrics(m,truth,support,case['amplitude']);row.update(metrics);row.update(noise_metrics(pred['Y_prediction'],public,case['inverse_crime']))
   row['same_grid_material_noise_success']=bool(metrics['material_image_success'] and row['noise_CN99_pass']) if case['inverse_crime'] else None
   row['cross_grid_material_and_raw_noise_pass']=bool(metrics['material_image_success'] and row['noise_CN99_pass']) if not case['inverse_crime'] else None
   row['residual_receipt_difference']=row['relative_residual']-float(rec['physical_relative_residual'])
   row['independent_Phi_from_saved_prediction']=.5*row['relative_residual']**2+fixed_prior(m);row['Phi_receipt_difference']=row['independent_Phi_from_saved_prediction']-float(rec['physical_objective'])
   if 'gradient' in pred:
    grad=pred['gradient'].ravel();row['projected_gradient_inf']=float(np.linalg.norm(m-np.clip(m-grad,.02,20),np.inf))
   pictures[(ident,task['arm'])]=m
   # Grade paid checkpoint incumbents only; material float32 precision is stated.
   with np.load(directory/'BUDGET_STATES.npz',allow_pickle=False) as z:columns=[str(v) for v in z['columns']];costs=z['costs'].copy();materials=z['m'].copy()
   if columns!=['CPU_seconds_including_process_prefix','LU','complex_RHS','original_Phi'] or costs.ndim!=2 or costs.shape[1]!=4 or materials.shape!=(len(costs),576):raise ValueError('Budget snapshot schema mismatch')
   if np.any(np.diff(costs[:,:3],axis=0)<-1e-9):raise ValueError('Nonmonotone cumulative costs')
   final_quantized=material_metrics(materials[-1],truth,support,case['amplitude'])
   budget_available.append(dict(case_id=ident,arm=task['arm'],checkpoints=len(costs),last_budget_cpu=float(costs[-1,0]),final_process_cpu=rec['cpu_seconds'],final_CPU_unobserved_tail=float(rec['cpu_seconds']-costs[-1,0]),final_material_float32_max_difference=float(np.max(abs(materials[-1].astype(float)-m))),final_nrmse_quantization_difference=final_quantized['material_nrmse']-metrics['material_nrmse']))
   for index,(c,mat) in enumerate(zip(costs,materials)):
    met=material_metrics(mat,truth,support,case['amplitude']);checkpoint_rows.append(dict(case_id=ident,arm=task['arm'],index=index,CPU=float(c[0]),LU=float(c[1]),complex_RHS=float(c[2]),original_Phi=float(c[3]),**met))
  rows.append(row)
 out.mkdir(parents=True);write(out/'SCORES.json',dict(rows=rows,objects=30,arms=list(ARMS),material_success_definition='amplitude-normalized RMSE<=.15, prediction-only256-bin Otsu IoU>=.7, legal [.02,20]',same_grid_joint_definition='material success AND independent CN99 noise discrepancy;24objects',cross_grid_definition='six36-to24 model-mismatch cases: raw residual/noise-only pass retained; sigma is not a full error budget',physical_scoring_LU=0,physical_scoring_RHS=0));csvwrite(out/'SCORES.csv',rows)
 write(out/'CHECKPOINT_SCORES.json',dict(rows=checkpoint_rows,scope='Already-paid original-Phi incumbent, every10 observed paid evaluations plus auxiliary/final points, m float32; no interpolation or GT selection',objects=30,new_physical_calls=0));write(out/'BUDGET_AUDIT.json',budget_available)
 write(out/'GT_FEASIBLE_WITNESS.json',feasible_witness(rows,cases,mapping,truths))
 plot_panels(ROOT/'FIGURES/primary',rows,truths,pictures,labels)
 write(out/'SCORE_RECEIPT.json',dict(status='COMPLETE',authorization='ALL_PRIMARY_SEALED explicitly received',seal=str(seal_path),created_utc=datetime.now(timezone.utc).isoformat(),wall_seconds=time.perf_counter()-started,cpu_seconds=time.process_time()-cpu,rows=len(rows),objects=30,truth_access='POST_COMPLETE180_TASK_PUBLIC_SEAL_ONLY',new_physical_calls=0,new_physical_LU=0,new_physical_RHS=0,selection='Fixed final paid public-Phi incumbent; no GT root/iterate selection',scorer_family_visibility='Generator source schema revealed deterministic family ordering during code inspection before seal; private mapping/GT unopened until this scoring pass; not used by any optimizer'))
 print(json.dumps(dict(status='SCORED',rows=len(rows),objects=30,out=str(out))))

def plot_panels(folder,rows,truths,pictures,labels):
 import matplotlib;matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 folder=Path(folder);folder.mkdir(parents=True,exist_ok=True);ids=sorted(truths)
 # All30 objects and all6 arms, six groups of5; no selection by score.
 for group in range(6):
  selected=ids[5*group:5*group+5];fig,axes=plt.subplots(5,7,figsize=(18,13),squeeze=False,layout='constrained')
  for i,ident in enumerate(selected):
   amplitude=next(r['amplitude'] for r in rows if r['case_id']==ident);vmax=max([.1+amplitude]+[float(np.max(pictures[(ident,arm)])) for arm in ARMS if (ident,arm) in pictures])
   axes[i,0].imshow(truths[ident][0].reshape(24,24),vmin=.02,vmax=vmax,cmap='viridis');axes[i,0].set_title(ident+'\n'+labels[ident]+'\nGT | row max %.2f'%vmax,fontsize=9)
   for j,arm in enumerate(ARMS,1):
    row=next(r for r in rows if r['case_id']==ident and r['arm']==arm);m=pictures.get((ident,arm))
    if m is not None:axes[i,j].imshow(m.reshape(24,24),vmin=.02,vmax=vmax,cmap='viridis');title=f"{arm}\nN={row['material_nrmse']:.3f} IoU={row['otsu_iou']:.2f}\nr={row['relative_residual']:.3g}"
    else:title=arm+'\n'+row['status']
    axes[i,j].set_title(title,fontsize=9)
   for ax in axes[i]:ax.set_xticks([]);ax.set_yticks([])
  fig.suptitle('All30 objects ×6 arms, group%d/6; shared row color range [.02,max(truth,all six predictions)]; raw576-variable material maps'%(group+1),fontsize=12)
  fig.savefig(folder/('all30_group_%02d.png'%(group+1)),dpi=150);plt.close(fig)

def main():
 p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='command',required=True)
 q=s.add_parser('seal');q.add_argument('--manifest',default=str(ROOT/'configs/primary_final_manifest.json'));q.add_argument('--output',required=True)
 q=s.add_parser('score');q.add_argument('--seal',required=True);q.add_argument('--output',required=True);q.add_argument('--authorization',required=True)
 a=p.parse_args()
 if a.command=='seal':seal(a.manifest,a.output)
 else:score(a.seal,a.output,a.authorization)
if __name__=='__main__':main()
