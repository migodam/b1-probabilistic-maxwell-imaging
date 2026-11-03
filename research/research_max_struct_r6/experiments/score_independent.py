"""Independent offline scorer. No solver imports or new physical calls.

seal: public outputs only; score: explicit parent authorization, complete unchanged seal.
"""
from pathlib import Path
from datetime import datetime,timezone
import argparse,csv,json,sys,time
import numpy as np
from scipy.stats import chi2
ROOT=Path(__file__).resolve().parents[1]; REPO=ROOT.parents[1]

def load(p): return json.loads(Path(p).read_text())
def resolve(p): return Path(p) if Path(p).is_absolute() else REPO/p
def write(p,v):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
 with p.open('x') as f: json.dump(v,f,indent=2,allow_nan=False)
def stamp(p):
 s=p.stat();return dict(path=str(p),size=s.st_size,mtime_ns=s.st_mtime_ns)
def unchanged(s):
 p=Path(s['path']);return p.exists() and stamp(p)==s

def seal(manifest,out):
 data=load(manifest); tasks=[];evidence=[stamp(manifest),stamp(ROOT/'data_public/cases.json')]
 if data.get('config'):evidence.append(stamp(resolve(data['config'])))
 for task in data['tasks']:
  directory=resolve(task['output_dir']); receipt=directory/'RECEIPT.json'
  if not receipt.exists():raise RuntimeError('No terminal receipt: '+task['task_id'])
  status=load(receipt)['status']
  if status not in ('SEALED_PUBLIC_PREDICTION','FAILED','WORKER_FAILED','CANCELLED','HARD_TIMEOUT','UNKNOWN'):
   raise RuntimeError('Nonterminal task: '+task['task_id'])
  if status=='SEALED_PUBLIC_PREDICTION' and not (directory/'prediction.npz').exists():raise RuntimeError('Missing sealed prediction')
  files=[stamp(p) for p in sorted(directory.iterdir()) if p.is_file()]
  public=[stamp(resolve(p)) for p in task['public_files']]
  tasks.append(dict(task=task,status=status,files=files,public_files=public))
 write(out,dict(status='COMPLETE_PUBLIC_SEAL',created_utc=datetime.now(timezone.utc).isoformat(),manifest=str(manifest),evidence=evidence,tasks=tasks,truth_access='NONE',hash_checks='NONE'))
 print(json.dumps(dict(status='SEALED',tasks=len(tasks),seal=str(out))))

def otsu(m):
 if np.ptp(m)<=1e-12:return np.zeros(m.size,bool),float(m[0])
 hist,edges=np.histogram(m,bins=256);centers=(edges[:-1]+edges[1:])/2
 w=np.cumsum(hist);m1=np.cumsum(hist*centers);den=w*(w[-1]-w);between=np.full(256,-np.inf);good=den>0
 between[good]=(m1[-1]*w[good]-m1[good]*w[-1])**2/den[good]
 threshold=float(centers[np.argmax(between)]);return m>threshold,threshold

def prior(m,p):
 d=np.r_[np.diff(m.reshape(24,24),axis=0).ravel(),np.diff(m.reshape(24,24),axis=1).ravel()]/float(p['prior_width'])
 return float(float(p['smooth_weight'])*np.mean(d*d)+float(p['tv_weight'])*np.mean(np.sqrt(d*d+float(p['tv_epsilon'])**2)-float(p['tv_epsilon'])))

def score(seal_path,out,batch,authorized):
 if not authorized:raise RuntimeError('Explicit parent instruction and --authorize-after-seal required before GT access')
 started,cpu=time.perf_counter(),time.process_time(); sealed=load(seal_path)
 assert sealed['status']=='COMPLETE_PUBLIC_SEAL'
 for entry in sealed['evidence']+[s for task in sealed['tasks'] for s in task['files']+task['public_files']]:
  if not unchanged(entry):raise RuntimeError('Sealed evidence changed: '+entry['path'])
 if Path(out).exists():raise FileExistsError('Write-once scoring output already exists')
 cases={r['case_id']:r for r in load(ROOT/'data_public/cases.json')['records']}
 historical=load(ROOT/'sealed_truth/scoring/historical_mapping_PRIVATE.json')
 fresh=load(ROOT/'sealed_truth/confirmation/mapping_PRIVATE.json')['records']
 mapping={r['case_id']:r for r in historical+fresh}
 legacy={}
 with (REPO/'research/research_fbg_r3_experiments/results/e3_all_scores.csv').open() as f:
  for row in csv.DictReader(f):
   if row.get('physical_data_gate'):legacy[row['case_id']]=float(row['physical_data_gate'])
 rows=[]; pictures={}; truths={};private_labels={}
 for entry in sealed['tasks']:
  task=entry['task'];ident=task['case_id'];case=cases[ident];meta=mapping[ident]
  isnew=case['cohort']=='NEW_CONFIRMATION_UNSCORED'
  label=meta.get('shape',meta.get('family',ident))
  stratum='NEW6_CONFIRMATION' if isnew else ('DEV9_DIGITS' if label.startswith('digit') else 'DEV6_DIFFICULT')
  private_labels[ident]=label
  directory=resolve(task['output_dir']);rec=load(directory/'RECEIPT.json')
  cost=rec.get('physical_model_cost') or {};corr=rec.get('corrector') or {}
  row=dict(batch=batch,task_id=task['task_id'],case_id=ident,arm=task['arm'],cohort=case['cohort'],stratum=stratum,label=label,status=entry['status'],denominator_included=True,wall_seconds=rec.get('wall_seconds'),cpu_seconds=rec.get('cpu_seconds'),LU=cost.get('LU'),source_RHS=cost.get('source_RHS'),adjoint_RHS=cost.get('adjoint_RHS'),actual_complex_RHS=cost.get('actual_complex_RHS'),stop=corr.get('stop'),corrector_success=corr.get('success'),corrector_iterations=corr.get('iterations'),corrector_evaluations=corr.get('evaluations'),auxiliary_status=rec.get('auxiliary_status'),auxiliary_iterations=rec.get('auxiliary_iterations'),original_Phi=rec.get('physical_objective'),relative_residual=rec.get('physical_relative_residual'),legacy_joint=None,confirmation_joint=None,discrepancy_gate=None,projected_gradient_inf=None,gradient_availability='NOT_SAVED_AT_FINAL_INCUMBENT',material_nrmse=None,RMSE=None,otsu_iou=None)
  # Truth is opened only after the complete public seal has been verified.
  with np.load(resolve(meta['truth_path']),allow_pickle=False) as z:truth=z['m'].ravel().copy();support=z['mask'].ravel().astype(bool)
  truths[ident]=truth
  if entry['status']=='SEALED_PUBLIC_PREDICTION':
   with np.load(directory/'prediction.npz',allow_pickle=False) as z: pred={k:z[k] for k in z.files}
   m=pred['m'].ravel();yp=pred['Y_prediction'];assert m.shape==(576,) and yp.shape==(48,12)
   if not np.isfinite(m).all() or not np.isfinite(yp).all():raise ValueError('Nonfinite prediction')
   public={**case}
   for name in task['public_files']:
    with np.load(resolve(name),allow_pickle=False) as z:
     for k in ('Y','sigma','noise_radius'): 
      if k in z.files:public[k]=z[k]
   residual=float(np.linalg.norm(yp-public['Y']));relative=residual/float(np.linalg.norm(public['Y']))
   rmse=float(np.linalg.norm(m-truth)/np.sqrt(576));nrmse=rmse/float(case['amplitude']);mask,threshold=otsu(m)
   union=np.count_nonzero(mask|support);overlap=float(np.count_nonzero(mask&support)/union) if union else 1.
   legal=bool(m.min()>=.02-1e-12 and m.max()<=20+1e-12);image=bool(nrmse<=.15 and overlap>=.7)
   phi=.5*relative**2+prior(m,public)
   row.update(RMSE=rmse,material_nrmse=nrmse,amplitude_normalized_nrmse=nrmse,otsu_iou=overlap,otsu_threshold=threshold,image_gate=image,material_gate=bool(nrmse<=.15),material_legal=legal,relative_residual=relative,residual_norm=residual,original_Phi=phi,Phi_receipt_difference=phi-float(rec['physical_objective']),residual_receipt_difference=relative-float(rec['physical_relative_residual']))
   sigma=float(public.get('sigma',0.));radius=float(np.sqrt(chi2.ppf(.99,2*yp.size)/2))
   row.update(public_sigma=sigma if sigma>0 else None,complex_observations=int(yp.size),CN99_radius=radius,whitened_residual=residual/sigma if sigma>0 else None,discrepancy_gate=bool(residual/sigma<=radius) if sigma>0 else None,discrepancy_status='AVAILABLE_PUBLIC_POSITIVE_SIGMA' if sigma>0 else 'NOT_APPLICABLE_CLEAN_OR_SIGMA_UNAVAILABLE')
   if isnew:row['confirmation_joint']=bool(image and legal and row['discrepancy_gate'])
   elif ident in legacy:
    row.update(legacy_data_gate_type='FBG-R3 relative residual threshold',legacy_data_threshold=legacy[ident],legacy_data_pass=bool(relative<=legacy[ident]),legacy_joint=bool(image and legal and relative<=legacy[ident]))
   else:
    oldradius=float(public['noise_radius']);passed=bool(sigma>0 and residual/sigma<=oldradius)
    row.update(legacy_data_gate_type='original R2/R4 whitened residual radius',legacy_data_threshold=oldradius,legacy_data_pass=passed,legacy_joint=bool(image and legal and passed))
   if 'gradient' in pred:
    grad=pred['gradient'].ravel();pg=m-np.clip(m-grad,.02,20.)
    row.update(projected_gradient_inf=float(np.linalg.norm(pg,np.inf)),gradient_inf=float(np.linalg.norm(grad,np.inf)),gradient_availability='SAVED_PAID_FINAL_INCUMBENT')
   pictures[(ident,task['arm'])]=m
  row['auxiliary_execution_receipt']=rec.get('auxiliary_execution_receipt')
  row['corrector_raw_summary']=corr
  rows.append(row)
 # Failed tasks remain in each denominator; missing metric values are not imputed.
 summaries=[]
 for stratum in ('DEV9_DIGITS','DEV6_DIFFICULT','NEW6_CONFIRMATION'):
  for arm in sorted({r['arm'] for r in rows}):
   group=[r for r in rows if r['stratum']==stratum and r['arm']==arm]
   if not group:continue
   summary=dict(batch=batch,stratum=stratum,arm=arm,objects=len(group),scored=sum(r['material_nrmse'] is not None for r in group),legacy_joint_passes=sum(r['legacy_joint'] is True for r in group),confirmation_joint_passes=sum(r['confirmation_joint'] is True for r in group))
   for field in ('material_nrmse','otsu_iou','relative_residual','original_Phi','wall_seconds','cpu_seconds','LU','actual_complex_RHS'):
    values=[r[field] for r in group if r.get(field) is not None]
    summary[field+'_median']=float(np.median(values)) if values else None
    summary[field+'_known_total']=float(sum(values)) if values else None
    summary[field+'_missing']=len(group)-len(values)
   summaries.append(summary)
 Path(out).mkdir(parents=True)
 write(Path(out)/'SCORES.json',dict(batch=batch,rows=rows,summaries=summaries,physical_scoring_LU=0,physical_scoring_RHS=0,selection='No GT-based candidate/root selection; score fixed final original-Phi incumbent only',legacy_scope='Inherited data contracts remain separate from positive-public-sigma CN99 discrepancy',new_confirmation_count=6 if batch=='primary' else 0,fft_scope='Engineering duplicate observations, never new independent objects' if batch=='fft' else None))
 with (Path(out)/'SCORES.csv').open('x',newline='') as f:
  writer=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in rows for k in r)));writer.writeheader();writer.writerows(rows)
 if batch=='primary':feasible_witness(Path(out),rows,cases,fresh)
 plot_all(Path(out),rows,truths,pictures,private_labels)
 write(Path(out)/'SCORING_RECEIPT.json',dict(status='COMPLETE',seal=str(seal_path),batch=batch,wall_seconds=time.perf_counter()-started,cpu_seconds=time.process_time()-cpu,truth_access='INDEPENDENT_POST_SEAL_ONLY',new_physical_calls=0))
 print(json.dumps(dict(status='SCORED',batch=batch,rows=len(rows),output=str(out))))

def feasible_witness(out,rows,cases,fresh):
 """Post-seal offline feasible witness; reuses saved complete generation forward."""
 records=[]
 for meta in fresh:
  ident=meta['case_id'];case=cases[ident]
  with np.load(resolve(meta['truth_path']),allow_pickle=False) as z:
   m=z['m'].ravel();clean=z['Y_clean'];noise=z['noise']
  with np.load(resolve(case['observation_file']),allow_pickle=False) as z:Y=z['Y'];sigma=float(z['sigma'])
  residual=float(np.linalg.norm(clean-Y));relative=residual/float(np.linalg.norm(Y));regularization=prior(m,case)
  phi=.5*relative**2+regularization
  records.append(dict(case_id=ident,GT_original_Phi=phi,GT_relative_residual=relative,GT_residual_norm=residual,GT_prior=regularization,GT_whitened_residual=residual/sigma,stored_noise_consistency_error=float(np.linalg.norm((Y-clean)-noise)),feasible_material=bool(m.min()>=.02 and m.max()<=20),predictions=[dict(arm=r['arm'],status=r['status'],prediction_original_Phi=r['original_Phi'],prediction_relative_residual=r['relative_residual'],prediction_minus_GT_Phi=None if r['original_Phi'] is None else r['original_Phi']-phi,prediction_to_GT_relative_residual_ratio=None if r['relative_residual'] is None else r['relative_residual']/relative) for r in rows if r['case_id']==ident]))
 write(out/'GT_FEASIBLE_WITNESS.json',dict(status='OFFLINE_DIAGNOSTIC_POST_PREDICTION_SEAL',generation_model_identity='research.research_z_alg_r4_diversity.src.physics.Model.forward; disk-pulse scalar-TM, original retained LU/all12source RHS',operator_source='research/research_z_alg_r3_digits/public_data/operators.npz',source='Previously generated sealed m/Y_clean/noise; no regeneration',interpretation='A larger prediction Phi than this feasible witness shows failure to attain a known feasible objective; this does not establish uniqueness or a recovery theorem. No root/prediction selection or algorithm changes use this diagnostic.',extra_LU=0,extra_complex_RHS=0,records=records))

def plot_all(out,rows,truths,pictures,labels):
 import matplotlib;matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 for stratum in ('DEV9_DIGITS','DEV6_DIFFICULT','NEW6_CONFIRMATION'):
  group=[r for r in rows if r['stratum']==stratum]
  ids=list(dict.fromkeys(r['case_id'] for r in group));arms=list(dict.fromkeys(r['arm'] for r in group))
  if not ids:continue
  ids.sort(key=lambda ident:(labels[ident],ident));fig,axes=plt.subplots(len(ids),len(arms)+1,figsize=(3*(len(arms)+1),2.8*len(ids)),squeeze=False)
  for i,ident in enumerate(ids):
   axes[i,0].imshow(truths[ident].reshape(24,24),vmin=.02,vmax=12.1,cmap='viridis');axes[i,0].set_title(labels[ident]+' GT\n'+ident)
   for k,arm in enumerate(arms,1):
    row=next((r for r in group if r['case_id']==ident and r['arm']==arm),None);m=pictures.get((ident,arm))
    if m is not None:axes[i,k].imshow(m.reshape(24,24),vmin=.02,vmax=12.1,cmap='viridis')
    axes[i,k].set_title(arm+'\n'+(f"N={row['material_nrmse']:.3f} IoU={row['otsu_iou']:.2f}\nr={row['relative_residual']:.3g}" if m is not None else str(row['status'] if row else 'MISSING')))
   for ax in axes[i]:ax.set_xticks([]);ax.set_yticks([])
  fig.tight_layout();fig.savefig(out/(stratum+'_all_objects.png'),dpi=140);plt.close(fig)

def main():
 parser=argparse.ArgumentParser();sub=parser.add_subparsers(dest='command',required=True)
 p=sub.add_parser('seal');p.add_argument('--manifest',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
 p=sub.add_parser('score');p.add_argument('--seal',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--batch',choices=('primary','fft'),required=True);p.add_argument('--authorize-after-seal',action='store_true')
 a=parser.parse_args()
 if a.command=='seal':seal(a.manifest,a.output)
 else:score(a.seal,a.output,a.batch,a.authorize_after_seal)
if __name__=='__main__':main()
