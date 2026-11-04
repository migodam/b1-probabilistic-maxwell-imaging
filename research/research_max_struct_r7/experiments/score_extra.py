"""Post-seal offline extra runs: init18, disclosed dev20, old digit9.
These repeated/disclosed objects never enlarge the R7 new30 denominator.
No physics imports or new forward calls; all endpoints fixed before GT access.
"""
from pathlib import Path
from datetime import datetime,timezone
import argparse,json,time
import numpy as np
try:from .score_primary import load,write,csvwrite,resolve,stamp,authorize,verify_seal,material_metrics,noise_metrics,public_observation,TERMINAL
except ImportError:from score_primary import load,write,csvwrite,resolve,stamp,authorize,verify_seal,material_metrics,noise_metrics,public_observation,TERMINAL
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
DEFAULT_MANIFESTS=['configs/initialization_manifest.json','configs/aplus_all_manifest.json','delegated/aplus/DEVELOPMENT_MANIFEST.json','delegated/aplus/GUARDED_DEVELOPMENT_MANIFEST.json','delegated/aplus/ARC_DEVELOPMENT_MANIFEST.json']

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--primary-seal',required=True);p.add_argument('--primary-scores',required=True);p.add_argument('--authorization',required=True);p.add_argument('--output',required=True);p.add_argument('--manifest',action='append');a=p.parse_args()
 authorize(a.authorization);primary=verify_seal(a.primary_seal);start,cpu=time.perf_counter(),time.process_time();out=Path(a.output)
 if out.exists():raise FileExistsError('Write-once extra-score output exists')
 primaryids={r['task']['task_id'] for r in primary['tasks']};seen={};stamps=[]
 for pth in a.manifest or [str(ROOT/m) for m in DEFAULT_MANIFESTS]:
  manifest=resolve(pth);stamps.append(stamp(manifest));data=load(manifest)
  if data.get('config'):stamps.append(stamp(resolve(data['config'])))
  for t in data['tasks']:
   if t['task_id'] in primaryids:continue
   directory=resolve(t['output_dir']);key=str(directory)
   if key in seen:
    if seen[key]['task_id']!=t['task_id']:raise RuntimeError('Conflicting endpoint tasks share directory')
    continue
   seen[key]=t
   if t.get('config'):stamps.append(stamp(resolve(t['config'])))
 tasks=list(seen.values())
 if a.manifest is None and len(tasks)!=47:raise RuntimeError('Expected47 extra tasks after primary exclusion and directory deduplication')
 # Seal all expected47 extras publicly before private mapping/GT opens.
 entries=[]
 for t in tasks:
  directory=resolve(t['output_dir']);r=load(directory/'RECEIPT.json')
  if r['status'] not in TERMINAL:raise RuntimeError('Nonterminal extra task '+t['task_id'])
  if r['status']=='SEALED_PUBLIC_PREDICTION' and not(directory/'prediction.npz').is_file():raise RuntimeError('Missing extra prediction')
  entries.append(dict(task=t,status=r['status']));stamps.extend(stamp(f) for f in sorted(directory.iterdir()) if f.is_file());stamps.extend(stamp(resolve(f)) for f in t['public_files'])
 old=REPO/'research/research_max_struct_r6';newcases={r['case_id']:r for r in load(ROOT/'data_public/cases.json')['records']};oldcases={r['case_id']:r for r in load(old/'data_public/cases.json')['records']}
 fresh=load(ROOT/'sealed_truth/mapping_PRIVATE.json')['records'];oldmap=load(old/'sealed_truth/scoring/historical_mapping_PRIVATE.json')+load(old/'sealed_truth/confirmation/mapping_PRIVATE.json')['records'];mapping={r['case_id']:r for r in fresh+oldmap}
 primaryrows=load(Path(a.primary_scores)/'SCORES.json')['rows'];primarylookup={(r['case_id'],r['arm']):r for r in primaryrows};oldrows=load(old/'RESULTS/scoring_primary/SCORES.json')['rows']
 rows=[];pictures={};truths={}
 for e in entries:
  t=e['task'];ident=t['case_id'];case=newcases.get(ident,oldcases.get(ident));meta=mapping[ident];directory=resolve(t['output_dir']);rec=load(directory/'RECEIPT.json');cost=rec.get('physical_model_cost') or {};corr=rec.get('corrector') or {}
  if t.get('initial_material') is not None:cohort='REPEATED_NEW_OBJECT_INITIALIZATION_DIAGNOSTIC'
  elif t.get('development_mode') is not None:cohort='DISCLOSED_FOUR_OBJECT_APLUS_DEVELOPMENT20'
  else:cohort='DISCLOSED_OLD_DIGIT9_REGRESSION'
  if ident not in truths:
   with np.load(resolve(meta['truth_path']),allow_pickle=False) as z:truths[ident]=(z['m'].ravel().copy(),z['mask'].ravel().astype(bool))
  truth,support=truths[ident];label=meta.get('parameters',{}).get('family',meta.get('shape',meta.get('family',ident)))
  row=dict(task_id=t['task_id'],case_id=ident,arm=t['arm'],initial_material=t.get('initial_material',.1),development_mode=t.get('development_mode'),diagnostic_cohort=cohort,status=e['status'],denominator_included=True,new_independent_R7_objects=0,label=label,material_image_success=False,material_nrmse=None,otsu_iou=None,relative_residual=None,original_Phi=rec.get('physical_objective'),cpu_seconds=rec.get('cpu_seconds'),wall_seconds=rec.get('wall_seconds'),LU=cost.get('LU'),actual_complex_RHS=cost.get('actual_complex_RHS'),corrector_stop=corr.get('stop'),corrector_success=corr.get('success'),projected_gradient_inf=None)
  if e['status']=='SEALED_PUBLIC_PREDICTION':
   with np.load(directory/'prediction.npz',allow_pickle=False) as z:pred={k:z[k].copy() for k in z.files}
   row.update(material_metrics(pred['m'],truth,support,case['amplitude']));row.update(noise_metrics(pred['Y_prediction'],public_observation(t),case.get('inverse_crime',True)))
   if 'gradient' in pred:
    m=pred['m'].ravel();g=pred['gradient'].ravel();row['projected_gradient_inf']=float(np.linalg.norm(m-np.clip(m-g,.02,20),np.inf))
   pictures[t['task_id']]=pred['m'].ravel()
  if cohort=='REPEATED_NEW_OBJECT_INITIALIZATION_DIAGNOSTIC':
   base=primarylookup[(ident,t['arm'])];row.update(primary_host_initialization=.1,primary_host_NRMSE=base.get('material_nrmse'),primary_host_material_success=base['material_image_success'],primary_host_original_Phi=base.get('original_Phi'),NRMSE_change_init1_minus_host=None if row['material_nrmse'] is None or base.get('material_nrmse') is None else row['material_nrmse']-base['material_nrmse'],Phi_change_init1_minus_host=None if row['original_Phi'] is None or base.get('original_Phi') is None else row['original_Phi']-base['original_Phi'],selection_scope='Both endpoints retained; no choose-best-initialization material score is promoted to primary')
  if cohort=='DISCLOSED_OLD_DIGIT9_REGRESSION':
   refs=[r for r in oldrows if r['case_id']==ident];row['old_R6_dense_references']=[{k:r.get(k) for k in ['arm','material_nrmse','otsu_iou','legacy_joint','relative_residual','original_Phi','cpu_seconds','LU','actual_complex_RHS']} for r in refs];row['comparison_limit']='New Aplus uses FFT actions versus old R6 dense historical endpoints; readout and backend/roundoff/path confounds prevent a pure single-factor causal interpretation'
  rows.append(row)
 cohorts=[]
 for cohort in sorted({r['diagnostic_cohort'] for r in rows}):
  group=[r for r in rows if r['diagnostic_cohort']==cohort];keys=sorted({(r['arm'],r['development_mode'] or '') for r in group})
  for arm,mode in keys:
   g=[r for r in group if r['arm']==arm and (r['development_mode'] or '')==mode];valid=[r['material_nrmse'] for r in g if r['material_nrmse'] is not None]
   cohorts.append(dict(diagnostic_cohort=cohort,arm=arm,development_mode=mode or None,runs=len(g),distinct_existing_objects=len({r['case_id'] for r in g}),new_independent_R7_objects=0,scored=len(valid),material_successes=sum(r['material_image_success'] for r in g),median_material_NRMSE=float(np.median(valid)) if valid else None,mean_material_NRMSE=float(np.mean(valid)) if valid else None,total_cpu_seconds=float(sum(r['cpu_seconds'] or 0 for r in g)),total_LU=int(sum(r['LU'] or 0 for r in g))))
 out.mkdir(parents=True);write(out/'EXTRA_PUBLIC_SEAL.json',dict(created_utc=datetime.now(timezone.utc).isoformat(),tasks=len(tasks),files=stamps,truth_access_at_seal='NONE',hash_checks='NONE'));write(out/'EXTRA_SCORES.json',dict(rows=rows,summaries=cohorts,primary_independent_object_denominator=30,new_independent_objects_added=0,material_success_rule='same frozen amplitudeNRMSE<=.15, OtsuIoU>=.7 and legalbox',selection='Score all fixed endpoints, never use GT to revise frozen Aplus or promote best repeated initialization'));csvwrite(out/'EXTRA_SCORES.csv',rows);write(out/'EXTRA_SUMMARY.json',cohorts)
 plot_extra(rows,pictures,truths,ROOT/'FIGURES/primary')
 lines=['# Offline extra-run recovery diagnostics','', 'These repeated/disclosed objects add zero independent objects to the R7 new30 denominator. All GT scoring follows complete primary sealing.','','|Cohort|Arm/readout|Runs|Existing objects|Material successes|Median NRMSE|CPU s|LU|','|---|---|---:|---:|---:|---:|---:|---:|']
 for r in cohorts:lines.append('|%s|%s|%d|%d|%d|%s|%.2f|%d|'%(r['diagnostic_cohort'],r['development_mode'] or r['arm'],r['runs'],r['distinct_existing_objects'],r['material_successes'],'NA' if r['median_material_NRMSE'] is None else '%.4f'%r['median_material_NRMSE'],r['total_cpu_seconds'],r['total_LU']))
 lines+=['','Init18 repeats six existing new30 objects at m=1 against their fixed primary m=.1 endpoints for three arms. Both predictions remain visible; no GT-selected multi-start replacement enters the primary table.','','Development20 is five readout variants on four disclosed historical objects. Recovery is scored only after the Aplus publicPhi/cost-based choice was frozen; post-hoc GT does not authorize algorithm reselection.','','The9 Aplus digit cases are disclosed historical regression. Old R6 dense endpoints are copied as labeled references without rerun. FFT action roundoff, readout changes, and resulting long-trajectory path changes prevent interpreting the new-vs-old difference as a pure readout causal effect. These old digits are not new confirmation.','','Scoring reads complete576-variable endpoint materials and paid Y_prediction only;0new physical forward/adjoint calls.']
 (out/'EXTRA_ANALYSIS.md').write_text('\n'.join(lines)+'\n');write(out/'SCORE_RECEIPT.json',dict(status='COMPLETE',created_utc=datetime.now(timezone.utc).isoformat(),authorization='ALL_PRIMARY_SEALED',primary_seal=a.primary_seal,rows=len(rows),new_independent_objects=0,new_physical_calls=0,new_physical_LU=0,new_physical_RHS=0,wall_seconds=time.perf_counter()-start,cpu_seconds=time.process_time()-cpu))
 print(json.dumps(dict(status='EXTRA_SCORED',runs=len(rows),out=str(out))))

def plot_extra(rows,pictures,truths,folder):
 import matplotlib;matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 folder=Path(folder);folder.mkdir(parents=True,exist_ok=True)
 for cohort in sorted({r['diagnostic_cohort'] for r in rows}):
  group=[r for r in rows if r['diagnostic_cohort']==cohort];ids=sorted({r['case_id'] for r in group});variants=list(dict.fromkeys((r['development_mode'] or r['arm']) for r in group));fig,axes=plt.subplots(len(ids),len(variants)+1,figsize=(2.7*(len(variants)+1),2.6*len(ids)),squeeze=False,layout='constrained')
  for i,ident in enumerate(ids):
   g=[r for r in group if r['case_id']==ident];truth=truths[ident][0];vmax=max([float(np.max(truth))]+[float(np.max(pictures[r['task_id']])) for r in g if r['task_id'] in pictures]);axes[i,0].imshow(truth.reshape(24,24),vmin=.02,vmax=vmax,cmap='viridis');axes[i,0].set_title(ident+'\n'+g[0]['label']+' GT',fontsize=8)
   for j,mode in enumerate(variants,1):
    r=next((r for r in g if (r['development_mode'] or r['arm'])==mode),None)
    if r and r['task_id'] in pictures:axes[i,j].imshow(pictures[r['task_id']].reshape(24,24),vmin=.02,vmax=vmax,cmap='viridis');title=mode+'\nNRMSE %.3f | IoU %.2f'%(r['material_nrmse'],r['otsu_iou'])
    else:title=mode+'\n'+(r['status'] if r else 'NOT_RUN')
    axes[i,j].set_title(title,fontsize=8)
   for ax in axes[i]:ax.set_xticks([]);ax.set_yticks([])
  titles={'DISCLOSED_OLD_DIGIT9_REGRESSION':'Historical digit9 regression\nDisclosed objects; no new confirmation n','REPEATED_NEW_OBJECT_INITIALIZATION_DIAGNOSTIC':'Initialization diagnostic: six existing objects at m=1\nPrimary m=.1 predictions retained; no new objects','DISCLOSED_FOUR_OBJECT_APLUS_DEVELOPMENT20':'Aplus readout development: four disclosed objects, five variants\nFrozen before new30 GT scoring'}
  fig.suptitle(titles[cohort],fontsize=11);fig.savefig(folder/(cohort.lower()+'.png'),dpi=150);plt.close(fig)
if __name__=='__main__':main()
