"""Explicitly authorized post64seal scoring only; all model error uses cached cleanY."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
from datetime import datetime,timezone
import argparse,json,time,csv
import numpy as np
from scipy.stats import chi2
HERE=Path(__file__).absolute().parent;ROOT=HERE.parents[2];REPO=ROOT.parents[1]
def read(p):return json.loads(Path(p).read_text())
def dump(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def arrays(p):
 with np.load(p,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}
def csvwrite(p,rows):
 keys=list(dict.fromkeys(k for r in rows for k in r))
 with Path(p).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows({k:json.dumps(v) if isinstance(v,(list,dict)) else v for k,v in r.items()} for r in rows)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--public-seal',type=Path,default=HERE/'public_seal');ap.add_argument('--output',type=Path,default=HERE/'scored');ap.add_argument('--authorize',required=True);a=ap.parse_args();wall=time.perf_counter()
 if a.authorize!='ROOT_EXPLICIT_POST64SEAL_SCORE':raise PermissionError('Explicit root post64seal scoring authorization required; private unopened')
 seal=read(a.public_seal/'RECEIPT.json');rule=read(a.public_seal/'SCORING_CONFIG_V2.json');selection=read(a.public_seal/'SELECTIONS.json')
 if seal['status']!='ALL64_PUBLIC_PAIRED_GRID_CHOICES_SEALED_BEFORE_GT' or seal['tasks']!=64 or seal['independent_materials']!=8 or seal['truth_access']!='NONE':raise RuntimeError('All64 public seal required')
 if rule['NRMSE']['denominator']!=11.9 or rule['NRMSE']['cells']!=576 or rule['support']['threshold']!=1. or rule['success']!={'NRMSE_max':.15,'IoU_min':.7,'material_box':[.02,20.],'box_tolerance':1e-12}:raise ValueError('Frozen score definition differs')
 sources=selection['sources'];cases=sorted({e['task']['case_id'] for e in sources})
 if len(sources)!=64 or len(cases)!=8 or len({e['task']['task_id'] for e in sources})!=64:raise ValueError('Original64 identities/8materials required')
 a.output.mkdir(parents=True,exist_ok=False);dump(a.output/'PUBLIC_GATE_BEFORE_PRIVATE.json',dict(status='PASS_ALL64_PUBLIC_SEAL',truth_opened=False,created_utc=datetime.now(timezone.utc).isoformat(),selection_seal_finished_utc=seal['finished_utc']))
 # First private read is below, after explicit root authorization and public gate.
 private=ROOT/'sealed_truth/crossgrid8_paired';truth={cid:arrays(private/(cid+'.npz')) for cid in cases};meta={r['case_id']:r for r in read(private/'generator_metadata.json')['records']};selected=arrays(a.public_seal/'SELECTED_MATERIALS.npz');rows=[];cache={}
 for r in selection['rows']:
  cid=r['case_id'];g=r['grid_condition'];gt=truth[cid];s=r['selected'];row={k:r[k] for k in ('case_id','grid_condition','method','task_id','view','pair_id','side','budget_B','task_status','actual_complete_child_CPU_seconds')};row.update(family=meta[cid]['family'],regime=meta[cid]['regime'],available=s is not None,material_success=False)
  if s is not None:
   m=selected['m'][s['material_index']];target=gt['m'];error=m-target
   if m.shape!=(576,) or np.iscomplexobj(m) or not np.isfinite(m).all():raise ValueError('576 finite real materials required')
   rmse=float(np.linalg.norm(error)/24);nrmse=rmse/11.9;union=(m>1)|(target>1);iou=float(np.count_nonzero((m>1)&(target>1))/np.count_nonzero(union)) if union.any() else 1.;legal=bool(m.min()>=.02-1e-12 and m.max()<=20+1e-12)
   row.update(NRMSE=nrmse,RMSE=rmse,fixed_support_IoU=iou,material_legal=legal,material_success=bool(legal and nrmse<=.15 and iou>=.7),original_Phi=s['original_Phi'],paid_CPU_seconds=s['paid_CPU_seconds'],paid_LU=s['paid_LU'],paid_complex_RHS=s['paid_complex_RHS'])
   masks=[gt[k].astype(bool) for k in ('pure_low','pure_high','pure_background','partial_volume_boundary')]
   if not np.all(np.sum(masks,axis=0)==1):raise ValueError('Fine-defined phase partition must cover576cells exactly')
   for name,mask in zip(('pure_low','pure_high','pure_background','partial_volume_boundary'),masks):row[name+'_cells']=int(mask.sum());row[name+'_RMSE']=float(np.sqrt(np.mean(error[mask]**2))) if mask.any() else None
   if r['view']=='FINAL_CAP_ELIGIBLE':
    key=(cid,g)
    if key not in cache:
     path=Path(r['public_files'][-1].replace('\\','/'));obs=arrays(path if path.is_absolute() else REPO/path);cache[key]=obs
    obs=cache[key];yp=selected['Y_prediction'][r['Y_prediction_index']];sigma=float(obs['sigma']);res=float(np.linalg.norm(yp-obs['Y'])/sigma);radius=float(np.sqrt(chi2.ppf(.99,2*obs['Y'].size)/2))
    row.update(residual_over_sigma=res,CN99_radius=radius,CN99_pass=bool(res<=radius),noise_over_sigma=float(np.linalg.norm(obs['Y']-gt['Y_clean_'+g])/sigma),data_scope='Coarse model prediction vs this condition observations, raw measurement sigma; CN failure allowed')
  rows.append(row)
 csvwrite(a.output/'ALL_SELECTION_SCORES.csv',rows);final=[r for r in rows if r['view']=='FINAL_CAP_ELIGIBLE'];csvwrite(a.output/'FINAL64_SCORES.csv',final)
 model=[]
 for cid in cases:
  gt=truth[cid];diff=gt['Y_clean_matched24']-gt['Y_clean_cross48to24'];norm=float(np.linalg.norm(diff));model.append(dict(case_id=cid,clean_model_error_Frobenius=norm,relative_to_fine_clean=norm/float(np.linalg.norm(gt['Y_clean_cross48to24'])),over_cross_raw_sigma=norm/float(gt['sigma_cross48to24']),over_matched_raw_sigma=norm/float(gt['sigma_matched24']),new_forward=0,scope='Cached cleanY, same projected target; fine/coarse disk-model mismatch, not PDE error bound'))
 csvwrite(a.output/'CACHED_MODEL_MISMATCH_8.csv',model)
 differences=[]
 for p in selection['pairs']:
  l,r=rows[p['left_selection_row']],rows[p['right_selection_row']];differences.append(dict(case_id=p['case_id'],pair_id=p['pair_id'],common_B=p['common_B'],available=p['both_available'],NRMSE_right_minus_left=r['NRMSE']-l['NRMSE'] if p['both_available'] else None,left_material_success=l['material_success'],right_material_success=r['material_success'],independence='One paired material, grid/methods not independent objects',spend_scope='Both original children fully charged, commonB trajectory readout only'))
 csvwrite(a.output/'ACTUAL_SPEND_PAIRED_DIFFERENCES.csv',differences)
 by={(r['case_id'],r['grid_condition'],r['method']):r for r in final};grid=[]
 for cid in cases:
  for m in rule['methods']:
   l,r=by[cid,'matched24',m],by[cid,'cross48to24',m];grid.append(dict(case_id=cid,method=m,available=l['available'] and r['available'],final_NRMSE_cross_minus_matched=r['NRMSE']-l['NRMSE'] if l['available'] and r['available'] else None,matched_material_success=l['material_success'],cross_material_success=r['material_success']))
 csvwrite(a.output/'FINAL_PAIRED_GRID_DIFFERENCES.csv',grid)
 summary=[]
 for condition in rule['grid_conditions']:
  for m in rule['methods']:
   rr=[r for r in final if r['grid_condition']==condition and r['method']==m];finite=[r['NRMSE'] for r in rr if r['available']];summary.append(dict(grid_condition=condition,method=m,independent_materials=8,available=len(finite),material_success=sum(r['material_success'] for r in rr),mean_NRMSE=float(np.mean(finite)) if finite else None,median_NRMSE=float(np.median(finite)) if finite else None))
 dump(a.output/'SUMMARY.json',summary);dump(a.output/'RECEIPT.json',dict(status='SCORED_POST64_PUBLIC_SEAL',independent_materials=8,grid_conditions=2,tasks=64,new_FW=0,new_adjoint=0,new_inverse=0,GT_choice=False,CPU_seconds_near_write=time.process_time(),wall_seconds=time.perf_counter()-wall,statistical_scope='Descriptive8 material pairs; no16/64independent-object or continuous-PDE claim',cost_scope='Future external complete child CPU replaces nested self timer'))
 print(json.dumps(summary))
if __name__=='__main__':main()
