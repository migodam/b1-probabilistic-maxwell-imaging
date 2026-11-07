"""Mechanical summary of the one explicitly authorized main48 scored CSV only."""
from pathlib import Path
from datetime import datetime,timezone
import csv,json,math,time,statistics
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[3]
SOURCE=ROOT/'delegated/twod/independent48/scored/SCORES_AND_FIXED_BUDGETS.csv'
METHODS=('cie_single','cie_dual','global_dual','fullG19','freezeA','contractedT')
REGIONS={'low':('true_low_cells','RMSE_true_low_region'),'high':('true_high_cells','RMSE_true_high_region'),'background':('true_background_cells','RMSE_true_background_region')}
def write(path,rows):
 keys=list(dict.fromkeys(k for r in rows for k in r))
 with path.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows({k:json.dumps(v) if isinstance(v,(list,dict)) else v for k,v in r.items()} for r in rows)
def val(x):return float(x) if x not in ('',None) else None
def basic(values):
 v=[x for x in values if x is not None]
 return dict(known=len(v),mean=statistics.mean(v) if v else None,median=statistics.median(v) if v else None,min=min(v) if v else None,max=max(v) if v else None)
def main():
 wall=time.perf_counter()
 with SOURCE.open(newline='') as f:
  reader=csv.DictReader(f);columns=list(reader.fieldnames);source=[r for r in reader if r['view']=='FINAL_CAP_ELIGIBLE']
 keys={(r['case_id'],r['method']) for r in source};cases={r['case_id'] for r in source}
 if len(source)!=288 or len(keys)!=288 or len(cases)!=48 or keys!={(c,m) for c in cases for m in METHODS}:raise ValueError('Exactly48x6 canonical FINAL_CAP_ELIGIBLE rows required; do not mix duplicateCPU480 view')
 mixed_cases={r['case_id'] for r in source if r['regime']=='mixed'}
 if len(mixed_cases)!=24:raise ValueError('24 independent mixed objects required')
 # No source column contains a low-phase mean. Never infer one from region RMSE.
 candidate_mean_fields=('low_region_mean_estimate','mean_estimate_true_low_region','mean_true_low_region','estimated_low_region_mean')
 mean_field=next((c for c in candidate_mean_fields if c in columns),None)
 records=[]
 for r in source:
  high=float(r['high_chi']);out=dict(case_id=r['case_id'],method=r['method'],family=r['family'],regime=r['regime'],high_chi=high,global_RMSE=val(r['RMSE']),global_NRMSE=val(r['NRMSE']),global_NRMSE_amplitude_denominator=high-.1,global_cell_denominator=576,material_success=r['material_success']=='True',global_NRMSE_units='dimensionless; RMSE/(high_chi-.1)',region_RMSE_units='absolute real contrast chi units',estimated_low_region_mean=val(r.get(mean_field)) if mean_field else None,estimated_low_region_mean_source=mean_field,estimated_low_region_mean_missing_reason=None if mean_field else 'NOT_RETAINED_IN_SOURCE_SCORED_CSV; no GT reconstruction')
  for name,(countfield,rmsefield) in REGIONS.items():
   out[name+'_cells']=int(r[countfield]) if r.get(countfield) not in ('',None) else None;out[name+'_region_RMSE']=val(r.get(rmsefield))
  records.append(out)
 mixed=[r for r in records if r['regime']=='mixed'];write(HERE/'MIXED24_BY_METHOD_OBJECT144.csv',mixed);write(HERE/'ALL48_FINAL_BY_METHOD_OBJECT288.csv',records)
 region_summary=[]
 for method in METHODS:
  rr=[r for r in mixed if r['method']==method];nr=basic([r['global_NRMSE'] for r in rr]);out=dict(method=method,independent_mixed_objects=24,rows=len(rr),global_NRMSE_mean=nr['mean'],global_NRMSE_median=nr['median'],material_success=sum(r['material_success'] for r in rr),estimated_low_region_mean=None,estimated_low_region_mean_known=sum(r['estimated_low_region_mean'] is not None for r in rr),estimated_low_region_mean_missing_reason='NOT_RETAINED_IN_SOURCE_SCORED_CSV' if mean_field is None else 'See per-object fields',units='RegionRMSE:absolute chi; globalNRMSE:dimensionless; object means/medians weight each of24materials equally')
  for name in REGIONS:
   stats=basic([r[name+'_region_RMSE'] for r in rr]);valid=[r for r in rr if r[name+'_region_RMSE'] is not None and r[name+'_cells'] is not None and r[name+'_cells']>0];cells=sum(r[name+'_cells'] for r in valid)
   out.update({name+'_objects_with_metric':stats['known'],name+'_object_mean_RMSE':stats['mean'],name+'_object_median_RMSE':stats['median'],name+'_object_min_RMSE':stats['min'],name+'_object_max_RMSE':stats['max'],name+'_total_true_cells':cells,name+'_pooled_cell_RMSE':math.sqrt(sum(r[name+'_cells']*r[name+'_region_RMSE']**2 for r in valid)/cells) if cells else None})
  region_summary.append(out)
 write(HERE/'MIXED24_REGION_SUMMARY.csv',region_summary)
 groups=[]
 for method in METHODS:
  for regime in ('high','mixed'):
   for high in (8.,12.,18.):
    rr=[r for r in records if r['method']==method and r['regime']==regime and r['high_chi']==high]
    if len(rr)!=8:raise ValueError('Each high/mixed x contrast group must have8 independent objects per method')
    n=basic([r['global_NRMSE'] for r in rr]);groups.append(dict(method=method,regime=regime,high_chi=high,independent_objects=8,NRMSE_amplitude_denominator=high-.1,NRMSE_mean=n['mean'],NRMSE_median=n['median'],NRMSE_min=n['min'],NRMSE_max=n['max'],material_success=sum(r['material_success'] for r in rr),case_ids=sorted(r['case_id'] for r in rr),sample_scope='Same8materials reused across6methods, not48independent observations per group'))
 write(HERE/'REGIME_CONTRAST_GROUP_SUMMARY36.csv',groups)
 text=['# Main48:24 mixed-object phase summary','', 'Only FINAL_CAP_ELIGIBLE rows from the already-authorized scored CSV are used. CPU_480 rows are not pooled or duplicated. There are24 independent mixed materials,144 method/object rows, plus24high materials. No GT array, prediction array, private file, forward, scoring or bootstrap is used.','', '|Method|Mixed materialPASS/24|Median globalNRMSE|Low-regionRMSE mean / median|High-regionRMSE mean / median|BackgroundRMSE mean / median|Estimated low mean|','|---|---:|---:|---:|---:|---:|---|']
 for r in region_summary:text.append(f"|{r['method']}|{r['material_success']}|{r['global_NRMSE_median']:.6f}|{r['low_object_mean_RMSE']:.6f} / {r['low_object_median_RMSE']:.6f}|{r['high_object_mean_RMSE']:.6f} / {r['high_object_median_RMSE']:.6f}|{r['background_object_mean_RMSE']:.6f} / {r['background_object_median_RMSE']:.6f}|MISSING|")
 text+=['', 'RegionRMSE is in absolute chi units. The low truth phase ischi2.5; each per-object regional RMSE uses only that object’s truth-region cell count from the source table. Summary means/medians give each material equal weight. The CSV also reports a separately labelled pooled-cell RMSE=sqrt(sum(n_cells*RMSE_region²)/sum(n_cells)); those cells are not independent scientific samples.','', 'GlobalNRMSE is dimensionless: globalRMSE over576cells divided byhigh_chi−.1, i.e.7.9/11.9/17.9 for the three declared contrasts. Thus global and regional numbers have different denominators/units. Regional errors do not change the inherited material gate.','', 'Estimated low-phase mean is absent from the original scored CSV. It remainsnull/MISSING for all144mixed rows; it cannot be reconstructed from RMSE alone. No GT/prediction arrays were opened to supply it, and the earlier16-object low means must not be substituted into this48-object cohort.','', 'Each high/mixed×high_chi8/12/18 group has8 independent objects per method. Continuous NRMSE and material-success counts for all36method/group rows are in REGIME_CONTRAST_GROUP_SUMMARY36.csv; all288per-object rows and the144mixed rows remain available. This is mechanical descriptive reporting only, with no new threshold or inferential claim.']
 text+=['','## High/mixed x contrast groups','', 'Each cell is mean dimensionless NRMSE (material successes/8); each row has8independent materials reused across methods.','', '|Regime / highchi|Objects per method|CIE|Dual|Global|G19|freezeA|contractedT|','|---|---:|---:|---:|---:|---:|---:|---:|']
 for regime in ('high','mixed'):
  for high in (8.,12.,18.):
   gg=[next(x for x in groups if x['method']==m and x['regime']==regime and x['high_chi']==high) for m in METHODS]
   text.append('|'+regime+' / '+str(int(high))+'|8|'+'|'.join(f"{x['NRMSE_mean']:.6f} ({x['material_success']}/8)" for x in gg)+'|')
 (HERE/'SUMMARY.md').write_text('\n'.join(text)+'\n')
 receipt=dict(status='COMPLETE_MECHANICAL_MIXED_PHASE_SUMMARY',created_utc=datetime.now(timezone.utc).isoformat(),explicit_input_whitelist=[str(SOURCE)],source_view='FINAL_CAP_ELIGIBLE',source_rows=288,mixed_independent_materials=24,mixed_method_object_rows=144,regime_contrast_groups=6,method_group_rows=36,estimated_low_mean_source_field=mean_field,estimated_low_mean_available=mean_field is not None,missing_not_reconstructed=True,new_forward=0,new_adjoint=0,new_scoring=0,new_bootstrap=0,GT_or_prediction_arrays_read=0,private_or_unsealed_paths_accessed=0,CPU_seconds_near_write=time.process_time(),wall_seconds=time.perf_counter()-wall)
 (HERE/'RECEIPT.json').write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n');print(json.dumps(region_summary))
if __name__=='__main__':main()
