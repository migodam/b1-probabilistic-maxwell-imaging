"""Offline paired-object statistics and actual-cost grading; no physics imports."""
from pathlib import Path
from datetime import datetime,timezone
import argparse,json,time
import numpy as np
try:from .score_primary import load,write,csvwrite,ARMS
except ImportError:from score_primary import load,write,csvwrite,ARMS
ROOT=Path(__file__).resolve().parents[1]
COMPARISONS=[('a_material','global_dual','PRESPECIFIED_A_VS_GLOBAL_DUAL'),('a_plus','a_material','PRESPECIFIED_APLUS_VS_A'),('a_material','cie_dual','SECONDARY_A_VS_SOURCEWISE_DUAL')]

def descriptive(values):
 a=np.asarray([v for v in values if v is not None and np.isfinite(v)],float)
 return dict(n=int(len(a)),mean=float(np.mean(a)) if len(a) else None,median=float(np.median(a)) if len(a) else None,q10=float(np.quantile(a,.1)) if len(a) else None,q25=float(np.quantile(a,.25)) if len(a) else None,q75=float(np.quantile(a,.75)) if len(a) else None,q90=float(np.quantile(a,.9)) if len(a) else None,min=float(np.min(a)) if len(a) else None,max=float(np.max(a)) if len(a) else None,total=float(np.sum(a)) if len(a) else None)

def aggregate(rows,groupkeys):
 groups={}
 for r in rows:groups.setdefault(tuple(r[k] for k in groupkeys),[]).append(r)
 out=[]
 for key,group in groups.items():
  row=dict(zip(groupkeys,key));row.update(objects=len(group),sealed=sum(r['status']=='SEALED_PUBLIC_PREDICTION' for r in group),material_successes=sum(r.get('material_image_success') is True for r in group),same_grid_objects=sum(r['inverse_crime'] for r in group),same_grid_joint_successes=sum(r.get('same_grid_material_noise_success') is True for r in group),cross_grid_objects=sum(not r['inverse_crime'] for r in group),cross_grid_material_and_raw_noise_passes=sum(r.get('cross_grid_material_and_raw_noise_pass') is True for r in group))
  for field in ['material_nrmse','otsu_iou','relative_residual','original_Phi','cpu_seconds','wall_seconds','LU','actual_complex_RHS','auxiliary_cpu_seconds']:row[field]=descriptive([r.get(field) for r in group])
  out.append(row)
 return out

def bootstrap_pair(rows,left,right,label,rng,reps=10000):
 by={(r['case_id'],r['arm']):r for r in rows};ids=sorted({r['case_id'] for r in rows});pairs=[(by.get((i,left)),by.get((i,right))) for i in ids]
 complete=[(a,b) for a,b in pairs if a is not None and b is not None];n=len(complete)
 result=dict(comparison=label,left=left,right=right,object_denominator=len(ids),paired_rows=n,confidence=.95,bootstrap='paired object percentile bootstrap; each resample uses same object indices for both arms; acquisitions are attributes, not extra samples',repetitions=reps,multiple_comparison_adjustment='NONE; two prespecified contrasts plus separately labeled secondary contrast',metrics={})
 # Success includes failed task as false; continuous metrics use explicitly named complete pairs.
 for field in ['material_image_success','same_grid_material_noise_success','material_nrmse','otsu_iou','relative_residual','cpu_seconds','LU','actual_complex_RHS']:
  valid=[(a,b) for a,b in complete if a.get(field) is not None and b.get(field) is not None and np.isfinite(a[field]) and np.isfinite(b[field])]
  if not valid:continue
  x=np.array([float(a[field]) for a,b in valid]);y=np.array([float(b[field]) for a,b in valid]);d=x-y;idx=rng.integers(0,len(d),size=(reps,len(d)))
  means=np.mean(d[idx],axis=1);medians=np.median(d[idx],axis=1)
  result['metrics'][field]=dict(pairs=len(d),mean_difference_left_minus_right=float(np.mean(d)),mean_difference_CI95=np.quantile(means,[.025,.975]).tolist(),median_paired_difference=float(np.median(d)),median_paired_difference_CI95=np.quantile(medians,[.025,.975]).tolist(),left_better=int(np.sum(d<0)) if field!='otsu_iou' and 'success' not in field else int(np.sum(d>0)),ties=int(np.sum(d==0)),direction='higher better' if field=='otsu_iou' or 'success' in field else 'lower better')
 return result

def matched_cost(rows,checkpoints):
 final={(r['case_id'],r['arm']):r for r in rows};groups={}
 for r in checkpoints:groups.setdefault((r['case_id'],r['arm']),[]).append(r)
 result=[]
 for ident in sorted({r['case_id'] for r in rows}):
  for cost in ['CPU','LU','complex_RHS']:
   available={arm:groups.get((ident,arm),[]) for arm in ARMS}
   if any(not values for values in available.values()):
    result.append(dict(case_id=ident,cost_axis=cost,status='UNAVAILABLE_ARM_CHECKPOINTS',denominator_included=True));continue
   maxbudget=min(max(r[cost] for r in values) for values in available.values())
   for fraction in [.25,.5,.75,1.]:
    budget=fraction*maxbudget
    for arm,values in available.items():
     eligible=[r for r in values if r[cost]<=budget+1e-10];base=final[(ident,arm)]
     row=dict(case_id=ident,arm=arm,condition=base['condition'],family=base['family'],inverse_crime=base['inverse_crime'],cost_axis=cost,fraction_of_common_max=fraction,common_maxbudget=maxbudget,target_budget=budget,denominator_included=True)
     if not eligible:row.update(status='NO_CHECKPOINT_BELOW_BUDGET',material_nrmse=None,otsu_iou=None,material_image_success=False)
     else:
      chosen=eligible[-1];row.update(status='SCORED_LAST_SAVED_AT_OR_BELOW_BUDGET',actual_checkpoint_cost=chosen[cost],unused_budget=budget-chosen[cost],checkpoint_index=chosen['index'],CPU=chosen['CPU'],LU=chosen['LU'],complex_RHS=chosen['complex_RHS'],original_Phi=chosen['original_Phi'],material_nrmse=chosen['material_nrmse'],otsu_iou=chosen['otsu_iou'],RMSE=chosen['RMSE'],material_image_success=chosen['material_image_success'],material_legal=chosen['material_legal'])
     result.append(row)
 return result

def pairwise_frontiers(rows,checkpoints,rng,reps=10000,comparisons=None):
 final={(r['case_id'],r['arm']):r for r in rows};groups={}
 for r in checkpoints:groups.setdefault((r['case_id'],r['arm']),[]).append(r)
 records=[];summaries=[];ids=sorted({r['case_id'] for r in rows})
 for left,right,label in (COMPARISONS if comparisons is None else comparisons):
  for axis in ['CPU','LU','complex_RHS']:
   for fraction in np.linspace(0.,1.,11):
    paired=[]
    for ident in ids:
     av=groups.get((ident,left),[]);bv=groups.get((ident,right),[])
     if not av or not bv:
      records.append(dict(case_id=ident,left=left,right=right,comparison=label,cost_axis=axis,window_fraction=float(fraction),status='UNAVAILABLE_CHECKPOINTS'));continue
     low=max(av[0][axis],bv[0][axis]);high=min(av[-1][axis],bv[-1][axis])
     if high<low:
      records.append(dict(case_id=ident,left=left,right=right,comparison=label,cost_axis=axis,window_fraction=float(fraction),status='NO_OVERLAPPING_OBSERVED_COST_RANGE',lower_cost=low,upper_cost=high));continue
     budget=low+fraction*(high-low);aa=[r for r in av if r[axis]<=budget+1e-10];bb=[r for r in bv if r[axis]<=budget+1e-10]
     a,b=aa[-1],bb[-1];record=dict(case_id=ident,left=left,right=right,comparison=label,cost_axis=axis,window_fraction=float(fraction),status='PAIRED_LAST_SAVED_WITHIN_OBSERVED_RANGE',lower_cost=low,upper_cost=high,target_budget=float(budget),left_checkpoint_cost=a[axis],right_checkpoint_cost=b[axis],left_unused_budget=float(budget-a[axis]),right_unused_budget=float(budget-b[axis]),left_NRMSE=a['material_nrmse'],right_NRMSE=b['material_nrmse'],NRMSE_difference_left_minus_right=a['material_nrmse']-b['material_nrmse'],left_IoU=a['otsu_iou'],right_IoU=b['otsu_iou'],left_material_success=a['material_image_success'],right_material_success=b['material_image_success'],left_actual_CPU=a['CPU'],right_actual_CPU=b['CPU'],left_actual_LU=a['LU'],right_actual_LU=b['LU'],left_actual_RHS=a['complex_RHS'],right_actual_RHS=b['complex_RHS']);records.append(record);paired.append(record)
    row=dict(comparison=label,left=left,right=right,cost_axis=axis,window_fraction=float(fraction),object_denominator=len(ids),paired_objects=len(paired),no_extrapolation=True)
    if paired:
     delta=np.array([r['NRMSE_difference_left_minus_right'] for r in paired]);ix=rng.integers(0,len(delta),size=(reps,len(delta)));boot=np.mean(delta[ix],axis=1)
     row.update(mean_NRMSE_difference=float(np.mean(delta)),mean_NRMSE_difference_CI95=np.quantile(boot,[.025,.975]).tolist(),median_NRMSE_difference=float(np.median(delta)),left_NRMSE_median=float(np.median([r['left_NRMSE'] for r in paired])),right_NRMSE_median=float(np.median([r['right_NRMSE'] for r in paired])),left_material_successes=sum(r['left_material_success'] for r in paired),right_material_successes=sum(r['right_material_success'] for r in paired),median_actual_target_budget=float(np.median([r['target_budget'] for r in paired])),median_left_unused_budget=float(np.median([r['left_unused_budget'] for r in paired])),median_right_unused_budget=float(np.median([r['right_unused_budget'] for r in paired])))
    summaries.append(row)
 return records,summaries

def plot_frontiers(summaries,folder):
 import matplotlib;matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 fig,axes=plt.subplots(len(COMPARISONS),3,figsize=(14,11),layout='constrained')
 for row,(left,right,label) in enumerate(COMPARISONS):
  for col,axis in enumerate(['CPU','LU','complex_RHS']):
   ax=axes[row,col];g=[r for r in summaries if r['comparison']==label and r['cost_axis']==axis and r.get('paired_objects',0)>0]
   if g:
    x=[r['window_fraction'] for r in g];y=[r['mean_NRMSE_difference'] for r in g];lo=[r['mean_NRMSE_difference_CI95'][0] for r in g];hi=[r['mean_NRMSE_difference_CI95'][1] for r in g]
    ax.plot(x,y,marker='.',linewidth=1.4);ax.fill_between(x,lo,hi,alpha=.2);ax.text(.02,.02,'Median actual budget %.2g -> %.2g'%(g[0]['median_actual_target_budget'],g[-1]['median_actual_target_budget']),transform=ax.transAxes,fontsize=8)
   ax.axhline(0,color='black',linewidth=.7);ax.set_title(left+' minus '+right+' | '+axis,fontsize=10);ax.set_xlabel("Relative budget in observed paired window");ax.set_ylabel('Mean paired material NRMSE difference');ax.grid(alpha=.2)
 fig.suptitle('Pairwise equal-cost frontiers;95% paired-object bootstrap CIs\nWindow=[max(first saved cost),min(last saved cost)]; no extrapolation; negative favors left',fontsize=12)
 fig.savefig(Path(folder)/'pairwise_actual_cost_frontiers.png',dpi=170);plt.close(fig)

def plot(rows,matched,folder):
 import matplotlib;matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 folder=Path(folder);folder.mkdir(parents=True,exist_ok=True)
 fig,axes=plt.subplots(1,3,figsize=(15,4.8),layout='constrained')
 for ax,field,title in zip(axes,['material_nrmse','otsu_iou','relative_residual'],['Material NRMSE','Prediction-only Otsu IoU','Original physical data residual']):
  values=[[r[field] for r in rows if r['arm']==arm and r.get(field) is not None] for arm in ARMS]
  ax.boxplot(values,tick_labels=ARMS,showfliers=False)
  for i,arm in enumerate(ARMS,1):
   group=[r for r in rows if r['arm']==arm and r.get(field) is not None];ax.scatter(i+np.linspace(-.16,.16,len(group)),[r[field] for r in group],c=[r['condition'] for r in group],cmap='viridis',vmin=0,vmax=4,s=18,alpha=.8)
  ax.set_title(title);ax.tick_params(axis='x',rotation=35);ax.grid(axis='y',alpha=.2)
  if field=='relative_residual':ax.set_yscale('log')
 axes[0].axhline(.15,color='grey',linestyle='--',linewidth=1);axes[1].axhline(.7,color='grey',linestyle='--',linewidth=1)
 fig.suptitle('30 independently generated objects; color=acquisition condition; six arms are paired, not independent samples')
 fig.savefig(folder/'continuous_material_distributions.png',dpi=170);plt.close(fig)
 fig,axes=plt.subplots(1,2,figsize=(10,4.5),layout='constrained');by={(r['case_id'],r['arm']):r for r in rows}
 for ax,(left,right,label) in zip(axes,COMPARISONS[:2]):
  pairs=[(by[(i,left)],by[(i,right)]) for i in sorted({r['case_id'] for r in rows}) if by[(i,left)].get('material_nrmse') is not None and by[(i,right)].get('material_nrmse') is not None]
  if not pairs:
   ax.set_title(label+' NO_SCORED_PAIRS');continue
  x=np.array([b['material_nrmse'] for a,b in pairs]);y=np.array([a['material_nrmse'] for a,b in pairs]);ax.scatter(x,y,c=[a['condition'] for a,b in pairs],cmap='viridis',vmin=0,vmax=4,s=28);lim=max(np.max(x),np.max(y))*1.03;ax.plot([0,lim],[0,lim],'k--',linewidth=.8);ax.set(xlabel=right+' NRMSE',ylabel=left+' NRMSE',title=label.replace('PRESPECIFIED_',''));ax.grid(alpha=.2)
 fig.savefig(folder/'paired_material_comparisons.png',dpi=170);plt.close(fig)
 fig,axes=plt.subplots(1,3,figsize=(15,4.7),layout='constrained')
 for ax,cost in zip(axes,['CPU','LU','complex_RHS']):
  subset=[r for r in matched if r.get('fraction_of_common_max')==1 and r['cost_axis']==cost and r.get('material_nrmse') is not None]
  values=[[r['material_nrmse'] for r in subset if r['arm']==arm] for arm in ARMS];ax.boxplot(values,tick_labels=ARMS,showfliers=False)
  for i,v in enumerate(values,1):ax.scatter(i+np.linspace(-.15,.15,len(v)),v,s=12,alpha=.65)
  ax.set_title('Per-object equal '+cost);ax.set_ylabel('Material NRMSE');ax.tick_params(axis='x',rotation=35);ax.grid(axis='y',alpha=.2)
 fig.suptitle('Common budget=min of six arms’ maximum saved costs; last paid checkpoint <=budget; every10-call/float32 approximation')
 fig.savefig(folder/'matched_actual_budget_material.png',dpi=170);plt.close(fig)

def secondary_aplus_dual(source,out):
 started,cpu=time.perf_counter(),time.process_time();source=Path(source);out=Path(out);rows=load(source/'SCORES.json')['rows'];checkpoints=load(source/'CHECKPOINT_SCORES.json')['rows'];rng=np.random.default_rng(2026110422)
 comparison=('a_plus','cie_dual','SECONDARY_USER_REQUESTED_APLUS_VS_MATURE_SOURCEWISE_DUAL')
 endpoint=bootstrap_pair(rows,*comparison,rng);frontiers,summary=pairwise_frontiers(rows,checkpoints,rng,comparisons=[comparison]);endstats=[]
 for axis in ['CPU','LU','complex_RHS']:
  g=[r for r in frontiers if r['cost_axis']==axis and r['window_fraction']==1 and r['status']=='PAIRED_LAST_SAVED_WITHIN_OBSERVED_RANGE'];paired_rows=[]
  for r in g:
   for side,arm in [('left','a_plus'),('right','cie_dual')]:paired_rows.append(dict(case_id=r['case_id'],arm=arm,material_nrmse=r[side+'_NRMSE'],otsu_iou=r[side+'_IoU'],material_image_success=r[side+'_material_success']))
  endstats.append(dict(cost_axis=axis,**bootstrap_pair(paired_rows,*comparison,rng)))
 by={(r['case_id'],r['arm']):r for r in rows};differences=[]
 for ident in sorted({r['case_id'] for r in rows}):
  a,b=by[(ident,'a_plus')],by[(ident,'cie_dual')];differences.append(dict(case_id=ident,condition=a['condition'],inverse_crime=a['inverse_crime'],Aplus_material_NRMSE=a['material_nrmse'],Dual_material_NRMSE=b['material_nrmse'],Aplus_success=a['material_image_success'],Dual_success=b['material_image_success'],Aplus_CPU=a['cpu_seconds'],Dual_CPU=b['cpu_seconds'],Aplus_LU=a['LU'],Dual_LU=b['LU'],Aplus_complex_RHS=a['actual_complex_RHS'],Dual_complex_RHS=b['actual_complex_RHS']))
 result=dict(classification='SECONDARY_USER_REQUESTED',comparison_scope='Aplus versus mature sourcewiseDual as complete algorithms, same30 paired objects; appended after primary scoring, original3 comparisons preserved',bootstrap_seed=2026110422,endpoint=endpoint,endpoint_object_differences=differences,pairwise_cost_frontiers=frontiers,pairwise_frontier_summaries=summary,common_max_cost_endpoint_bootstrap=endstats,noninferiority='NOT_TESTED: no posthoc margin introduced',new_physical_calls=0,wall_seconds=time.perf_counter()-started,cpu_seconds=time.process_time()-cpu)
 write(out/'SECONDARY_APLUS_VS_DUAL.json',result)
 lines=['# Secondary requested comparison: Aplus versus mature sourcewiseDual','', 'SECONDARY_USER_REQUESTED. This append does not alter the original three comparison ledgers. It compares complete algorithms on the same30 objects. No new threshold, noninferiority margin, optimizer run or physical call is introduced.','','|Endpoint metric|Mean Aplus-minus-Dual|95% paired object bootstrap CI|','|---|---:|---|']
 for f in ['material_nrmse','material_image_success','cpu_seconds','LU','actual_complex_RHS']:
  q=endpoint['metrics'][f];lines.append('|%s|%.6g|[%.6g, %.6g]|'%(f,q['mean_difference_left_minus_right'],*q['mean_difference_CI95']))
 lines+=['','|Paired common maximum axis|Aplus material successes|Dual material successes|Mean NRMSE delta|95% CI|','|---|---:|---:|---:|---|']
 for q in summary:
  if q['window_fraction']==1:lines.append('|%s|%d/30|%d/30|%.6g|[%.6g, %.6g]|'%(q['cost_axis'],q['left_material_successes'],q['right_material_successes'],q['mean_NRMSE_difference'],*q['mean_NRMSE_difference_CI95']))
 lines+=['','Each pair/object cost window is [max(first observed cost), min(last observed cost)] separately for actual process-prefix CPU, physical LU and complex RHS. Eleven fixed window fractions use the last saved paid original-Phi incumbent at or below budget. No stopped-method extrapolation. Auxiliary work is counted in CPU; equal LU is not equal CPU. Float32/every-ten-record plus sparse auxiliary-stage resolution limitations remain.','','Endpoint continuous means and paired medians answer different questions. Success-rate intervals include finite-object uncertainty; no formal noninferiority or population-wide superiority claim follows.']
 (out/'SECONDARY_APLUS_VS_DUAL.md').write_text('\n'.join(lines)+'\n');print(json.dumps(dict(status='SECONDARY_COMPLETE',out=str(out))))

def main():
 p=argparse.ArgumentParser();p.add_argument('--scores',required=True);p.add_argument('--output',required=True);p.add_argument('--secondary-aplus-dual',action='store_true');a=p.parse_args()
 if a.secondary_aplus_dual:secondary_aplus_dual(a.scores,a.output);return
 started,cpu=time.perf_counter(),time.process_time();source=Path(a.scores);out=Path(a.output)
 if out.exists():raise FileExistsError('Write-once analysis output exists')
 receipt=load(source/'SCORE_RECEIPT.json')
 if receipt['status']!='COMPLETE' or receipt['truth_access']!='POST_COMPLETE180_TASK_PUBLIC_SEAL_ONLY':raise RuntimeError('Need independent complete post-seal scores')
 rows=load(source/'SCORES.json')['rows'];checkpoints=load(source/'CHECKPOINT_SCORES.json')['rows'];rng=np.random.default_rng(2026110421)
 summaries=aggregate(rows,['arm']);strata=aggregate(rows,['condition','arm']);families=aggregate(rows,['family','arm']);pairs=[bootstrap_pair(rows,*c,rng) for c in COMPARISONS];matched=matched_cost(rows,checkpoints);frontier_rows,frontier_summary=pairwise_frontiers(rows,checkpoints,rng)
 matched_summaries=[];matched_pairs=[]
 for cost in ['CPU','LU','complex_RHS']:
  group=[r for r in matched if r.get('fraction_of_common_max')==1 and r['cost_axis']==cost]
  for arm in ARMS:
   g=[r for r in group if r.get('arm')==arm];matched_summaries.append(dict(cost_axis=cost,arm=arm,objects=len(g),scored=sum(r.get('material_nrmse') is not None for r in g),material_successes=sum(r.get('material_image_success') is True for r in g),material_nrmse=descriptive([r.get('material_nrmse') for r in g]),otsu_iou=descriptive([r.get('otsu_iou') for r in g]),unused_budget=descriptive([r.get('unused_budget') for r in g])))
  for left,right,label in COMPARISONS:matched_pairs.append(dict(cost_axis=cost,**bootstrap_pair(group,left,right,label,rng)))
 out.mkdir(parents=True);write(out/'ARM_SUMMARY.json',summaries);write(out/'CONDITION_SUMMARY.json',strata);write(out/'FAMILY_SUMMARY.json',families);write(out/'PAIRED_BOOTSTRAP.json',pairs);write(out/'MATCHED_COST_SCORES.json',matched);csvwrite(out/'MATCHED_COST_SCORES.csv',matched);write(out/'MATCHED_COST_SUMMARY.json',matched_summaries);write(out/'MATCHED_COST_PAIRED_BOOTSTRAP.json',matched_pairs);write(out/'PAIRWISE_COST_FRONTIERS.json',frontier_rows);csvwrite(out/'PAIRWISE_COST_FRONTIERS.csv',frontier_rows);write(out/'PAIRWISE_COST_FRONTIER_SUMMARY.json',frontier_summary)
 costs=[dict(arm=r['arm'],objects=r['objects'],**{k:r[k]['total'] for k in ['cpu_seconds','wall_seconds','LU','actual_complex_RHS','auxiliary_cpu_seconds']}) for r in summaries];write(out/'COST_TOTALS.json',dict(rows=costs,wall_scope='sum of per-task wall is work duration, not concurrent queue elapsed time; nested auxiliary CPU already in total, never add twice',failed_attempt_costs='Final receipt costs included; parent separate attempt ledger must add superseded/failed task attempts if any'))
 plot(rows,matched,ROOT/'FIGURES/primary');plot_frontiers(frontier_summary,ROOT/'FIGURES/primary')
 lines=['# Independent R7 primary analysis','', '30 independent object realizations, six paired arms; acquisition strata are object attributes. Final images use all576 material variables. No forward/adjoint calls are made by scoring or analysis.','','|Arm|Material successes /30|Same-grid material+noise /24|Median NRMSE|Median IoU|CPU total s|LU total|Complex RHS total|','|---|---:|---:|---:|---:|---:|---:|---:|']
 for r in summaries:
  display=lambda v,fmt: 'NA' if v is None else format(v,fmt)
  lines.append('|%s|%d/%d|%d/%d|%s|%s|%s|%s|%s|'%(r['arm'],r['material_successes'],r['objects'],r['same_grid_joint_successes'],r['same_grid_objects'],display(r['material_nrmse']['median'],'.4g'),display(r['otsu_iou']['median'],'.4g'),display(r['cpu_seconds']['total'],'.2f'),display(r['LU']['total'],'.0f'),display(r['actual_complex_RHS']['total'],'.0f')))
 lines+=['','## Paired object bootstrap95% confidence intervals','', 'Percentile bootstrap10000 paired object resamples, fixed seed2026110421. Mean deltas are left minus right. Two prespecified contrasts and one separately labeled secondary contrast; no multiplicity-adjusted inference. Failed tasks remain failures in success denominators; continuous metrics state their finite-pair count.','','|Comparison|Metric|Pairs|Mean delta|95% CI|','|---|---|---:|---:|---|']
 for p0 in pairs:
  for field in ['material_nrmse','material_image_success']:
   q=p0['metrics'].get(field)
   if q:lines.append('|%s vs %s|%s|%d|%.5g|[%.5g, %.5g]|'%(p0['left'],p0['right'],field,q['pairs'],q['mean_difference_left_minus_right'],*q['mean_difference_CI95']))
 lines+=['','## Scope and cost matching','', 'The24 same-grid objects permit the frozen independent CN99 noise discrepancy interpretation. For six36-to24 cross-grid objects, raw residual and noise-only pass are shown separately: physical model mismatch is unbudgeted, and sigma is not a complete data-error budget. Material recovery is always separate from data fit.','','Equal-cost grading uses each object’s common budget=min(maximum saved cost across six arms), separately for process-prefix CPU, physical LU, and total complex RHS. It takes the last saved original-Phi incumbent at or below budget; no interpolation or GT-based iterate selection. Snapshots every ten corrector observation events plus sparse auxiliary stage endpoints, and float32 material storage, make this a checkpoint-resolution approximation. Auxiliary CPU gaps can be large; unused-budget distributions are reported. CPU includes auxiliary work; equal LU/RHS do not imply equal CPU. Final serialization after last budget snapshot is outside that last checkpoint; BUDGET_AUDIT.json quantifies the tail.','','Pairwise A-versus-globalDual and A+-versus-A frontiers are primary cost comparisons; A-versus-sourcewiseDual is secondary. For each object and pair, the observed common window is [max(first saved cost),min(last saved cost)] on CPU/LU/RHS separately, sampled at11 fixed fractions. No stopped method is extrapolated. The all-six common-budget table is retained but can collapse when a method stops early at a poor local point, so it is not the sole cost comparison. Mean paired NRMSE differences and95% object-bootstrap intervals accompany exact per-object budgets and unused gaps. All30×6 images are shown in six fixed groups of five, without score-based selection. Continuous distributions, condition/family summaries and matched-budget material scores are retained. Acquisition conditions jointly vary contrast, noise, sources, aperture and generation grid; their descriptive differences do not isolate one factor. Bootstrap intervals across this fixed heterogeneous acquisition mixture do not establish population-wide universality.','','Historical R6 results may be quoted separately as already-public development/old-confirmation evidence; they are not added to this30-object denominator and are not rerun.']
 (out/'ANALYSIS.md').write_text('\n'.join(lines)+'\n')
 write(out/'ANALYSIS_RECEIPT.json',dict(status='COMPLETE',created_utc=datetime.now(timezone.utc).isoformat(),source=str(source),wall_seconds=time.perf_counter()-started,cpu_seconds=time.process_time()-cpu,new_physical_calls=0,GT_access='Only independent post-seal score arrays; no optimizer imports',bootstrap_replicates=10000,bootstrap_seed=2026110421,independent_objects=30,arms=6,conditions=5))
 print(json.dumps(dict(status='ANALYSED',out=str(out),objects=30)))
if __name__=='__main__':main()
