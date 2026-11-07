"""Pure reindex of already scored commonB/final rows; no GT array or FW."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
import csv,json,argparse,time,html
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];D=ROOT/'delegated/twod/crossgrid8_paired';REFS=['cie_single','cie_dual','global_dual'];LABELS=['CIE','Dual CIE','Global Dual'];GRIDS=['matched24','cross48to24']
def rows(p):
 with Path(p).open(newline='') as f:return list(csv.DictReader(f))
def table(p,rr):
 with Path(p).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in rr for k in r)));w.writeheader();w.writerows(rr)
def yes(v):return v=='True'
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();wall=time.perf_counter();allrows=rows(D/'scored/ALL_SELECTION_SCORES.csv');final={(r['case_id'],r['grid_condition'],r['method']):r for r in allrows if r['view']=='FINAL_CAP_ELIGIBLE'};selected=[r for r in allrows if r['view']=='ACTUAL_CPU_PAIRED'];pairs={}
 for r in selected:pairs.setdefault((r['case_id'],r['pair_id']),{})[r['side']]=r
 if len(final)!=64 or len(pairs)!=128 or any(set(p)!= {'left','right'} for p in pairs.values()):raise ValueError('Original64/128pairs required')
 a.output.mkdir(exist_ok=False);wide=[];g19=[];transitions=[]
 for (case,pid),p in sorted(pairs.items()):
  l,r=p['left'],p['right'];out=dict(case_id=case,pair_id=pid,common_B=l['budget_B'],left_grid=l['grid_condition'],right_grid=r['grid_condition'],left_method=l['method'],right_method=r['method'],selection='Presealed min originalPhi within fixed commonB; no GT choice',exact_equal_execution_spend_claimed=False)
  if l['budget_B']!=r['budget_B']:raise ValueError('Different pair budgets')
  for side,member in p.items():
   for key in ['NRMSE','fixed_support_IoU','original_Phi','paid_CPU_seconds','paid_LU','paid_complex_RHS','material_success','actual_complete_child_CPU_seconds','task_status']:out[side+'_'+key]=member[key]
  wide.append(out)
  if pid.startswith('METHOD_PAIR_') and r['method']=='fullG19' and l['method'] in REFS:
   end=final[case,r['grid_condition'],'fullG19'];z=dict(out,grid_condition=r['grid_condition'],reference_method=l['method'],G19_final_NRMSE=end['NRMSE'],G19_final_IoU=end['fixed_support_IoU'],G19_final_original_Phi=end['original_Phi'],G19_final_material_success=end['material_success']);g19.append(z)
   if yes(r['material_success'])!=yes(end['material_success']):transitions.append(dict(z,transition='COMMON_B_PASS_FINAL_FAIL' if yes(r['material_success']) else 'COMMON_B_FAIL_FINAL_PASS',Phi_final_minus_B=float(end['original_Phi'])-float(r['original_Phi']),IoU_final_minus_B=float(end['fixed_support_IoU'])-float(r['fixed_support_IoU'])))
 table(a.output/'COMMON_B_ALL128_PAIRED_METRICS.csv',wide);table(a.output/'COMMON_B_G19_COMPARATOR_PER_OBJECT.csv',g19);table(a.output/'G19_SUCCESS_VIEW_TRANSITIONS.csv',transitions)
 summary=[];cost=[]
 for g in GRIDS:
  for ref in REFS:
   rr=[r for r in g19 if r['grid_condition']==g and r['reference_method']==ref];summary.append(dict(grid_condition=g,reference_method=ref,independent_materials=8,pairs=len(rr),common_B_rule='min(actual valid-child endCPU of reference, actual valid-child endCPU of G19), pair-specific',common_B_min=min(float(r['common_B']) for r in rr),common_B_max=max(float(r['common_B']) for r in rr),reference_common_B_success=sum(yes(r['left_material_success']) for r in rr),G19_common_B_success=sum(yes(r['right_material_success']) for r in rr),G19_full_endpoint_success=sum(yes(r['G19_final_material_success']) for r in rr),primary_endpoint_not_replaced=True))
  for method in REFS+['fullG19']:
   rr=[r for r in final.values() if r['grid_condition']==g and r['method']==method];total=sum(float(r['actual_complete_child_CPU_seconds']) for r in rr);cost.append(dict(grid_condition=g,method=method,independent_materials=8,total_actual_child_CPU_seconds=total,mean_actual_child_CPU_seconds=total/8,material_success=sum(yes(r['material_success']) for r in rr),CN99_pass=sum(yes(r['CN99_pass']) for r in rr),CN99_denominator=8))
 table(a.output/'COMMON_B_G19_COMPARATOR_SUMMARY.csv',summary);table(a.output/'COST_CN99_BY_GRID_METHOD.csv',cost);artifacts=[]
 for metric,label in [('NRMSE','NRMSE'),('fixed_support_IoU','Support IoU')]:
  fig,axs=plt.subplots(2,3,figsize=(18,9),squeeze=False)
  for i,g in enumerate(GRIDS):
   for j,(ref,name) in enumerate(zip(REFS,LABELS)):
    rr=sorted([r for r in g19 if r['grid_condition']==g and r['reference_method']==ref],key=lambda r:r['case_id']);ax=axs[i,j];x=np.arange(8);ax.plot(x,[float(r['left_'+metric]) for r in rr],'o-',label=name,c='#2166ac');ax.plot(x,[float(r['right_'+metric]) for r in rr],'s-',label='Full G19',c='#d95f02');ax.axhline(.15 if metric=='NRMSE' else .7,c='.5',ls='--',lw=.8);ax.set_xticks(x,[r['case_id'][-4:]+'\nB='+f"{float(r['common_B']):.1f}" for r in rr],fontsize=8);ax.set_ylabel(label);ax.set_title(g+' / '+name+' vs G19');ax.grid(alpha=.2);ax.legend(fontsize=8)
  fig.suptitle('All three presealed commonB rules; B labels rounded for display, exact B in CSV; main endpoint remains6/8');fig.tight_layout(rect=[0,0,1,.95],h_pad=3)
  for ext in ('png','svg'):
   p=a.output/('COMMON_B_'+metric+'.'+ext);fig.savefig(p,dpi=180,bbox_inches='tight');artifacts.append(p.name)
  plt.close(fig)
 lines=['# Complete commonB comparisons and endpoint transition audit','','Pure reindex of already authorized ALL_SELECTION_SCORES.csv. No GT array, new scoring, FW, optimizer, iteration or budget selection. All three reference-vs-G19 rules are shown under each grid condition. B is pair-specific min(actual reference endCPU, actual G19 endCPU); both original full children remain charged. Main final endpoint stays6/8 G19, and no7/8 commonB view replaces it.','','|Grid|Reference|Reference atB|G19 atB|G19 final|B min..max CPU seconds|','|---|---|---:|---:|---:|---|']
 for r in summary:lines.append('|'+r['grid_condition']+'|'+r['reference_method']+'|'+str(r['reference_common_B_success'])+'/8|'+str(r['G19_common_B_success'])+'/8|'+str(r['G19_full_endpoint_success'])+'/8|'+str(r['common_B_min'])+' .. '+str(r['common_B_max'])+'|')
 lines+=['','Exact B for every independent material (not a common fixed budget across pairs):','','|Case|matched CIE|matched Dual|matched Global|cross CIE|cross Dual|cross Global|','|---|---:|---:|---:|---:|---:|---:|']
 for case in sorted({r['case_id'] for r in g19}):lines.append('|'+case+'|'+'|'.join(next(r['common_B'] for r in g19 if r['case_id']==case and r['grid_condition']==g and r['reference_method']==ref) for g in GRIDS for ref in REFS)+'|')
 lines+=['','Transitions below describe lower-originalPhi endpoint selection, never a GT-selected budget or iterate:','','|Case/grid/reference|B|G19 NRMSE atB→final|G19 IoU atB→final|G19 originalPhi atB→final|Transition|','|---|---:|---|---|---|---|']
 for r in transitions:lines.append('|'+r['case_id']+'/'+r['grid_condition']+'/'+r['reference_method']+'|'+r['common_B']+'|'+r['right_NRMSE']+' → '+r['G19_final_NRMSE']+'|'+r['right_fixed_support_IoU']+' → '+r['G19_final_IoU']+'|'+r['right_original_Phi']+' → '+r['G19_final_original_Phi']+'|'+r['transition']+'|')
 lines+=['','All128 pair metrics and all48 reference-vs-G19 per-object material/cost/Phi records are retained in the adjacent CSVs. COST_CN99_BY_GRID_METHOD.csv supplies total/mean actual CPU and CN99 counts for all8 grid-method groups.','','![CommonB NRMSE](COMMON_B_NRMSE.png)','','![CommonB IoU](COMMON_B_fixed_support_IoU.png)','']
 (a.output/'REPORT.md').write_text('\n'.join(lines));(a.output/'RECEIPT.json').write_text(json.dumps(dict(status='COMPLETED_PURE_REINDEX_ALL128_PAIRS',source=str(D/'scored/ALL_SELECTION_SCORES.csv'),all_pairs=128,G19_comparison_pairs=48,independent_materials=8,all_three_reference_rules=True,summary=summary,transitions=transitions,artifacts=artifacts,GT_access='NO_ARRAY_READ_CSV_REINDEX_ONLY',new_FW=0,new_score=0,CPU_seconds=time.process_time(),wall_seconds=time.perf_counter()-wall),indent=2)+'\n');report=D/'plots/REPORT.md';report.write_text(report.read_text()+'\n## All three commonB views, separately from primary endpoint\n\nSee [complete reindex and exact pair-specific B](../common_b_reindex/REPORT.md). MainG19 endpoint remains6/8 under both grids. All CIE/Dual/Global comparison rules and all128 paired records are retained; a7/8 intermediate commonB view is never substituted for this endpoint.\n\n'+ '\n'.join(lines[4:]) .replace('(COMMON_B_NRMSE.png)','(../common_b_reindex/COMMON_B_NRMSE.png)').replace('(COMMON_B_fixed_support_IoU.png)','(../common_b_reindex/COMMON_B_fixed_support_IoU.png)'))
 htmlfile=D/'plots/REPORT.html';htmlfile.write_text(htmlfile.read_text()+'<h2>All three commonB views</h2><p>Main G19 endpoint remains6/8. See <a href="../common_b_reindex/REPORT.md">complete exact-B tables and endpoint-transition audit</a>; no7/8 intermediate view replaces the primary endpoint.</p><img src="../common_b_reindex/COMMON_B_NRMSE.png"><img src="../common_b_reindex/COMMON_B_fixed_support_IoU.png">');print(json.dumps(dict(summary=summary,transitions=transitions)))
if __name__=='__main__':main()
