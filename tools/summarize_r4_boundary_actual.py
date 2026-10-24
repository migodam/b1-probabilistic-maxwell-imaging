#!/usr/bin/env python3
"""Read-only per-beta facts and shared-system allocation overlap afterclosure."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import sys,json,csv,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from p3r4.core import OUT,Book,write

ARMS=('full','richardson','global_gmres','global_gmres_accuracy')
def load(p):return json.loads(Path(p).read_text())
def csv_write(p,data):
 with p.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
def main():
 root=OUT/'results/II/interval_chain_boundary_v1';dest=OUT/'AUDIT/interval_boundary'
 if (dest/'ACTUAL_FACTS.json').exists():raise FileExistsError('preserve closedfacts')
 complete=load(root/'COORDINATOR_COMPLETE.json');comparison=load(root/'COMPARISON.json')
 if any(t.get('status')!='COMPLETE' for t in comparison['costs']):raise ValueError('this supplement expects actualall4complete grid')
 b=Book('II','interval_boundary_actual_fact_audit_v1',max_wall=600,analysis_only=True);status='COMPLETE';error=None
 try:
  cases=comparison['cases'];beta_rows=[];parent_rows=[];profiles=[]
  cfm_start=load(OUT/'I/jobs/geometry_cfm_DEV_infer_v1/START.json');cfm_receipt=OUT/'I/jobs/geometry_cfm_DEV_infer_v1/CLOSED.json'
  cfm_closed=load(cfm_receipt) if cfm_receipt.exists() else None;cfm_end=cfm_start['start_unix']+cfm_closed['controlled_wall_seconds'] if cfm_closed else time.time();overlap=[]
  for arm in ARMS:
   job='interval_chain_boundary_v1_'+arm;start=load(OUT/'II/jobs'/job/'START.json');receipt=load(OUT/'II/jobs'/job/'CLOSED.json')
   end=start['start_unix']+receipt['controlled_wall_seconds'];seconds=max(0.,min(end,cfm_end)-max(start['start_unix'],cfm_start['start_unix']))
   overlap.append(dict(arm=arm,arm_start_unix=start['start_unix'],arm_end_derived_unix=end,CFM_CPU_analysis_active_allocation_overlap_seconds=seconds,overlap_fraction_of_controlled_arm=seconds/receipt['controlled_wall_seconds'],
    CFM_receipt_closed_at_audit=cfm_closed is not None,scope='job allocation-window overlap from START+controlledwall,notCPU utilization measurement'))
   for beta in (.005,.01,.02):
    selected=[r for r in cases if r['arm']==arm and r['beta']==beta];n=sum(r['completed_steps'] for r in selected);accept=sum(r['accepts'] for r in selected)
    row=dict(arm=arm,beta=beta,independent_parents=4,seed_replicates=3,planned_cases=12,completed_cases=sum(r['coverage']=='COMPLETE' for r in selected),transitions=n,
     accepts=accept,acceptance=accept/n,Full_fallbacks=sum(r['Full_fallbacks'] for r in selected),Full_attempts=sum(r['Full_attempts'] for r in selected),
     total_RHS=sum(r['case_RHS'] for r in selected),cold_R0_RHS=sum(r['setup_RHS'] for r in selected),reduced_complex_LS_RHS=sum(r['reduced_LS_RHS'] for r in selected),
     mean_proposal_material_RMS=sum(r['mean_proposal_material_RMS']*r['completed_steps'] for r in selected)/n,
     mean_actual_material_move_RMS=sum(r['mean_actual_material_move_RMS']*r['completed_steps'] for r in selected)/n,
     mean_move_conditional_accept=sum(r['mean_actual_material_move_RMS']*r['completed_steps'] for r in selected)/accept if accept else None,
     max_proposal_material_RMS=max(r['max_proposal_material_RMS'] for r in selected),mean_terminal_RMS_from_anchor=float(np.mean([r['terminal_state_RMS_from_originalanchor'] for r in selected])),
     sum_case_wall_seconds=sum(r['case_wholewall_seconds'] for r in selected),recorded_setup_seconds=sum(r['recorded_setup_seconds'] for r in selected))
    beta_rows.append(row)
    for parent in ('R4DEV_00056','R4DEV_00057','R4DEV_00058','R4DEV_00059'):
     s=[r for r in selected if r['parent']==parent];steps=sum(r['completed_steps'] for r in s)
     parent_rows.append(dict(arm=arm,beta=beta,parent=parent,seeds=3,transitions=steps,accepts=sum(r['accepts'] for r in s),RHS=sum(r['case_RHS'] for r in s),
      Full_fallbacks=sum(r['Full_fallbacks'] for r in s),mean_proposal_material_RMS=sum(r['mean_proposal_material_RMS']*r['completed_steps'] for r in s)/steps,
      mean_actual_material_move_RMS=sum(r['mean_actual_material_move_RMS']*r['completed_steps'] for r in s)/steps))
   counts={};proofhist={};eta_invalid=0;rolecounts={};ls_shapes=set();small_numeric=0.;small_callback=0.;cachehits=0
   for case in [r for r in cases if r['arm']==arm]:
    folder=root/arm/case['case_id'];summary=load(folder/'SUMMARY.json')
    for key,value in summary['matrix_actions'].items():counts[key]=max(counts.get(key,0),value) if key.endswith('_peak_bytes') else counts.get(key,0)+value
    cachehits+=int(summary['setup'].get('shared_original_R0_cache_hit',False))
    for line in (folder/'ENDPOINT_PROOFS.jsonl').read_text().splitlines():
     p=json.loads(line);key=str(p['accuracy_iteration']);proofhist[key]=proofhist.get(key,0)+1;eta_invalid+=int(not p['exact'] and p.get('inverse_upper') is None)
     if p['event']=='paid_Full':rolecounts[p['role']]=rolecounts.get(p['role'],0)+1
    p=folder/'GLOBAL_GMRES_REDUCED_PROOFS.jsonl'
    if p.exists():
     for line in p.read_text().splitlines():
      r=json.loads(line)
      if r['event']=='global_GMRES_refined' and not r['solution_cache_hit'] and r['iteration']>0:
       ls_shapes.add((r['reduced_rows'],r['reduced_columns']));small_numeric+=r['reduced_numeric_solver_wall_seconds'];small_callback+=r['reduced_solver_wall_seconds']
   profiles.append(dict(arm=arm,recorded_refinement_proof_histogram=proofhist,invalid_inverse_bound_proof_events=eta_invalid,paid_Full_roles=rolecounts,original_R0_pool_reuse_case_hits=cachehits,
    global_reduced_LS_shapes=sorted(ls_shapes),reduced_LS_numeric_seconds=small_numeric,reduced_LS_ledger_callback_seconds=small_callback,matrix_actions=counts))
  # Same innovations + U, pathflags and storedattempts were checked in sealed
  # comparison; never change those checks or select a favorable subset here.
  facts=dict(status='COMPLETE_READONLY_PERBETA_FACTS',comparison_status=comparison['status'],numerical_misses=len(comparison['misses']),allpotentialproofs=comparison['all_potential_proofs_checked'],
   costs=comparison['costs'],perbeta=beta_rows,perparent_beta=parent_rows,numerical_profiles=profiles,CFM_CPU_allocation_overlap=overlap,
   coordinator_whole_wall_seconds=complete['coordinator_whole_wall_seconds'],coordinator_scope=complete['scope'],combined_complex_RHS=comparison['combined_complex_RHS'],
   observed_shared_system=True,quiet_system_wall_claim=False,independent_parents=4,settings_not_independent_parents=True,new_audit_RHS=0,new_audit_LU=0,new_audit_states=0,NN_updates=0,
   no_new_posterior_reference=True,float64_directed_interval_certificate=False,source_release=str(OUT/'configs/INTERVAL_CHAIN_BOUNDARY_RELEASE.json'))
  write(dest/'ACTUAL_FACTS.json',facts);csv_write(dest/'PER_BETA.csv',beta_rows);csv_write(dest/'PER_PARENT_BETA.csv',parent_rows);csv_write(dest/'CPU_ANALYSIS_OVERLAP.csv',overlap)
  report=['# 较大步长四臂：完整实测事实','',f"{comparison['status']}；36单元/臂全部完成，4个原父对象、9个设置/种子不当作36独立对象。实际复RHS合计 {comparison['combined_complex_RHS']} /300000，无重跑/扩额。全部共同随机、候选、flags/x路径、势proof和固定精度条件检查错误 {len(comparison['misses'])}；势proof {comparison['all_potential_proofs_checked']} 条。",'',
   '|方法|beta|接受/3072|Full回退|RHS|其中冷R0|小LS|平均proposal材料RMS|平均实际moveRMS（含拒绝0）|','|---|---:|---:|---:|---:|---:|---:|---:|---:|']
  for r in beta_rows:report.append(f"|{r['arm']}|{r['beta']}|{r['accepts']}|{r['Full_fallbacks']}|{r['total_RHS']}|{r['cold_R0_RHS']}|{r['reduced_complex_LS_RHS']}|{r['mean_proposal_material_RMS']:.6g}|{r['mean_actual_material_move_RMS']:.6g}|")
  report += ['', '所有approx臂各自只在首个beta/seed设置给4个原父对象冷付2304RHS，后32单元复用同一已付R0；单元basis不跨材料/设置重用。冷费在实际支付的beta=.005行分列，没有把摊销当成免费构造。所有Full回退失败与末段也纳入保守预约；本次全部完成。', '',
   '|臂|整外部进程秒|受控秒|物理调用秒|记录设置秒（含warmcopy）|总RHS|','|---|---:|---:|---:|---:|---:|']
  for r in comparison['costs']:report.append(f"|{r['arm']}|{r['whole_process_seconds']:.3f}|{r['controlled_seconds']:.3f}|{r['physical_operation_seconds']:.3f}|{r['recorded_setup_including_warm_copy_seconds']:.3f}|{r['RHS']}|")
  report += ['',f"协调器整段（含子进程、冷/热设置、所有状态/场/proof数组封存和内置postclose/报告）为 {complete['coordinator_whole_wall_seconds']:.3f} 秒。这里另做的只读beta汇总另有独立分析收据，不回写该已关闭时钟，也不重复计费子臂墙时。",'',
   '共享系统条件必须保留：释放时 geometry_cfm_DEV_infer_v1 为并行 CPU analysis，无并行授权物理。CPU_ANALYSIS_OVERLAP.csv 用START+受控作业区间给出每臂与CFM的分配墙时重叠；这不是CPU利用率测量，也不是隔离/静默机器基准。壁时比值仅为本次观测，未来quiet locked复制另计，不能用这里的重叠结果保证独立部署加速。','',
   'PER_PARENT_BETA.csv 保留每父对象/半径的三种子事实；原逐单元/逐步数组和拒绝全部保留。没有根据较大beta结果调小幅度、加修正次数或补物理；低接受率也原样保留，不转化成后验、混合或MCSE声明。普通同缓存GMRES的权利及小复LS费用均披露，不强求Richardson胜出。']
  (dest/'ACTUAL_REPORT_ZH.md').write_text('\n'.join(report)+'\n')
  with (root/'REPORT_ZH.md').open('a') as f:f.write('\n\n逐beta、逐父对象、共享CPU系统条件与外部协调器全时钟事实补充：AUDIT/interval_boundary/ACTUAL_REPORT_ZH.md；不是quiet机器壁时结论。\n')
  print(json.dumps(dict(status=facts['status'],combined_RHS=facts['combined_complex_RHS'],perbeta=beta_rows,CFM_overlap=overlap)),flush=True)
 except BaseException as exc:status='ERROR';error=repr(exc);raise
 finally:b.close(status,error)

if __name__=='__main__':main()
