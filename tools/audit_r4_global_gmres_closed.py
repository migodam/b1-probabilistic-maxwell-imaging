#!/usr/bin/env python3
"""Postclose cost/randomness inspection only. Never restarts a physical arm."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import sys,json,csv
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from p3r4.core import OUT,Book,write

ARMS=('full','richardson','global_gmres','global_gmres_accuracy')
PARENTS=tuple(f'R4DEV_{i:05d}' for i in range(48,56))
def main():
 root=OUT/'results/II/interval_chain_global_gmres_v1';dest=OUT/'AUDIT/global_gmres';comparison=json.loads((root/'COMPARISON.json').read_text())
 if (dest/'EXECUTED_AUDIT.json').exists():raise FileExistsError('preserve executed audit')
 gates={arm:json.loads((root/arm/'ARM_COMPLETE.json').read_text()) for arm in ARMS}
 if any(g['status']!='COMPLETE' for g in gates.values()):raise ValueError('all four complete first')
 b=Book('II','global_gmres_postclose_audit_v1',max_wall=600,analysis_only=True);status='COMPLETE';error=None
 try:
  rows=[];random_misses=[];array_misses=[];small=[];numerical=[]
  for arm in ARMS:
   gate=gates[arm];receipt=json.loads(Path(gate['receipt']).read_text());result=json.loads((root/arm/'RESULTS.json').read_text())
   if receipt['status']!='COMPLETE':raise ValueError('closed receipt required')
   counts={};setup=0.;fallback=0;accepts=0;proofevents=0;invalid=0;accuracy_iterations={};actual_ls=[]
   for parent in PARENTS:
    summary=json.loads((root/arm/parent/'SUMMARY.json').read_text());fallback+=summary['fallbacks'];accepts+=summary['accepts'];setup+=summary['setup'].get('wall_including_R0_boundcopy_and_persistence',0.)
    for key,v in summary['matrix_actions'].items():
     counts[key]=max(counts.get(key,0),v) if key.endswith('_peak_bytes') else counts.get(key,0)+v
    refs=[json.loads(line) for line in (root/arm/parent/'CHAIN_REFERENCES.jsonl').read_text().splitlines()]
    if arm!='full':
     gold=[json.loads(line) for line in (root/'full'/parent/'CHAIN_REFERENCES.jsonl').read_text().splitlines()]
     if len(refs)!=1024 or len(gold)!=1024:raise ValueError('1024 steps required')
     for k,(r,g) in enumerate(zip(refs,gold)):
      if r['logu']!=g['logu'] or r['step']!=g['step']:random_misses.append(dict(arm=arm,parent=parent,step=k))
      # Actual stored rejected candidates as well as accepted paths are checked.
      with np.load(root/arm/parent/'attempts'/f'c{k:04d}.npz',allow_pickle=False) as z:a={key:z[key] for key in ('x','q','chi','logu')}
      with np.load(root/'full'/parent/'attempts'/f'c{k:04d}.npz',allow_pickle=False) as z:
       if any(not np.array_equal(a[key],z[key]) for key in a):array_misses.append(dict(arm=arm,parent=parent,step=k))
    for line in (root/arm/parent/'ENDPOINT_PROOFS.jsonl').read_text().splitlines():
     p=json.loads(line);proofevents+=1;invalid+=not p['exact'] and (p.get('inverse_upper') is None)
     key=str(p['accuracy_iteration']);accuracy_iterations[key]=accuracy_iterations.get(key,0)+1
    path=root/arm/parent/'GLOBAL_GMRES_REDUCED_PROOFS.jsonl'
    if path.exists():
     for line in path.read_text().splitlines():
      p=json.loads(line)
      # Each refinement proof is emitted once, solution-cache replays charge0.
      if p['event']=='global_GMRES_refined' and not p['solution_cache_hit'] and p['iteration']>0:actual_ls.append(p)
   primary=gate['costs']['states']*8+(4608 if arm!='full' else 0);reduced=gate['costs']['RHS']-primary
   if reduced<0:raise AssertionError('negative reduced solve debit')
   if arm.startswith('global_gmres') and reduced!=counts.get('global_gmres_reduced_LS_complex_RHS',0):raise AssertionError('small LS debit/counter mismatch')
   if not arm.startswith('global_gmres') and reduced:raise AssertionError('unexpected reduced solve cost')
   rows.append(dict(arm=arm,total_complex_RHS=gate['costs']['RHS'],primary_full_grid_RHS=primary,reduced_complex_LS_RHS=reduced,Full_states=gate['costs']['states'],Full_LU=gate['costs']['LU'],
    Full_fallbacks=fallback,accepts=accepts,controlled_seconds=receipt['controlled_wall_seconds'],CPU_seconds=receipt['profile']['CPU_seconds'],RSS_peak_bytes=receipt['profile']['RSS_peak_bytes'],
    controlled_physics_seconds=None,whole_process_seconds=comparison['costs'][arm]['whole_process_seconds'],cold_R0_bound_setup_persist_seconds=setup,
    reduced_LS_numeric_seconds=counts.get('global_gmres_reduced_LS_numeric_nanoseconds',0)/1e9,reduced_LS_callback_seconds=counts.get('global_gmres_small_solver_nanoseconds',0)/1e9,
    proof_events=proofevents,invalid_inverse_bound_events=invalid))
   numerical.append(dict(arm=arm,matrix_actions=counts,proof_refinement_histogram=accuracy_iterations))
   if actual_ls:small.append(dict(arm=arm,actual_small_LS_records=len(actual_ls),charged_reduced_complex_RHS=reduced,
    shapes=sorted(set((p['reduced_rows'],p['reduced_columns']) for p in actual_ls)),happy_breakdowns=sum(p['happy_breakdown'] for p in actual_ls),
    maximum_preconditioned_initial_residual=max(p['preconditioned_initial_residual_F'] for p in actual_ls)))
  # Read complete ledger records only after arms close; no mutation/recharge.
  selected={g['receipt']:arm for arm,g in gates.items()};jobs={'interval_chain_global_gmres_v1_'+a:a for a in ARMS};physical={a:0. for a in ARMS};rhs={a:0 for a in ARMS};reduced_ledger={a:0 for a in ARMS}
  with (OUT/'budget.lock').open('a') as lock:
   import fcntl
   fcntl.flock(lock,fcntl.LOCK_SH)
   with (OUT/'RESOURCE_LEDGER.jsonl').open('rb') as f:f.seek(0,2);end=f.tell()
   fcntl.flock(lock,fcntl.LOCK_UN)
  with (OUT/'RESOURCE_LEDGER.jsonl').open('rb') as f:
   while f.tell()<end:
    line=f.readline()
    if not line.endswith(b'\n'):break
    # Skip JSON decoding unrelated campaign jobs.
    if b'interval_chain_global_gmres_v1_' not in line:continue
    entry=json.loads(line);arm=jobs.get(entry['job'])
    if arm is None:continue
    physical[arm]+=entry.get('physics_wall_seconds',0.);rhs[arm]+=entry.get('RHS',0)
    if entry.get('role')=='samecache_global_GMRES_reduced_LS':reduced_ledger[arm]+=entry['RHS']
  for r in rows:
   a=r['arm'];r['controlled_physics_seconds']=physical[a]
   if rhs[a]!=r['total_complex_RHS'] or reduced_ledger[a]!=r['reduced_complex_LS_RHS']:raise AssertionError('closed receipt vs ledger cost mismatch')
  result=dict(status='PASS_POSTCLOSE_PATH_AND_COST_RECONCILIATION' if not random_misses and not array_misses else 'NO_GO_RANDOM_OR_PROPOSAL_MISMATCH',
   comparison_status=comparison['status'],comparison_misses=comparison['total_misses'],uniform_mismatches=random_misses,stored_candidate_x_q_chi_logu_mismatches=array_misses,
   costs=rows,numerical_profiles=numerical,actual_reduced_solve_profiles=small,combined_complex_RHS=sum(r['total_complex_RHS'] for r in rows),
   bounds='ordinary float64 enclosures, not directed certificates',posterior_qualification='NOT_ESTABLISHED',new_audit_RHS=0,new_audit_states=0,new_audit_LU=0,
   ledger_complete_prefix_bytes=end,source_release=str(OUT/'configs/INTERVAL_CHAIN_GLOBAL_GMRES_RELEASE.json'))
  write(dest/'EXECUTED_AUDIT.json',result)
  with (dest/'EXECUTED_COSTS.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
  lines=['# 同权限 global GMRES 四臂：实际结果附录','',f"闭合比较：{comparison['status']}，全部判决、路径和中间势包络错误 {comparison['total_misses']}。额外逐项检查已存的候选 x/q/chi/logU（包含拒绝候选）与 logU，差异 {len(array_misses)} / {len(random_misses)}。",'',
   '|方法|总复 RHS|全网格 RHS|小 LS RHS|Full 状态/LU|Full 回退|外部整进程秒|受控秒|物理秒|','|---|---:|---:|---:|---:|---:|---:|---:|---:|']
  for r in rows:lines.append(f"|{r['arm']}|{r['total_complex_RHS']}|{r['primary_full_grid_RHS']}|{r['reduced_complex_LS_RHS']}|{r['Full_states']}|{r['Full_fallbacks']}|{r['whole_process_seconds']:.3f}|{r['controlled_seconds']:.3f}|{r['controlled_physics_seconds']:.3f}|")
  lines += ['',f"组合费用 {result['combined_complex_RHS']} / 300,000 RHS，所有方法按预分配子预算运行一次，没有重跑或增加迭代。每种近似方法各自冷付 8×576=4,608 逆列 RHS，同等享有 R0/G/S 八源矩阵运算；R0 矩阵乘法不是新 LU 求解。GMRES 是把 576×8 电流视作单个复向量的普通 left-preconditioned global GMRES，非八套独立或 block GMRES。小型复 LS 逐次先保守计 1 RHS，再作实际 ≤5×4 一右端解；小 LS 与全网格项均列出。",'',
   '接受近似点没有变成新的精确锚点，原始 R0 保持不变；两个端点的当前精度/basis 可保留，跨材料 basis 不复用。所有尝试、拒绝、Full 回退、最终电流/场和版本/中间 proof 标量保留。共同 stored/canonical endpoint-once 存储权限适用于全部方法；外部整进程时间包含导入、冷 R0、写入、proof 和审计前封存。Book 受控时间与物理调用时间分列，小 LS/矩阵计算不属于 Full 物理时间。', '',
   '普通 Krylov 同样具备已付缓存权限，因此结果不能自动归因于 A17/current-space 理论，也没有要求 Richardson 胜出。decision 方法与固定精度方法比较的是不同停止准则的实际代价，不能将固定精度更贵泛化为所有科学任务均无需精度。该固定 pCN 核的数值一致性不证明混合、独立样本、MCSE、后验参考或 float64 定向区间证书。', '',
   '原协议 PROTOCOL_ZH.md 和全部封账结果未改写。机器证据 EXECUTED_AUDIT.json、EXECUTED_COSTS.csv；正式逐项比较 results/II/interval_chain_global_gmres_v1/COMPARISON.json。float32 材料框只读补充见 AUDIT/raw_risk_actions_DEV_float32/REPORT_ZH.md。最终科学获益和普通 Krylov 是否解释观察收益由父线程判断。']
  (dest/'EXECUTED_REPORT_ZH.md').write_text('\n'.join(lines)+'\n')
  print(json.dumps(dict(status=result['status'],combined_RHS=result['combined_complex_RHS'],costs=rows)),flush=True)
 except BaseException as exc:status='ERROR';error=repr(exc);raise
 finally:b.close(status,error)

if __name__=='__main__':main()
