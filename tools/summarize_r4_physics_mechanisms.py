#!/usr/bin/env python3
"""Factual closed-result extraction plus inert locked wrapper preparation."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import sys,json,csv,subprocess,ast
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from p3r4.core import OUT,Book,write
from p3r4.operator_surrogate import BASELINES

def load(p):return json.loads(Path(p).read_text())
def closed(stage,job):
 p=OUT/stage/'jobs'/job/'CLOSED.json';d=load(p)
 if d['status']!='COMPLETE':raise ValueError('closed COMPLETE only: '+str(p))
 return p,d
def csv_write(path,rows):
 with path.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def main():
 dest=OUT/'AUDIT/physics_mechanism_tables';summary=OUT/'AUDIT/PHYSICS_MECHANISM_SUMMARY.md'
 if dest.exists() or summary.exists():raise FileExistsError('new frozen factual summary only')
 b=Book('II','physics_mechanism_summary_and_locked_prepare_v1',max_wall=600,analysis_only=True);status='COMPLETE';error=None
 try:
  dest.mkdir();dino=[];operator=[];geometry=[];owner=[]
  # Preparation is explicit and inert; neither subcommand can call Book/Target.
  p=subprocess.run([sys.executable,'tools/run_r4_interval_chain_locked.py','prepare'],cwd=ROOT,capture_output=True,text=True)
  (dest/'LOCKED_PREPARE.log').write_text(p.stdout+p.stderr)
  if p.returncode:raise RuntimeError('inert locked prepare failed')
  p=subprocess.run([sys.executable,'tools/run_r4_interval_chain_locked.py','preflight'],cwd=ROOT,capture_output=True,text=True)
  (dest/'LOCKED_PREFLIGHT.json').write_text(p.stdout)
  if p.returncode:raise RuntimeError('locked prereq inspection failed')
  for source in ('tools/run_r4_interval_chain_locked.py','tools/run_r4_interval_chain_fast.py'):
   ast.parse((ROOT/source).read_text())
  test=subprocess.run([sys.executable,'-m','unittest','discover','-s','tests/p3r4','-p','test_locked_fourarm_preparation.py','-v'],cwd=ROOT,env={**os.environ,'PYTHONPATH':str(ROOT/'src')},capture_output=True,text=True)
  (dest/'LOCKED_WRAPPER_TESTS.log').write_text(test.stdout+test.stderr)
  if test.returncode:raise RuntimeError('zero-physics wrapper checks failed')
  for mode in ('value','jvp'):
   for seed in (2026102501,2026102502,2026102503):
    job=f'forward_{mode}_66_s{seed}';folder=OUT/'results/II/forward_jvp_evaluation'/(job+'_evaluation')
    e=load(folder/'RESULTS.json');t=load(OUT/'results/II/forward_jvp_training'/job/'RESULTS.json');rp=load(OUT/'results/II/forward_jvp_training'/job/'RUN_PROTOCOL.json')
    receipt,c=closed('II',job);er,ec=closed('II',job+'_evaluation');m=e['metrics'];white=[];phi=[]
    for parent in m['parents']:
     with np.load(folder/(parent['parent']+'.npz'),allow_pickle=False) as z:white.extend(z['white_L2_error']);phi.extend(abs(z['Phi_error']))
    row=dict(mode=mode,seed=seed,parents=len(m['parents']),pairs=len(white),optimizer_updates=rp['updates'],selected_update=t['selected_update'],
     parent_balanced_white_MSE=m['parent_balanced_white_MSE'],max_white_L2=m['max_white_L2_error'],white_L2_q95_descriptive=float(np.quantile(white,.95)),
     max_absolute_Phi_error=m['max_absolute_Phi_error'],absolute_Phi_error_q95_descriptive=float(np.quantile(phi,.95)),
     original_validation_all_point_gate=e['validation_gate'],operator_test_all_point_gate=m['empirical_validation_gate'],
     warm_query_seconds=m['single_query_compute_seconds_mean'],cold_checkpoint_and_torch_seconds=e['cold_checkpoint_and_torch_load_seconds'],cold_model_seconds=e['cold_model_build_and_state_load_seconds'],
     training_controlled_seconds=c['controlled_wall_seconds'],training_CPU_seconds=c['profile']['CPU_seconds'],evaluation_controlled_seconds=ec['controlled_wall_seconds'],
     shared_teacher_owner_RHS_not_additive=459712,cold_anchor_FullJ_RHS=40,actual_new_physics_RHS=e['actual_new_physics_RHS'],
     required_Full_fallback_RHS_per_candidate=m['Full_fallback_complex_RHS_per_candidate'],counterfactual_all_fallback_RHS_not_paid=m['counterfactual_all_Full_fallback_RHS'],
     source_result=str(folder/'RESULTS.json'),training_receipt=str(receipt),evaluation_receipt=str(er))
    dino.append(row)
  for seed in (2026102501,2026102502,2026102503):
   job=f'operator_ranker66_s{seed}';folder=OUT/'results/II/operator_evaluation'/(job+'_evaluation');e=load(folder/'RESULTS.json');t=load(OUT/'results/II/operator_training'/job/'RESULTS.json')
   receipt,c=closed('II',job);er,ec=closed('II',job+'_evaluation')
   history=[json.loads(s) for s in (OUT/'results/II/operator_training'/job/'VALIDATION.jsonl').read_text().splitlines()]
   row=dict(seed=seed,parents=len(e['parents']),pairs=sum(p['pairs'] for p in e['parents']),optimizer_updates=max(h['update'] for h in history),selected_update=t['selected_update'],
    selected_material_action_error_parent_mean=float(np.mean([p['selected_observed_error_mean'] for p in e['parents']])),
    training_controlled_seconds=c['controlled_wall_seconds'],training_CPU_seconds=c['profile']['CPU_seconds'],evaluation_controlled_seconds=ec['controlled_wall_seconds'],
    shared_teacher_owner_RHS_not_additive=459712,cold_anchor_FullJ_RHS=40,candidate_Full_feature_RHS=8,actual_new_physics_RHS=e['new_physics_RHS'],
    inference_scope='candidate Full response/current-field features required; action error ranker, not forward-only saving',source_result=str(folder/'RESULTS.json'),training_receipt=str(receipt),evaluation_receipt=str(er))
   for i,name in enumerate(BASELINES):row[name+'_material_action_error_parent_mean']=float(np.mean([p['fixed_baseline_error_mean'][i] for p in e['parents']]))
   operator.append(row)
  gpath=OUT/'results/III/frozen_geometry_probability_v1/RESULTS.json';g=load(gpath);greceipt,gc=closed('III','frozen_geometry_probability_v1')
  first=load(OUT/'results/III/dynamic_geometry_first/RESULTS.json');heun=load(OUT/'results/III/dynamic_geometry_heun_r1/RESULTS.json')
  with np.load(OUT/'results/III/frozen_geometry_probability_v1/MATERIAL_FUNCTIONS.npz',allow_pickle=False) as z:function_names=z['names'].astype(str).tolist()
  for i in range(8):
   parent=f'R4DEV_{i:05d}';context=load(OUT/'results/contexts/dev'/(parent+'.json'))
   ownfirst=sum(p['actual_new_costs']['RHS'] for p in first['configurations'] if p['parent']==parent)
   ownheun=sum(p['actual_new_cost']['RHS'] for p in heun['configurations'] if p['parent']==parent)
   owner.append(dict(parent=parent,actual_shared_online_context_RHS=context['new_RHS'],actual_shared_online_context_states=context['new_states'],
    actual_shared_online_context_seconds=context['wall_seconds'],anchor_Phi=context['Phi'],
    paid_first_geometry_all_configurations_RHS=ownfirst,paid_Heun_all_configurations_RHS=ownheun,selected_rank4_66_centralFD_shared_RHS=480,
    ownership_scope='context and geometry owner once perparent, shared by all four kernels; first+Heun are full research acquisition, not marginal subset bill'))
  for p in g['parents']:
   diag=p['material_diagnostics'];ess=diag['ESS_initial_positive_pair'];finite=[v for v in ess if v is not None and np.isfinite(v)]
   row=dict(parent=p['parent'],kernel=p['kernel'],steps=p['steps'],accepts=p['accepts'],RHS=p['costs']['RHS'],states=p['costs']['states'],LU=p['costs']['LU'],wall_seconds=p['wall_seconds'],
    mean_material_move_RMS=p['mean_material_move_RMS'],expected_material_linear_RMS=p['kernel_normalization']['material_linear_expected_RMS'],cap_binds=p['kernel_normalization']['cap_binds'],
    finite_function_ESS_min=min(finite),finite_function_ESS_median=float(np.median(finite)),constant_or_undefined_functions=sum(v is None for v in ess),
    qualification=diag['qualification'],common_reference=diag['common_reference'],source_result=str(gpath))
   for name,v in zip(function_names,ess):row['ESS_'+name]=v
   geometry.append(row)
  csv_write(dest/'DINO_SIX_FITS.csv',dino);csv_write(dest/'OPERATOR_THREE_FITS.csv',operator);csv_write(dest/'FROZEN_GEOMETRY_PER_PARENT_KERNEL.csv',geometry);csv_write(dest/'SHARED_PARENT_ACQUISITION.csv',owner)
  facts=dict(status='CLOSED_RESULT_FACT_EXTRACTION',DINO=dino,operator_ranker=operator,frozen_geometry=geometry,shared_parent_acquisition=owner,
   shared_geometry_teacher_owner_RHS=459712,shared_teacher_J_state_RHS=414040,teacher_unique_pairs=10063,
   teacher_fold='parents0..15train,16..23validation,24..31operator_test;EXPOSED development',
   DINO_original_gate='ALL evaluated points maxwhiteL2<=.05 AND maxabsPhierror<=.05;95thpercentiles descriptive only, no new95%gate',
   frozen_geometry_total_new_RHS=g['costs']['RHS'],frozen_geometry_controlled_seconds=gc['controlled_wall_seconds'],frozen_geometry_receipt=str(greceipt),
   frozen_geometry_common_reference='NONE',frozen_geometry_qualification='UNQUALIFIED finite-run ESS; no pooling acrossparents, no MCSE certificate',
   locked_preparation='AUDIT/locked_samecache/FROZEN_PROTOCOL.json',locked_new_physics=0,summary_new_physics_RHS=0,summary_new_neural_updates=0,no_GT=True,
   source_snapshot_scope='future default-preserving fast ledgerstage/context/cache hooks; closed as-run snapshots preserved; pairedphysicsVI sources not modified')
  write(OUT/'AUDIT/PHYSICS_MECHANISM_FACTS.json',facts)
  lines=['# 物理机制控制：封账事实表','',
   '本页只汇总已完成的 6 次局部前向/JVP 拟合、3 次动作误差排序器拟合和前 8 个 DEV 的 4 种冻结几何核。此次提取不作物理查询、不拟合、不读取真值，也不作最终 GO 判断。细表和来源/收据路径见 physics_mechanism_tables/ 与 PHYSICS_MECHANISM_FACTS.json。','',
   '三个模块共享先前几何教师拥有成本 459,712 RHS（Full-J 状态部分 414,040 RHS；10,063 独特配对），不能在每个种子重复计入，也不能因复用将冷标签说成免费。原在线锚点取得另列 SHARED_PARENT_ACQUISITION.csv；若部署时未缓存，锚点 Full+接收器 J 的冷构造为 8+32 RHS。','',
   '## 局部前向/JVP：六次真实拟合','',
   '|监督|种子尾号|更新/所选|最大白化误差|白化误差95分位（描述）|最大绝对 Phi 误差|原验证/算子测试门限|','|---|---:|---:|---:|---:|---:|---|']
  for r in dino:lines.append(f"|{r['mode']}|{r['seed']%1000}|{r['optimizer_updates']}/{r['selected_update']}|{r['max_white_L2']:.4g}|{r['white_L2_q95_descriptive']:.4g}|{r['max_absolute_Phi_error']:.4g}|{r['original_validation_all_point_gate']}/{r['operator_test_all_point_gate']}|")
  lines += ['', '实际保存门限是所有点的最大白化 L2≤.05 且最大绝对 Phi 误差≤.05；这里不增造“95%通过”规则。95 分位只是误差分布描述。六次都未通过实际门限，因此短推理时间不构成前向加速资格；每个候选需要 8 RHS Full 回退的要求保留，反事实回退量不再次计费。所有测试对象仍是已暴露开发对象。', '',
   '## 动作误差排序器：三个种子','',
   '|种子尾号|所选动作误差|frozenJ|decoder exact|field updated|Broyden|','|---|---:|---:|---:|---:|---:|']
  for r in operator:lines.append(f"|{r['seed']%1000}|{r['selected_material_action_error_parent_mean']:.6f}|{r['frozenJ_material_action_error_parent_mean']:.6f}|{r['decoder_exact_material_action_error_parent_mean']:.6f}|{r['field_updated_material_action_error_parent_mean']:.6f}|{r['one_secant_Broyden_material_action_error_parent_mean']:.6f}|")
  lines += ['', '误差为 H=I+JᵀJ 的实际 GN 材料动作教师误差，先按父对象再平均。排序器本身需要候选点已付 Full 响应/场特征（8 RHS），不是无需候选前向的替代器；误差标签或选择准确性不是部署证书。', '',
   '## 冻结几何核：每父对象的有限运行 ESS','',
   '4 核×8 对象×1,024 步实际增加 262,144 RHS、32,768 Full 状态/LU。每个核/对象 8,192 RHS；此前教师与公共上下文取得单独披露，不在各核重复计费。ESS 是 8 个预定材料函数的 initial-positive-pair 有限运行估计，常量函数为未评估，不当作无限 ESS 或 MCSE=0；没有共同后验参考，不能跨不同父对象合并。下表给各对象/核的有限函数 ESS 最小值，完整向量见 CSV/JSON。','',
   '|父对象尾号|isotropic|weak|bracket enriched|random enriched|公共锚点 RHS|几何/Heun 拥有 RHS|','|---|---:|---:|---:|---:|---:|---:|']
  for own in owner:
   match={r['kernel']:r['finite_function_ESS_min'] for r in geometry if r['parent']==own['parent']}
   lines.append(f"|{int(own['parent'].split('_')[-1])}|{match['isotropic']:.3f}|{match['weak']:.3f}|{match['weak_bracket']:.3f}|{match['weak_random']:.3f}|{own['actual_shared_online_context_RHS']}|{own['paid_first_geometry_all_configurations_RHS']+own['paid_Heun_all_configurations_RHS']}|")
  lines += ['', '几何/Heun 拥有 RHS 是同一父对象原始研究的全部配置，不是某个冻结核的边际建造费；选用的 rank4/66 中央 FD 480 RHS 已在其中，不加第二遍。核的线性期望材料 RMS=.002；不是本次任务尺度 .02 的新实验。所有核依然 UNQUALIFIED，非后验参考。', '',
   '## 前瞻 locked 四臂包装器：只准备','',
   '`tools/run_r4_interval_chain_locked.py`，名称 interval_chain_locked_v1，未来 Book IV，R4L1_00000–00007，seed20261025171，8×1,024 步 pCN beta=.002；Full/Richardson/globalGMRES/固定精度 globalGMRES 保留相同缓存、存储与小 LS 收费权利及 300,000 组合 RHS 子预算。每臂 1,200 秒/3GiB。', '',
   '只有 FROZEN_LOCKED、在线 locked 上下文/缓存/观测、已关闭 StageIV 上下文收据及显式匹配父释放全部具备才可运行；当前包装器没有等待后自动启动逻辑。缺失输入在 Book 前停止。fast 引擎仅增加默认仍为原 II/DEV 的 stage/context/cache hooks；已封账 source_snapshot 与旧 global 入口不改写。准备与缺失门禁测试通过不等于真实 locked 实验完成。此次物理 RHS、状态、LU 与神经更新都是 0。']
  summary.write_text('\n'.join(lines)+'\n');print('Factual summary and inert locked wrapper prepared; physics0,training0.')
 except BaseException as exc:status='ERROR';error=repr(exc);raise
 finally:b.close(status,error)

if __name__=='__main__':main()
