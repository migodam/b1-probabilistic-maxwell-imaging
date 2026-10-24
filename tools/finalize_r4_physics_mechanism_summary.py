#!/usr/bin/env python3
"""Finish prose from saved fact extraction; no fit, physics or extraction replay."""
from pathlib import Path
import sys,json
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3r4.core import OUT,Book,write
from p3r4.evaluation import checkpoint_identity

def main():
 summary=OUT/'AUDIT/PHYSICS_MECHANISM_SUMMARY.md'
 if summary.exists():raise FileExistsError('preserve report')
 b=Book('II','physics_mechanism_summary_finalize_v1',max_wall=60,analysis_only=True);status='COMPLETE';error=None
 try:
  facts=json.loads((OUT/'AUDIT/PHYSICS_MECHANISM_FACTS.json').read_text());dino=facts['DINO'];operator=facts['operator_ranker'];geometry=facts['frozen_geometry'];owner=facts['shared_parent_acquisition']
  prepared=OUT/'AUDIT/locked_samecache/FROZEN_PROTOCOL.json';template=OUT/'configs/INTERVAL_CHAIN_LOCKED_RELEASE_TEMPLATE.json'
  if template.exists():raise FileExistsError('preserve inert release template')
  write(template,dict(released=False,experiment='interval_chain_locked_v1',campaign_id='P3LR4-20261025',prepared_protocol_identity=checkpoint_identity(prepared),
   context_receipt='<absolute StageIV contexts_locked CLOSED.json after acquisition COMPLETE>',
   source='tools/run_r4_interval_chain_locked.py',roster=[f'R4L1_{i:05d}' for i in range(8)],reference='independent Full with common fixed pCN innovations/U;after-close compare only',
   scope='inert template only; parent must explicitly create matching released:true document after FROZEN_LOCKED and online acquisition'))
  lines=['# 物理机制控制：封账事实表','',
   '本页汇总已完成的 6 次局部前向/JVP 拟合、3 次动作误差排序器拟合，以及前 8 个 DEV 的 4 种冻结几何核。此次提取没有物理查询、拟合或真值读取，也不作最终 GO 判断。逐种子/父对象机器事实与来源路径见 PHYSICS_MECHANISM_FACTS.json 和 physics_mechanism_tables/ 的四份 CSV。','',
   '共享几何教师拥有成本是 459,712 RHS（保存 Full-J 状态部分 414,040 RHS；10,063 个独特配对），种子之间不重复计费，复用也不意味着冷标签免费。公共在线锚点取得另外列出。若部署未缓存，单锚点 Full+接收器 J 需 8+32 RHS。','',
   '## 六次局部前向/JVP 拟合','',
   '|监督|种子尾号|更新/所选|最大白化 L2|95分位 L2（描述）|最大绝对 Phi 误差|原验证/算子测试门限|','|---|---:|---:|---:|---:|---:|---|']
  for r in dino:lines.append(f"|{r['mode']}|{r['seed']%1000}|{r['optimizer_updates']}/{r['selected_update']}|{r['max_white_L2']:.4g}|{r['white_L2_q95_descriptive']:.4g}|{r['max_absolute_Phi_error']:.4g}|{r['original_validation_all_point_gate']}/{r['operator_test_all_point_gate']}|")
  lines += ['', '实际保存的通过条件是所有点最大白化 L2≤.05 且最大绝对 Phi 误差≤.05；六次原验证和算子测试门限全部未通过。这里没有增加“95%通过”规则，95 分位仅描述误差分布。每次算子测试为 8 个已暴露开发父对象、790 配对；热查询约 79–83 微秒，但失败门限仍要求每候选 8 RHS Full 回退。每模型反事实全回退 6,320 RHS 未在此执行或计费；不能用短推理时间宣称前向加速。精确概率用途仍需 Full 校正。','',
   '## 三次动作误差排序器拟合','',
   '|种子尾号|所选材料动作误差|frozenJ|decoder exact|field updated|Broyden|','|---|---:|---:|---:|---:|---:|']
  for r in operator:lines.append(f"|{r['seed']%1000}|{r['selected_material_action_error_parent_mean']:.6f}|{r['frozenJ_material_action_error_parent_mean']:.6f}|{r['decoder_exact_material_action_error_parent_mean']:.6f}|{r['field_updated_material_action_error_parent_mean']:.6f}|{r['one_secant_Broyden_material_action_error_parent_mean']:.6f}|")
  lines += ['', '每次完成 6,000 更新，误差先按测试父对象再平均。教师是 H=I+JᵀJ 对应的实际 GN 材料动作。排序器输入需要候选点 Full 响应/场特征，冷特征费为每候选 8 RHS；它不是不查询候选前向的响应替代器，经验误差或刷新选择不是部署证书。该表仅陈述与四个实际基线的数值关系。','',
   '## 冻结几何核：逐父对象有限 ESS','',
   '4 核×8 对象×1,024 步实际增加 262,144 RHS、32,768 Full 状态/LU；每个核/对象 8,192 RHS。实际整作业受控时间 '+f"{facts['frozen_geometry_controlled_seconds']:.3f}"+' 秒。ESS 是 8 个预定材料函数的 initial-positive-pair 有限运行估计；常量函数为未评估，不是无限 ESS 或 MCSE=0。无共同后验参考，不跨父对象合并。下表为各核有限函数 ESS 的最小值，完整向量见 FROZEN_GEOMETRY_PER_PARENT_KERNEL.csv。','',
   '|父对象尾号|isotropic|weak|weak bracket|weak random|公共锚点 RHS|原几何+Heun 拥有 RHS|','|---|---:|---:|---:|---:|---:|---:|']
  for own in owner:
   match={r['kernel']:r['finite_function_ESS_min'] for r in geometry if r['parent']==own['parent']}
   lines.append(f"|{int(own['parent'].split('_')[-1])}|{match['isotropic']:.3f}|{match['weak']:.3f}|{match['weak_bracket']:.3f}|{match['weak_random']:.3f}|{own['actual_shared_online_context_RHS']}|{own['paid_first_geometry_all_configurations_RHS']+own['paid_Heun_all_configurations_RHS']}|")
  lines += ['', '公共上下文、原几何与 Heun 成本每父对象只拥有一次，四核共享；几何拥有 RHS 是同对象研究的全部配置，不是单冻结核的边际构造费。所用 rank4/66 中央 FD 的 480 RHS 已在拥有费用内，不再计第二遍。线性期望材料 RMS=.002；这些冻结核没有提供任务尺度 .02 的新证据。所有 ESS/MCSE 和核均保持有限运行 UNQUALIFIED。','',
   '## 前瞻 locked 四臂：已准备，未释放','',
   '入口 tools/run_r4_interval_chain_locked.py，未来名称 interval_chain_locked_v1、Book IV，R4L1_00000–00007，seed20261025171，8×1,024 步 pCN beta=.002。四臂 Full/Richardson/globalGMRES/固定精度 GMRES 共享原协议的缓存/存储权利与小 LS 收费；组合上限 300,000 RHS（65,536/78,154/78,155/78,155），各臂 1,200 秒、3GiB。目标与门限不变。','',
   'AUDIT/locked_samecache/FROZEN_PROTOCOL.json 已冻结；configs/INTERVAL_CHAIN_LOCKED_RELEASE_TEMPLATE.json 为 released:false。只有 FROZEN_LOCKED、在线 locked 上下文/缓存/原观测、已关闭 StageIV 上下文收据和精确匹配协议身份的父释放全部具备才可运行。没有等待后自动启动逻辑，缺失输入在 Book 前停止；当前 preflight 明确缺 FROZEN_LOCKED。','',
   'fast 引擎只添加默认仍为 StageII/DEV 的 stage/context/cache hooks；旧 global 入口与已封账源码快照不改。两项零物理包装器检查通过不意味着 locked 物理实验完成。paired physics VI 的执行源码没有改动。','',
   '事实提取首个分析作业在写正文时使用了错误的核展示键 bracket_enriched，报 KeyError；实际 CSV/JSON 提取与零物理测试已完成且保留，原 ERROR 收据不改。最终整理只读取这些事实，将展示键改为实际 weak_bracket/weak_random，没有重做任何物理或训练。相关成本见 II/jobs/physics_mechanism_summary_and_locked_prepare_v1/CLOSED.json 与 physics_mechanism_summary_finalize_v1/CLOSED.json。此次新增 RHS、状态、LU、神经优化均为 0。']
  summary.write_text('\n'.join(lines)+'\n')
  write(OUT/'AUDIT/physics_mechanism_tables/FINALIZATION.json',dict(status='COMPLETE',prior_partial_analysis_status='ERROR_REPORT_KEY_ONLY',saved_facts_reused=True,no_closed_results_overwritten=True,new_RHS=0,new_training_steps=0,release_status='NOT_RELEASED',summary=str(summary)))
  print('Summary finalized from saved facts; locked wrapper inert;physics0,training0.')
 except BaseException as exc:status='ERROR';error=repr(exc);raise
 finally:b.close(status,error)

if __name__=='__main__':main()
