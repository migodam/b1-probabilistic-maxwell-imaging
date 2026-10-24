#!/usr/bin/env python3
"""Final descriptive diagnosis; reads sealed results only, no physics or GT files."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
import sys,json,csv
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from p3r4.core import OUT,Book,write
from p3r4.evaluation import checkpoint_identity

def read(path):return json.loads(Path(path).read_text())
def main():
 dest=OUT/'AUDIT/L1_optimization_final_diagnosis'
 if dest.exists():raise FileExistsError('preserve existing diagnosis')
 b=Book('IV','L1_optimization_final_diagnosis_v1',max_wall=900,analysis_only=True);status='COMPLETE';error=None
 try:
  paths=dict(scores=OUT/'results/offline_only/L1_optimization_addenda_GT_v1/RESULTS.json',plan=OUT/'results/offline_only/L1_optimization_addenda_GT_v1/CHECKED_GATE_AND_FIXED_ANALYSIS_PLAN_BEFORE_GT.json',nonGT=OUT/'AUDIT/L1_optimization_addenda_v2/NON_GT_FACTS.json',slsqp=OUT/'AUDIT/ordinary_slsqp_diagnosis/RESULTS.json',callback=OUT/'AUDIT/ordinary_slsqp_diagnosis/MIN_PAID_QUERY_AND_CALLBACK_REVIEW.json',receipt=OUT/'IV/jobs/L1_optimization_addenda_isolated_GT_v1/CLOSED.json')
  scores,plan,ng,slsqp,callback,receipt=(read(paths[k]) for k in ('scores','plan','nonGT','slsqp','callback','receipt'))
  if receipt['status']!='COMPLETE' or scores['choices_changed'] or scores['new_model_or_algorithm_selection']:raise ValueError('only closed fixed-choice scores')
  for ident in plan['choices_seals']+[plan['release_identity']]:
   if checkpoint_identity(ident['path'])!=ident:raise ValueError('source changed')
  if plan['scoring_pid'] in plan['source_pids']:raise ValueError('not isolated scoring PID')
  with np.load(OUT/'results/offline_only/L1_optimization_addenda_GT_v1/PER_PARENT_LOSSES.npz',allow_pickle=False) as z:
   actions=z['actions'].tolist();loss=z['full_mse'].copy();parents=z['parent'].tolist()
  assert loss.shape==(64,17) and len(set(parents))==64
  rows=[]
  for comp in scores['comparisons']:
   d=loss[:,actions.index(comp['reference'])]-loss[:,actions.index(comp['candidate'])]
   rows.append(dict(comp,parent_improved=int((d>0).sum()),parent_tied_exactly=int((d==0).sum()),parent_worse=int((d<0).sum())))
  dest.mkdir(parents=True)
  with (dest/'FIXED_PAIRED_COMPARISONS.csv').open('w',newline='') as f:
   w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
  facts=dict(status='FINAL_DESCRIPTIVE_DIAGNOSIS_COMPLETE',parents=64,actions=17,fixed_comparisons=48,new_RHS=0,new_states=0,new_LU=0,new_fits=0,direct_GT_file_reads=False,
   fixed_choice_scores=scores,nonGT_diagnosis=ng,saved_trial_diagnosis=slsqp,callback_source_review=callback,comparisons_with_parent_counts=rows,
   source_identities={k:checkpoint_identity(p) for k,p in paths.items()},scoring_pid=plan['scoring_pid'],source_pids=plan['source_pids'],
   choices_sealed_before_truth=True,primary45risk26density_unchanged=True,exploratory_only=True,global_optimum_claim=False,
   prior_metadata_failure_retained='L1_optimization_addendum_nonGT_report_v1 ERROR was duplicate report-key merge only; original choices/sources untouched; v2 report succeeded')
  write(dest/'FACTS.json',facts)
  lines=['# 第一64个锁定对象：优化失败诊断与封存补充的最终事实','',
   '本报告只读封存结果和已完成的隔离评分，不读取真值文件、不重选动作。两项补充均在主方法冻结后、首次GT释放前登记并封存，是探索性诊断；原320个代表元、45 risk与26 density主输出不变。统计单位始终是64个对象。','',
   '原普通SLSQP不是强额外优化证据。64个对象均返回原锚点，但1024次已付费试探中，每个对象都曾查询到更低Phi；返回值与已付费最小Phi之差中位0.00207002。198条callback记录均不满足原噪声球，且没有一条比锚点Phi更低。已安装SciPy Python包装器只在major-iteration计数增加时触发callback，没有下降或可行性谓词；因此callback不能被当作已接受下降步的证明。此处未对底层C实现的精确线搜索时序作额外断言。','',
   '首试探重放与负原始数据梯度一致（最大误差1.423e-13），步长中位2.6903，首试探Phi中位4632.77，超过原阈值282.873792。这与未缩放方向的数值困难相符；不能由预算停止推断没有更好的物理拟合。best-query普通臂保留已验证的低Phi查询，64/64降低返回Phi；四个mu臂共256例没有更低的可行材料目标候选。它是BEST_VALIDATED_QUERY_CANDIDATE，不是优化收敛或线搜索接受证书。','',
   'Full-J信赖域是独立、预登记的正确度控制。两个臂分别最小化Phi与Phi+.5||x||²；Phi=.5||r||²，r为正确复数打包和白化后的残差。每次接受后付费刷新J32，拒绝试探Full8也保留。λ0在64/64降低Phi，变化中位−1.149918；GN decrement中位5.52331→4.89802，仍不支持数据目标已经驻定。λ1目标在64/64下降，变化中位−5.434e-7；初/末GN decrement中位.00107859/.0000383092。','',
   'λ1小decrement表示多数点在GN曲率尺度上变化很小；不能泛化到全部对象：初/末最大decrement为1.68388/.542270，末精确总梯度无穷范数最大2.54760。GN矩阵不是完整非线性Hessian，所有128臂都由同一局部预算预检停止，不证明全局MAP、最优性或MMSE。λ0末精确梯度中位105.942亦不应隐藏；目标下降、梯度范数和GN decrement是不同量。','',
   '|固定动作|全场MSE均值|全场MAE均值|','|---|---:|---:|']
  for r in scores['means']:lines.append(f"|{r['action']}|{r['full_mse']:.10g}|{r['full_mae']:.10g}|")
  lines+=['','下表仅摘录已登记的两个TR及best-query普通臂对原锚点比较；正改进表示原锚点MSE减候选MSE。其余48个固定比较及全部64对象保留在配套CSV/原评分中。区间为按对象配对bootstrap的探索性描述，不是多重检验控制或主结论。','','|固定候选|MSE平均改进|描述性95%区间|改善/相同/恶化对象|','|---|---:|---|---|']
  for r in rows:
   if r['reference']=='original_anchor' and r['candidate'] in ('TR_lambda0','TR_lambda1','best_query_ordinary'):
    lines.append(f"|{r['candidate']}|{r['improvement_absolute']:.9g}|[{r['absolute_CI'][0]:.9g}, {r['absolute_CI'][1]:.9g}]|{r['parent_improved']}/{r['parent_tied_exactly']}/{r['parent_worse']}|")
  lines+=['','在这个已固定的64对象子集中，λ0降低数据失配却提高平均材料误差；λ1与原锚点的平均材料误差接近。此事实不能支持继续调参、选择lambda或改变已冻结的输出。材料目标代表元与raw/clipped估计量的不同角色和实现约束保留，不能仅依据这里的GT排序作新的在线选择。','',
   '成本：原代表元搜索80640 RHS；新TR30752 RHS（Full8为9184，接受后J32为21568，合计1148状态/LU），不超过32768上限；best-query新增物理0。公共前64上下文实际拥有240304 RHS，初始cached Full+J冷参考40 RHS/对象，不重计也不隐去。TR受控22.423865秒、CPU21.036236秒；隔离评分受控.884132秒、CPU.422254秒、GPU0。当前报告新增Maxwell、RHS、LU和训练均0。原raw/PPC诊断另有独立receipt，不归入优化收益。','',
   '首次评分PID83920与TR82674、best-query82788、原代表元80928不同；全部封存与release身份在首次truth访问前核验并写固定计划。元数据报告v1的重复键合并错误及ERROR账保留，v2只修复报告生成，不改动候选、算法或原结果。','',
   '证据入口：`../L1_optimization_addenda_v2/DIAGNOSIS_REPORT_ZH.md`、`../ordinary_slsqp_diagnosis/MIN_PAID_QUERY_AND_CALLBACK_REVIEW.json`、`../L1_OPTIMIZATION_ADDENDA_GT_APPEND_REPORT_ZH.md`。全部来源绝对路径和身份记于FACTS.json。最终科学解释由主线程保留。']
  (dest/'REPORT_ZH.md').write_text('\n'.join(lines)+'\n')
  print(json.dumps(dict(status=facts['status'],dest=str(dest),new_RHS=0,anchor_comparisons=[r for r in rows if r['reference']=='original_anchor' and r['candidate'] in ('TR_lambda0','TR_lambda1','best_query_ordinary')]),ensure_ascii=False),flush=True)
 except BaseException as exc:status='ERROR';error=repr(exc);raise
 finally:b.close(status,error)
if __name__=='__main__':main()
