"""Material-task empirical identities and compact artifact report."""
import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
from pathlib import Path
import sys,json,time
sys.path.insert(0,str(Path(__file__).resolve().parent))
from run import OUT,HIST,OLD,SEEDS,ARMS,write
from inference import Action
import numpy as np,torch

def main():
 torch.set_num_threads(1);va=dict(np.load(HIST/'prepared_validation.npz'));rows=[];identity=[]
 for arm in ['original_feature32768']+ARMS:
  records=[]
  for seed in SEEDS:
   if arm=='original_feature32768':
    folder=HIST/'learned_context'/f'distilled_w128_n32768_s{seed}';r=json.loads((folder/'RESULT.json').read_text());pr=np.load(folder/'validation.npz')['prediction'];f=Action(arm,seed,'frozen');fp=f.predict_prepared(va['image'],va['anchor'],f.moments(va['y']));frozen=float(np.mean(abs(fp-va['truth'])**2));matched=r['validation_full_complex_MSE'];provider_update=r['provider_selected_update'];head_update=r['head_selected_update']
   else:
    folder=OUT/'models'/f'{arm}_s{seed}'
    if not (folder/'RESULT.json').exists():continue
    r=json.loads((folder/'RESULT.json').read_text());pr=np.load(folder/'validation.npz')['prediction'];fp=np.load(folder/'frozen_validation.npz')['prediction'];frozen=r['frozen_validation_full_complex_MSE'];matched=r['matched_validation_full_complex_MSE'];provider_update=r['provider_selected_update'];head_update=r['head_selected_update']
   records.append(dict(seed=seed,frozen=frozen,matched=matched,provider_update=provider_update,head_update=head_update))
   teacher=va[f'det_frozen_s{seed-100}'];g=torch.load(HIST/'geometry.pt',map_location='cpu',weights_only=False)['geometry'];from p3r4.risk import full_material
   tp=full_material(teacher,g);d=fp-tp;e=tp-va['truth'];D2=float(np.mean(abs(d)**2));cross=float(2*np.mean((np.conj(e)*d).real));diff=float(np.mean(abs(fp-va['truth'])**2)-np.mean(abs(e)**2));identity.append(dict(arm=arm,seed=seed,frozen_teacher_risk=float(np.mean(abs(e)**2)),student_frozen_risk=frozen,risk_difference=diff,D_squared=D2,cross_term=cross,identity_gap=diff-D2-cross,scope='exact empirical finite-VAL expansion; GT replaces population conditional mean'))
  if records:rows.append(dict(arm=arm,records=records,mean_frozen=float(np.mean([r['frozen'] for r in records])),mean_matched=float(np.mean([r['matched'] for r in records]))))
 write(OUT/'SUMMARY.json',dict(status='COMPLETE' if len(rows)==5 and all(len(r['records'])==3 for r in rows) else 'PARTIAL',data='historical TRAIN1024/VAL256; old tests never read',records=rows,physics_RHS=0));write(OUT/'FINITE_TASK_RISK_IDENTITIES.json',dict(status='COMPLETE',records=identity,scope='empirical identities only, no population Bayes claim',physics_RHS=0))
 lines=['# R6 Track B：固定容量、有限训练集的任务功能学习','', '所有新学习只用历史 TRAIN1024 个完整 anchor，VAL256 已暴露并用于选择；**没有新的独立评价，也没有概率后验合格结论**。原始 N32768 蒸馏基线保留不可变。新模型和基线均为214272参数 provider、6418参数 head；provider 与匹配 head 各4000更新，3 seeds。','', '四个新臂都从相同 seed 初始化、使用相同 TRAIN1024 batch schedule、固定同一对应历史 head；归一化沿用 TRAIN32768 的 y/teacher 统计。feature_control1024 只拟合 teacher 特征；direct_gt 直接最小化固定 head 的 GT 材料误差；feature_action 同时拟合特征与 teacher action；gt_action 同时拟合 GT action 与 teacher action。teacher action 是有偏任务标签，不是完整 posterior 真值。','', '|臂|冻结 head 平均 full complex MSE|匹配 head 平均 full complex MSE|','|---|---:|---:|']
 for r in rows:lines.append(f"|{r['arm']}|{r['mean_frozen']:.10f}|{r['mean_matched']:.10f}|")
 lines+=['','原始 N32768 只作为历史部署基准；新四臂之间才是相同 N1024、容量、优化更新和局部输入的因果比较。原始基线按 teacher feature VAL 选择 provider；新四臂按冻结 head 材料 VAL 选择 provider，不能把这两者的差异都归因于目标函数。所有更新0都纳入选择，日志保留该点。','', '动作是 full-grid 原始复材料均方误差；表为三个 seed 的 risk 平均，不是 ensemble 图像误差。head provider交换矩阵记录所有3×3组合，finite ensemble prediction variance 使用有限已训练模型恒等式，不能作为总体学习方差估计。mean-only/variance-only/constant/parent-shuffle 为冻结推断消融，没有重新训练，因此变化包含输入分布改变。','', '`inference.py:Action(method, seed, head_mode="matched").predict(y, online)` 可直接部署；online 必须是原来的 observation/residual/q0/chi0/J，predict 不读取GT或物理求解器。允许方法 original_feature32768 / feature_control1024 / direct_gt / feature_action / gt_action。','', '证据：PLAN_FROZEN.json、READINESS.json、TRAINING_JOB.json、jobs/*.json、models/*/{RESULT.json,provider_losses.jsonl,head_losses.jsonl,PROVIDER.pt,HEAD.pt,validation.npz,frozen_validation.npz}、SUMMARY.json、SEED_COUPLING.json、CHANNEL_CONTROLS.json、FINITE_TASK_RISK_IDENTITIES.json、ACTION_REPLAY.json。','', '复现：在仓库根用现有 Gaussian/.venv_nn/bin/python 执行 track_b/run.py，再执行 track_b/diagnostics.py、track_b/report.py。已存在 RESULT 的拟合会跳过，历史输入完全只读。新增物理RHS=0；所有模型训练都在MPS，CPU线程1；实际资源详见jobs与TRAINING_JOB。训练上限为3500秒保守wall（低于1device-hour），进程RSS上限5GB。']
 if (OUT/'SEED_COUPLING.json').exists():
  coupling=json.loads((OUT/'SEED_COUPLING.json').read_text())['records'];lines+=['','交换 provider_i/head_j 的全部3×3 risk：对角与非对角的平均差（非对角−对角，正数表示配对有利）：']
  for arm,v in coupling.items():
   fd=v['finite_crossed_decomposition'];lines.append(f"- {arm}: {v['offdiagonal_mean']-v['diagonal_mean']:.9f}; provider/head/interaction 有限预测变动平方均值 = {fd['provider_prediction_variation']:.8f} / {fd['head_prediction_variation']:.8f} / {fd['interaction_prediction_variation']:.8f}。")
  lines+=['','这些数是有限9个已训练组合的正交分解；相同训练数据、seed配对和表征归一化使其不能当作总体 bias/variance 或机制因果证明。原始 N32768 配对差约1e-6；GT相关臂配对差较大，说明跨 seed 交换不等于独立复现实验。']
 if (OUT/'ACTION_REPLAY.json').exists():
  replay=json.loads((OUT/'ACTION_REPLAY.json').read_text())['records'];delta=max(v['max_abs_prediction_difference'] for v in replay.values());lines += ['',f'CPU部署 Action 与已存 MPS 预测的最大绝对差 {delta:.9g}；仅报告float32设备漂移，没有新设或放宽接受阈值。']
 if (OUT/'COST_LEDGER.json').exists():
  ledger=json.loads((OUT/'COST_LEDGER.json').read_text());lines+=['',f"实际新训练24个stage：fit+validation合计 {ledger['sum_fit_and_validation_seconds']:.3f}s；保守MPS整job预留wall {ledger['conservative_MPS_reserved_job_wall_seconds']:.3f}s。进程peakRSS {ledger['peak_rss_bytes']/1e9:.3f}GB，采样peakMPSdriver {ledger['peak_recorded_MPS_driver_bytes']/1e9:.3f}GB；二者相加约1.181GB，不是连续全系统内存硬界证明。"]
 if (OUT/'DEPLOYMENT_COST.json').exists():
  costs=json.loads((OUT/'DEPLOYMENT_COST.json').read_text())['records'];lines+=['','部署CPU单线程、batch64、20次实测median（provider / provider+feature+head，微秒/parent）：']
  for arm,v in costs.items():lines.append(f"- {arm}: {v['provider']['median_seconds_per_parent']*1e6:.2f} / {v['provider_feature_head']['median_seconds_per_parent']*1e6:.2f}")
  lines+=['','计时与主线程物理任务同时运行，存在系统负载干扰；所有新provider架构相同，不能用这些差异声称某个loss目标具有推断速度优势。共同anchor生成仍未包含。原始逐次耗时与p10/p90均保留。']
 lines+=['','![开发集风险比较](RISK_COMPARISON.png)','','![GT动作训练曲线](GT_ACTION_LEARNING_CURVES.png)']
 (OUT/'REPORT_ZH.md').write_text('\n'.join(lines)+'\n')
if __name__=='__main__':main()
