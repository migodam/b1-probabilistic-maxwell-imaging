#!/usr/bin/env python3
"""Close preparation evidence only: no neural fit, model forward or Maxwell call."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[name]='1'
from pathlib import Path
import sys,json,subprocess,ast
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3r4.core import OUT,Book,write

def main():
 audit=OUT/'AUDIT/forward_jvp_surrogate';b=Book('II','forward_jvp_prepare_finalize_v1',max_wall=900,analysis_only=True)
 status='COMPLETE';error=None
 try:
  for source in ('src/p3r4/forward_jvp_surrogate.py','tools/train_forward_jvp.py','tools/prepare_r4_forward_jvp.py'):
   ast.parse((ROOT/source).read_text())
  test=subprocess.run([sys.executable,'-m','unittest','discover','-s','tests/p3r4','-p','test_forward_jvp_algebra.py','-v'],cwd=ROOT,env={**os.environ,'PYTHONPATH':str(ROOT/'src')},capture_output=True,text=True)
  (audit/'FINAL_ALGEBRA_TEST_LOG.txt').write_text(test.stdout+test.stderr)
  if test.returncode:raise RuntimeError('pure packing/JVP algebra tests failed')
  prep=json.loads((OUT/'II/jobs/forward_jvp_prepare_v1/CLOSED.json').read_text())
  manifest=json.loads((OUT/'results/II/forward_jvp_pack_v1/MANIFEST.json').read_text())
  scope=dict(status='READY_PREPARED_NOT_TRAINED',syntax_files=3,algebra_tests=3,Maxwell_RHS_new=0,Maxwell_states_new=0,physical_LU_new=0,NN_optimizer_steps=0,
   model_forward_executed=False,NN_accuracy='NOT_RUN',inference_wall='NOT_MEASURED_UNTRAINED',teacher_geometry_owner_RHS=459712,
   teacher_J_state_owner_RHS=414040,paid_pairs=manifest['total_pairs'],initial_preparation_receipt=prep,
   checks=['complex packing reversible and Phi units','fixed raw-coordinate unit Rademacher columns and parameter count','normalization chain rule using literal linear tensors, no trained network'],
   unexecuted=['training backward through neural JVP','trained response/JVP errors','inference timing','online Full fallback'])
  write(audit/'TESTED_SCOPE.json',scope)
  (audit/'REPORT_ZH.md').write_text('''# 局部前向／JVP 代理：准备完成，尚未训练

本模块直接预测局部非线性观测响应，已与此前的动作误差排序器分开。现有资产包含 10,063 个已付费配对：66 维 probit 坐标 5,029 个，固定属性的 64 维场坐标 5,034 个。准备过程中没有新 Maxwell RHS、状态或物理 LU，也没有神经优化步骤。历史 Stage III 几何教师的完整拥有成本仍为 459,712 RHS，其中保存 Full-J 状态对应 414,040 RHS；公共上下文取得成本另外保留。复用不意味着冷启动标签免费。

## 模型与观测单位

令 $w(F)=\\sqrt2\\,\\mathrm{pack}(F)/\\sigma$，原观测势为 $\\Phi=\\tfrac12\\|w(F)-w(y)\\|^2$。网络预测 $R(h)=w(F(x_0+h))-w(F_0)-J_0h$，总预测为 $w(F_0)+J_0h+\\widehat R(h)$。输入依次为 $[x_0,h,w(F_0),w(F_0)-w(y),J_0h]$。部署输入不含候选点 Full 观测、场、Jacobian、误差标签或真值。这里保留了精确锚点 $J_0$ 的前提：若未缓存，取得锚点前向和接收器导数本身需 8+32 RHS，不能被一次网络查询的计时掩盖。

结构是宽度 256 的两层 LayerNorm/SiLU MLP，输出 512 维。用 $s_R[N(u(h))-N(u(0))]$ 保证零位移响应为零。66 维有 625,664 参数，64 维有 624,640 参数，均小于 0.8M。

只对原坐标候选位移 $h$ 求导，锚点上下文固定并从计算图分离。原始输入切向矩阵为 $A=[0;I;0;0;J_0]$，归一化后为 $\\mathrm{diag}(s_u^{-1})A$。网络响应导数为 $\\mathrm{diag}(s_R)DN\\mathrm{diag}(s_u^{-1})A$；减去的零位移网络分支与 $h$ 无关。两列全局固定 Rademacher 方向 $V_{ik}=\\pm1/\\sqrt d$ 的种子为 2026102518。教师为 $(J_1-J_0)V$，总响应 JVP 要再加 $J_0V$。保存的 Full-J 已足够构造标签；未保存的 $T_1$ 不需要补求解。

## 冻结比较与数据所有权

父对象 0–15 为训练，16–23 为验证，24–31 为算子测试，全部属于已暴露开发数据，不能称盲测。每个父对象文件保存来源 Full-J、端点文件路径和在线观测。标准化只用训练父对象等权矩，输入标准差下限 1e-6，输出关于零的 RMS 下限 1e-4。链式法则矩阵示例、两列 V 和标准化数组均已落盘。

值监督与导数监督使用相同模型、父对象均衡采样、批量 64、学习率 3e-4 和种子 2026102501/2026102502。损失为标准化响应均方误差；导数组另加权重 1 的两个标准化 JVP 均方误差平均，值组权重为 0。初始协议 1,000 更新，显式 CLI 上限 3,000；每作业受 Book 的 900 秒受控时间限制。检查点均用同一验证响应 MSE 选择，训练阶段不读取算子测试折。CPU 单进程，未授权启动训练。

未来评估保存预测、白化响应 L2 误差、原始 Phi 误差、两方向 JVP 误差和相同样本上的 frozen-J 基线。两项经验门限是最大响应 L2 误差不超过 .05、最大绝对 Phi 误差不超过 .05。失败时明确要求每候选 8 RHS Full 回退；代码不会偷偷执行该回退。通过也只是开发样本上的经验筛查，不是逐点误差证书，更不能推出后验或采样正确性。

计时器覆盖热模型和固定锚点下的单候选位移组装、归一化、两次锚定网络计算和 float64 线性重建；冷加载、上下文准备及 Full 回退分别列出。未训练时没有报告推理时间，也没有速度收益结论。

## 实际验证与文件

三项纯代数测试通过：复数打包可逆且 Phi 单位一致、固定方向单位范数及参数上限、用字面线性张量检查标准化链式法则。没有执行神经网络拟合、训练后的反向 JVP、在线回退或新物理验证。训练器语法检查通过；未来运行需显式 --train 或 --evaluate 与唯一 job。未训练权重不存在，误差与推理速度均为 NOT_RUN。

主要资产：`results/II/forward_jvp_pack_v1/`，冻结协议 `FROZEN_PROTOCOL.json`，模块 `src/p3r4/forward_jvp_surrogate.py`，有界 CLI `tools/train_forward_jvp.py`。`TESTED_SCOPE.json` 与两份测试日志限定已测试范围。初次准备已关闭收据：受控 37.60088654095307 秒，CPU 18.782833 秒，峰值 RSS 280,461,312 字节，新增物理 0，优化步骤 0；本次最终整理另记独立零物理分析收据。
''')
 except BaseException as exc:status='ERROR';error=repr(exc);raise
 finally:b.close(status,error)
 print('Forward/JVP preparation finalized; physics0, optimizer0.')

if __name__=='__main__':main()
