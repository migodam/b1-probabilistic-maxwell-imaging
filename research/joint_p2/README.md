# B1 × C1 Joint P2：probabilistic material hypothesis package

实验阶段20261013；公开交付2026-10-14。Fresh主概率参考与范围以 [ACCEPTED_REFERENCE_AUDIT](results/ACCEPTED_REFERENCE_AUDIT.json) 和 [全部fresh细化记录](results/ALL_FRESH_REFERENCE_REFINEMENT.json) 为准；旧13阶与中间结果保留为历史近似。

本目录保存同心球3Dvector-Mie材料prior中的多假设posterior研究。输出是互不重叠materialevents、事件质量、代表材料、条件不确定性、完整weightedparticles、剩余已建模质量与independentcontextupdate，不给continuousimagepoint非零probability。

先读 [START_HERE](START_HERE.md) 和 [独立中文报告](RESEARCH_REPORT_ZH.md)。可浏览图库入口为 [index.html](index.html)；概率package在 [results/packages](results/packages/)。full/proxytarget、模型外未知状态和数值误差均由JSON保留。C1density代码、checkpoint与seed-separatedresults位于 [现有C1仓库p2](https://github.com/migodam/c1-physics-state-generative-imaging/tree/main/p2)。

## 代码与资产

| 资产 | 内容 |
|---|---|
| `src/b1_joint/material.py` | continuousprior、strictvacuumdecoder、eventpartition、materialdistance与截面 |
| `likelihood.py` | realpacking、propertruncatedsharedgain、independentgeometrylikelihood |
| `reference.py` | 同密度evidence/moments、directCDF/quantile、same-priorphysicaladaptiveintegration |
| `correction.py` | full-minus-proxyeventintegralcontrolvariate、完整knownq与误差 |
| `smc.py` | generictemperedSMC、full-targetdelayedacceptance与原始particles |
| `proposal_correction.py` | matched64callMDN/full-MieSNIS与cold/warmproposal计时 |
| `package.py` | adaptiveTop-K、eventmedoid、forwardprediction与wholeweightedcontextupdate |
| `interface.py` | 无topology输入的`ImagingSession.decode/predict/infer/update_context`共同接口 |
| `accounting.py` / `forward.py` | 唯一物理预约、failure-inclusiveledger、cacheprovenance、历史账本隔离 |
| `tools/p2_*` | 必要probabilitycontrols、continuousaudit、handoff、gallery及受控jobmeter |

Python依赖为NumPy、SciPy、Matplotlib，C1需要CPU PyTorch；无CUDA基础设施、大型NN或DeepSeek。P2不存在locked-data入口。所有科研运行必须在批准的账本与wall预算内，不能新checkout后重置4000。

## 阅读和零物理重放

公开代码与summary保存在Git；较大teacher/rawarrays与checkpoint证据按Release归档。恢复证据包的原相对目录后，以B1与C1clone根目录设置路径，示例中的`/path/to`需要换成真实本地路径。

```bash
export B1_P2=/path/to/b1-probabilistic-maxwell-imaging
export C1_P2=/path/to/c1-physics-state-generative-imaging
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONPATH="$B1_P2:$C1_P2/p2"

# 冻结模型的raw评价写到新replay目录，不覆盖原始fresh结果。
python -m c1p2.evaluate "$B1_P2/research/joint_p2/results/C1_FINAL_DATA.npz" \
  --models "$C1_P2"/p2/results/models/*-2610132[123].pt --output "$C1_P2/p2/replay/evaluation"
python -m c1p2.summarize "$C1_P2/p2/replay/evaluation/EVALUATION_SUMMARY.json" \
  --output "$C1_P2/p2/replay/summary"
python -m c1p2.continuous_audit "$B1_P2/research/joint_p2/results/FRESH_CONTINUOUS_REFERENCE_ACCEPTED.npz" \
  --evaluation "$C1_P2/p2/replay/evaluation" --output "$C1_P2/p2/replay/continuous"
```

这些命令复用已支付physicsbank，不产生新增Maxwell调用，但仍消耗CPUwall；重放不是独立新validation。必要概率检查可读 [E0_AUDIT.md](E0_AUDIT.md) 或用`tools/p2_probability_controls.py`与C1的`p2/tests/test_density_probability.py`重放；不要求pytest、SHA或smoke。

## 从观测到多材料输出

`ImagingSession`直接接受32个实化观测通道及固定noise；不接受truth、topology或teacher作为onlinefeature。下面用已公开的measurement说明入口，新观测同样走该接口。

```python
import json
from src.b1_joint.interface import ImagingSession

observed = json.load(open('research/joint_p2/results/packages/solid.json'))
session = ImagingSession('research/joint_p2/results/COMMON_RESPONSE_BANK.npz')
package = session.infer(observed['measurement'], observed['noise_sigma'])
for hypothesis in package['hypotheses']:
    print(hypothesis['event_id'], hypothesis['posterior_mass'])
    material = session.decode(hypothesis['representative_latent'])
```

该通用入口用已支付bank做连续proxy积分，并明确标为新观测尚未自动full-Miecorrection/convergence认证。`predict`仅在传入显式预算控制的`MaterialForward`时调用新fullsolve，不提供免费physics后门。已有七个研究package的selectedfullendpoint与两个fullresidualreference按原ledger保存，不由此示例重新计费。

## Whole posterior context接口

```python
import numpy as np
from src.b1_joint.package import update_context
from src.b1_joint.likelihood import context_loglik

a = np.load('research/joint_p2/results/packages/large_cavity_samples.npz')
result = update_context(a['latent'], a['weights'],
                        lambda x: context_loglik(x, c=0.02, sd=0.01))
```

此接口接受声明的independentcontextlikelihood，返回全部eventquality与weights，而非只对展示Top-K强行归一化。上例`c=.02`仅演示API，具体科研sensor值以冻结manifest和下游报告为准。若context与原measurement不条件独立，必须提供正确`p(c|x,y)`，不能重复使用y的证据。

## 原始physicscampaign与预算

`src.b1_joint.run`分phase保存bank、mechanism、continuousreference、teacher、fullcorrection、SMC与freshtruth；相同completedjobID拒绝重用。历史old384/384保持readonly，P2总上限4000、controlledwall7200s、GPUtraining0。完整reproduction涉及物理调用，不应无授权自动重跑；已完成data/ledger足以审查概率与成像结果。原失败、negativecases和数值修正保存于 [FAILURES.md](FAILURES.md)，不靠重跑遮盖历史缺口。

概率保证范围见 [POSTERIOR_REFERENCE.md](POSTERIOR_REFERENCE.md)、[THEORY_STATUS.md](THEORY_STATUS.md)和 [THEORY_HANDOFF.md](THEORY_HANDOFF.md)。最后一份证明为何eventmass和代表图不足以一般更新，何时条件样本/密度需要保留，以及TV/regret界与KL链式分解的前提；不是本项目已获得uniform认证或原创定理。仅凭normalization、exactMH或高ESS不能声明全模态、现实calibration或一般三维支持。
