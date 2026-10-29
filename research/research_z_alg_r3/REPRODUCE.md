# R3复现：运行、验收、公开评分

## 环境与目录

使用Python、NumPy、SciPy、threadpoolctl；画图另需Matplotlib。实际版本列于`experiments/requirements.txt`和`audit/EXPERIMENT_RUNTIME.json`。本轮用现有Apple Silicon CPU环境，BLAS每个worker限制为一线程，没有安装新数值软件、使用GPU或训练网络。

从既有B1仓库进入：

```bash
git clone https://github.com/migodam/b1-probabilistic-maxwell-imaging.git
cd b1-probabilistic-maxwell-imaging/research/research_z_alg_r3
```

也可解压Release的`FULL_REPRODUCIBILITY.zip`后进入`research_z_alg_r3/`。完整包包含原始输出、冻结快照、公开观测及评分后公开的材料。Reading Kit是阅读/代码摘要，不含全部原始结果。

## 从一个未知观测输入求材料

```bash
python run_portfolio.py \
  --public public_data/confirm/NEWCONFIRM12_00.npz \
  --method guarded_globaldual \
  --output /tmp/r3_observation_only_prediction.json
```

同一输入、同一默认3秒总预算的强普通对照：

```bash
python run_portfolio.py \
  --public public_data/confirm/NEWCONFIRM12_00.npz \
  --method guarded_reciprocal \
  --output /tmp/r3_ordinary_prediction.json
```

默认就是冻结V7配置；显式配置可通过`--config`传JSON。两臂均包含同样的原材料Full-J保护阶段。所有输出文件必须为新路径。新观测只需相同NPZ输入字段，不需要GT、seed或评分器。运行在已公开数据上属于暴露后重放，不能再叫新盲测。

基础入口`run_solver.py`支持`fullj`、`coordinate`、`algebraic`、`globaldual`。完整V6配置在`experiments/HOLDOUT12_CONTRACT.json`，最终组合配置在`experiments/CONFIRM_V7_CONTRACT.json`。不要用默认基础入口重新解释历史各臂的参数。

## 零物理费用重现原冻结评分

```bash
python experiments/score_predictions.py \
  --split confirm \
  --freeze-manifest audit/CONFIRM_V7_FINAL_FREEZE.json \
  --truth-root truth_after_freeze \
  --disclosed-truth-replay \
  --output /tmp/r3_v7_replay_scores.json
```

V6对应把split换成`holdout`，manifest换成`audit/HOLDOUT_V6_FINAL_FREEZE.json`。本轮已经执行公开材料重放，全部V7评分行与原结果完全一致，新增物理RHS=0。该flag明确放宽Git克隆后的mtime身份要求，并检查已记录字节长度；不使用SHA，也不冒充严格的密码学数据认证。

原独立验证器为`experiments/validate_predictions.py`：只读取公开算子及冻结候选，独立NumPy解完整原矩阵，不调用求解器和约化核，不读取GT。若重新运行，会新增真实物理RHS；原4517次独立验证已经完整保存，无需为了读报告重算。

## 冻结身份与不可比的重放

V6八份源码以及20份预测、V7十二份源码/配置以及八份预测均有字节对照快照。V6所有结果保留，随后用于V7设计；V7四目标是在其独立冻结以后新生成。

最终方法受wall预算限制，机器负载、BLAS、解释器和硬件变化会改变预算结束时已到达的起点/候选，所以重放不保证逐位或逐候选相同。原预测、费用、时间轨迹和评分是主证据。只读评分重放则可精确再现同一冻结输出的指标。

开发阶段17次输出和当时配置全部保留；最早几版源码没有独立版本快照，因此不声称能用最终代码逐条重建早期重复修正等实现行为。主结论的V6/V7源码均独立冻结。首个生成失败消耗0RHS，其wall及早期部分评分wall未留存，账本明示缺口。

## 发布边界

`private_truth/`不入Git、不入压缩包。`truth_after_freeze/`仅包含10份评分后公开的`m`数组，不含随机目标生成seed。公开优化器seed属于算法配置，和私有目标seed不同。

新数据生成脚本需要同级冻结R2的`src/physics.py`；原仓库已有该文件，完整ZIP另包含最小同级依赖。它只复用物理实现，不读R2材料真值或旧实验结果。

更细的原始输入/输出合同、数值wall与子进程wall边界、候选去重、phase计费及全部命令见[EXPERIMENT_TECHNICAL_NOTES.md](EXPERIMENT_TECHNICAL_NOTES.md)。
