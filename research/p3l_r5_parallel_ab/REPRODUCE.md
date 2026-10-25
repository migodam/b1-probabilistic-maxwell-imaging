# 复现与原始证据

## 已封账，不隐式重跑物理

科学起点为B1提交85b7bcbf39f0eeab7062ecce17531538151082bf。本轮新增源码、固定配置、逐case结果、失败及checkpoint独立保存。`shared/cost_ledger/CLOSURE.json` 已关闭A/B付费入口；旧run ID和恢复逻辑不能重置余额。原始积分/优化脚本保留，用于审查和未来独立授权复现，不能因下载仓库就自动重新支付历史计算。

本机运行环境为Apple M4、24GiB、Python3.13.13、NumPy2.5.3、SciPy1.18.1、Torch2.14。CPU物理单worker，BLAS单线程；训练使用已有MPS。没有安装CUDA或使用外部云端。

## 路径与依赖

保留Gaussian下的原布局：B1/b1-p3l-r5-parallel-ab、只读B1/b1-p3l-r4-conditional-transport、B1/b1-p3l-r2-risk-posterior和A30/a30-maxwell-material-manifold。原prior源码/配置在本仓库src/p1_material及research/p1_material_target_20261018；A30使用已冻结vendor及原weak背景数组。

[具体路径与尺寸](audit/EVIDENCE_INPUT_REFERENCES.json) · [目录布局说明](audit/REPRO_PATH_LAYOUT.md) · [原始资产清单](audit/HISTORICAL_ASSET_AUDIT.json)。这些引用是实际文件，不是声称精简Pro ZIP含全部权重。没有复制175GB历史R4大包。

## A源码与结果

- `stage_a_response_coordinates/code/math_units.py`：解析数学测试。
- `code/run_maxwell.py`与`v1_failed_preserved/run_maxwell_v1_as_run.py`：初版与失败留痕。
- `code/run_v2.py`：区间求根、真实局部QR-RTO。
- `code/run_v3_strata.py`：rank0/1/2分层、原坐标/响应积分。
- 其他actual-run脚本及固定参数见[A复现](stage_a_response_coordinates/REPRODUCE.md)。
- `audit/curvature_bounds.py`与`audit/seam_maxwell_probe.py`：父线程计费的范数半径/Full-Maxwell接缝测试；不得在本封账目录重跑。

已保存18份观测、材料与预测、数值节点、权重、局部RTO尝试、原始和修复前后结果。任何NOT_RUN不通过补造文件假装完成。

## B推理与数据身份

主入口是 `stage_b_material_risk/inference.py:Action`。finalists为`deterministic_matched`与`distilled_w128_n32768`，各seed2026102601/02/03。`y`为[N,512]实化原测量，online仅包含相同strong anchor的observation/residual/q0/chi0/J。函数本身不读GT、不进行Maxwell求解、不在线更新。

```python
from inference import Action
action = Action('distilled_w128_n32768', 2026102601, device='cpu')
material = action.predict(y, online, batch_size=64)
```

模型及路径按[冻结选择](stage_b_material_risk/SELECTION_FROZEN.json)和[B复现说明](stage_b_material_risk/REPRODUCE.md)恢复。在线compact部署只需该seed的PROVIDER.pt/HEAD.pt及冻结的特征源码；不运行原Gaussian encoder或K32，也不需要teacher cache。但原teacher的离线训练与数据不是免费成本。

原CPU/MPS3e-6容差FAIL保留；同CPU两次独立加载的精确重复PASS是另一项检查，不能覆盖原FAIL。旧VAL重放只验证推理/保存一致性，不是新盲测。

## 原始费用与精度

物理权威为`shared/cost_ledger/LEDGER.jsonl`；B训练权威为`stage_b_material_risk/COST_LEDGER.json`及实际作业收据。根`RESOURCE_LEDGER.jsonl`仅索引/汇总，禁止双重计费。Git保留源码、紧凑指标与图；原始数组、checkpoint和大日志在本地及本轮独立证据包。交付/打包另计，不能写进物理算法节省。

配置身份仅使用canonical JSON的BLAKE2b标识；没有SHA-256检查或额外物理smoke campaign。无损NPZ保留原dtype，不降精度、不删除失败字段。图库为已有VAL固定前4个parent，不是新确认或best-of-K挑图。
