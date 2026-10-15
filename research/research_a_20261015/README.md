# Research A — Sequential Maxwell Material Discovery

2026-10-15 独立CPU研究阶段。研究对象是“同一观测下保留不同、合法且full-Maxwell相容的材料解释”，不将候选评分伪装成posterior probability。

[从这里开始](START_HERE.md) · [中文研究报告](RESEARCH_REPORT_ZH.md) · [Top-K图库](index.html) · [方法比较](METHOD_COMPARISON.md) · [消融](ALGORITHM_ABLATIONS.md) · [失败](FAILURE_CASES.md)

实际结果：Mie确认32对象、每对象3–10个静态archive世界，合计193；Best4比普通multi-start改善17.1558%，但在线Selected无改善，history adaptive未胜静态。独立同结构区间重启复现端点，收益主要来自结构分层。A30困难二维8条件中5个无合法解，新增结构/Best4收益0。完整状态和资源以本目录裁决、资源汇总和原账本为准。

## 代码与数据

- 仓库 `src/research_a/`：独立预算、冻结目标、受保护结构编辑、普通重启/population/adaptive、离线评分、A30反证与交接格式。
- `configs/protocol_v1.json`、`configs/CONFIRM_FREEZE.json`：原物理/噪声/表示及确认冻结规则；`EXECUTION_AUTHORIZATION.json`：本次执行overlay，不修改收到的Pro原配置。
- `results/CONFIRM_CAMPAIGN.json`：32个确认对象所有方法真实搜索轨迹、完整预测与修正记录，不含离线真值。
- `results/SUMMARY_E1.json`、`results/CONFIRM_METRICS.json`：方法汇总和逐对象离线指标。
- `results/E2/`：历史A30复杂材料下实际预测、审计、失败和独立stage费用；`RESULT_SUMMARY.json`给无解状态和全条件分母。
- `MaterialHypothesisPackage/confirm_00...confirm_31/`：32份实际静态候选的原生profile及独立common exchange；精确radial interval复ε数组唯一确定完整球对称三维材料。
- `contracts/`：原A schema/validator/test保留、no-digest及未暴露本构缺陷execution overlay、独立common v1和本地验证器说明。测试fixture只在临时目录，不是正式发现。
- `figures/`、`index.html`：统一材料色标的候选图、真值仅离线评分展示、散射拟合、结构差异及E2负结果。
- `COST_LEDGER.jsonl`：E1 append-only attempt/job账本；E2独立账本在其结果目录。不能重置run ID或重复累加32包global cost snapshot。

大数组如由主线程放入现有仓库Release，使用该次证据资产中的相同相对路径；不得将缺失缓存当成免费物理重新生成。公开页面与Release固定入口见文末。

## 复查已经完成的结果

以下从仓库根目录执行。使用已安装NumPy/SciPy/PyTorch的现有Python环境；本轮没有为了小测试另外安装工具。

```bash
export PYTHONPATH="$PWD/src"
python -m research_a.contracts research/research_a_20261015/MaterialHypothesisPackage/confirm_00/manifest.json
python -m research_a.contracts research/research_a_20261015/MaterialHypothesisPackage/confirm_00/common_exchange_v1.json
```

这两条只重放schema/语义、数组、白化预测和残差，不调用solver、不计算或核验摘要、不改变预算。它们不是连续物理认证。

本轮 `tools/analyze.py`、`tools/build_gallery.py`、`tools/export_actual_packages.py` 的实际输出均已保存。导出脚本在 `Budget.job('export-actual-packages')` 中运行过；随后decoder元数据修复在独立 `export-actual-packages-semantic-repair` job重建，均新物理0；命名job拒绝重复运行是为了保留执行费用。**不要把重新执行run/campaign作为读取结果的必要步骤**。需要另一个campaign时必须有新授权且另建protocol，保留旧账本。

## 单次新观测如何使用

实际API在 `src/research_a/target.py`、`search.py`：由已授权Budget构造Target，执行coarse reconstruction及选择的搜索arm；只有 `endpoint(...)` 的原full-Mie验算通过且material-space去重通过，才能进入archive。proxy搜索不导入真值/teacher/reference。已有bank与旧模型路径属于明确历史资产，需从证据包或原工作目录按provenance提供。

交接可使用 `make_mie_context(...)` 和 `write_packages(...)`；输入必须是真实已验证候选。详细字段见 [INPUT_GUIDE](contracts/INPUT_GUIDE.md)。`evaluate_new_action(...)` 要求单独budget token和实际记账的冻结物理callback；默认不执行任何新增测量。B原生格式兼容性未假定，也没有访问B工作区。

`weights`、`posterior_weight`、represented/unshown modeled mass均为null；未探索结构和模型外世界unquantified。可以比较候选材料差异或未来测量的区分力，不能进行带伪概率的posterior information gain。

## 物理和统计边界

Mie是三维vector的球对称分层模型，中心切面不是一般非同心三维成像。E2是二维scalar TM，不能混称vector三维证据。两平台噪声不同：Mie每实通道SD=σ；A30复噪声实虚SD=σ/√2，故白化sqrt(2)/σ。gain在Mie有界profile、在E2固定1，没有通过放宽噪声制造多解。

Best-of-K只离线评估；在线选择不读真值。确认对象来自冻结后同prior的32个新parent，A30的4parent×2条件仍是历史已暴露development。没有新大型训练、旧locked truth、旧P2预算借用、A17冻结稿修改或新非同心3D campaign。


公开入口：[可浏览材料图库](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_a_20261015/) · [本次证据Release](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-a-20261015)。所有历史阶段在同一B1仓库保留。

## 原执行命令与可移植路径

以下是本轮已付费执行顺序的记录，**不是要求在封口账本上重跑**：

```bash
# Python工作环境已有numpy/scipy；MDN提案还使用既有torch。
export PYTHONPATH="$PWD/src"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
# B1_P2_RESULTS可指向Release中的旧P2共用bank/资格metadata目录。
# C1_P2可指向旧C1的p2目录，其results/models中提供固定checkpoint。
python -m research_a.run --phase E0
python -m research_a.run --phase prepare
python -m research_a.run --phase audit
python -m research_a.run --phase dev1
python -m research_a.run --phase score-dev1
python -m research_a.run --phase dev2
python -m research_a.run --phase score-dev2
python -m research_a.run --phase dev3
python -m research_a.run --phase score-dev3
# CONFIRM_FREEZE.json在此处生成，先于truth。
python -m research_a.run --phase truth
python -m research_a.run --phase confirm
python -m research_a.run --phase score-confirm
python -m research_a.e2 --e1-gate research/research_a_20261015/results/E1_ENTRY_DECISION.json
```

必要的原始响应/真实材料与全部脚本见同一证据Release。A30通过`A30_ROOT`指向冻结A30代码与所列8个cache。实际重现求解属于新物理消费，需要新授权；`RUN_CLOSED.json`会拒绝在本轮预算上继续求解。已保存结果、接口和统计可零物理读取。

新测量的实际物理callback实现为 `research_a.handoff.make_mie_action_evaluator`：绑定准确径向edges、原世界共享gain和由调用方提供的独立授权enforcing ledger；导出包没有执行该callback。它不以当前A剩余额度为B自动充值。

每轮执行时Git HEAD与未提交代码状态如实记录在`ITERATION_SOURCE_PROVENANCE.json`；各轮未提前建立独立committed snapshot是provenance限制，不伪造已经冻结的历史commit。最终代码、原协议、原始轨迹和明确修复均随本次版本发布。
