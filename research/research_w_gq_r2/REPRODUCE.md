**发布说明（2026-10-29）：** 此为用户授权后的公开版本。[网页](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_w_gq_r2/) · [完整原始证据 Release](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-w-gq-r2-20261029) · [匿名 raw](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/research_w_gq_r2/REPORT_ZH.md)。执行结束时的本地未发布说明及打包回执保留为历史状态；当前发布记录单列在 `publication/`，科研数组、判决和费用不变。

# 复现与文件入口

## 已保存证据回放：零新增 PDE

在本目录运行，使用 `ENVIRONMENT.json` 所记的 Python 或安装 `requirements.txt` 的同版本科学环境：

```sh
python tools/run_job.py E3 saved-replay-new-label python tools/verify_saved.py
```

`run_job.py` 会自动设置单线程并加入本包的 `reuse/fixed_vendor`、`reuse/previous_gq`、`src`。job 名必须唯一；已有回执不能覆盖。库路径优先采用本包随附源码，旧目录只作原工作区 fallback。

这会重算已保存物理 pair 的距离和材料 RMS、检查全部必要文件与图、解析源码并核对 ledger cap。它不是新的 Maxwell 资格，也不计算任何哈希。

## 实际执行命令与顺序

全部精确命令、stdout、CPU/RSS/exit code 都在 `results/jobs/<name>/{RECEIPT.json,OUTPUT.log}`。关键脚本：

|顺序|脚本|物理 scope|
|---|---|---|
|0|`freeze.py`、`preflight.py`|合同和最高分辨率成本；preflight 的一条真实 response 被复用|
|1|`microstructure.py`|两parent×三尺度×两纹理，固定目标面板|
|2|`microstructure.py --validate`|最差probe自己的独立细化|
|3|`symmetry.py`|镜像/非对称/仪器列；空气合法性修复保留|
|4|`imaging.py bank --phase E1`|f₀公共宏观模型|
|5|`imaging.py bank --phase E2`|新增频率与micro已知墙诊断表|
|6|`imaging.py all`|generate→online封存→score，零新增PDE|
|7|`continuous_fwi.py`|离网真值生成、相同两起点有限FWI|
|8|`continuous_fwi_third.py`|预声明一次方向补充，同一ROI中心起点给三臂|
|9|`compensation.py`、`--validate`|双向Jacobian、A0有限曲线和端点|
|10|`compensation_extended.py`、`validate_extended.py`|新增观测下重新优化，不只比较固定pair|
|11|`compensation_last_two.py`、`final_pair_refinement.py`|限定2步后、新pair三层一致性|
|12|`return_operator.py --execute-frozen`|原返回表示不足，保留|
|13|`return_transpose_repair.py --execute-frozen`|共轭输入/转置输出配对，只修正一个因素|
|14|`e3_logic.py`|790逐原子受控逻辑；封存后全表评分|
|15|`e2_micro_reference.py`、`--consensus-only`|在线封存后的4×12离线reference和共识诊断|
|16|`theory_checks.py`、`build_evidence.py`、`plot_results.py`|必要新公式和保存数据整理|

原始错误及修复代码/receipt都保留；完整顺序应以作业时间戳和账本为准，不能将表内顺序当并发时序。

## 重新运行全波

当前 campaign 已封账，`PHYSICS_CLOSED.json` 会拒绝任何新的物理预约。这不是数据损坏。要做一次全新复现，应建立**新的输出 scope**，只复制源码、reuse、原始输入和冻结配置，建立该 scope 自己的预算与空账本；不要删除本目录历史或把它的已付支出归零。原阶段 cap 同样适用。

在新scope中按上表运行，默认 `return_operator.py` / `return_transpose_repair.py` 是 PLAN_ONLY，必须明确加 `--execute-frozen`。不要直接逐脚本重跑现有目录，因为文件封存/唯一job标签会有意阻止覆盖。

本包含旧 finite exposed bank 的复现数据，位于 `inputs/old_finite`；其原2950FE/11800RHS/9288秒构建成本是历史成本，回放不得当作新PDE部署节省。所有locked数据均未包含。

## 原始结果定位

- E1：`results/micro/*.npz` 为各wall/target/frequency真实复数据，部分含P3场与网格；`micro_validation`为独立层。
- 在线：`results/imaging` 中测量、在线封存、finite候选响应与逐条件表；truth仅在独立评分入口使用。
- 连续：两个`continuous_fwi*`目录有每次实际call、梯度、材料和预测。最终测量规则选解在`results/CONTINUOUS_COMBINED_SELECTED.json`。
- E2：`compensation*`、`final_pair_refinement`保存每个合法端点/试步/多频观测；`BIDIRECTIONAL_SCHUR_SUPPLEMENT.npz`保存两方向Schur。
- 后验：`e2_micro_reference`保存完整joint、条件和像素方差、parent配对及2000bootstrap。
- E3：两个`e3_return*`目录有原始actions、sparse matrices、exactself模板、逐rank current/field/response；`e3_logic`有受控查询和封存后的reference。
- 根目录三个`E*_RESULTS.csv`是索引汇总，异构行由experiment/status区分；原始表不可丢弃。
- `QUERY_LEDGER.jsonl`是24×790条主端口逻辑查询；与FE/RHS物理账分开。

## 打包身份

`Research_W_GQ_R2_Results.zip`包含主要源码、复现依赖、附件理论、必要历史组件、实际结果、表图和回执。打包完成后的最终费用保存在相邻`PHYSICS_COST.json`，同时作为`delivery/FINAL_ACCOUNTING.json`追加到ZIP；原报告中的费用段不需要假装预知打包过程自身费用。

本轮本地交付，未创建新commit或自动发布。固定历史链接与本轮新实验身份见START_HERE.md。
