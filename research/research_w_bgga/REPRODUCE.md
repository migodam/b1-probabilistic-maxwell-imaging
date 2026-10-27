# 复现、数据恢复与计费

本目录在既有B1仓库research/research_w_bgga。Python3.13、CPU单线程、complex128/float64；依赖精确版本见requirements.txt，环境证据见source_review/PREFLIGHT.json。固定vendor来自指定旧提交，不从main替换科学来源。

## 恢复已保存证据，无新增PDE

下载[本轮Release](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-w-bgga-20261028)的G0、G1与收据档案，在本目录解压，内部路径以results/、figures/为根。Git保留compact结果、代码和图片；大型NPZ、逐原子JSON/动作/区间在Release。具体文件/字节清单见PUBLICATION_ASSETS.json，无文件摘要校验。

```
python -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python tools/run_job.py G1_LOGIC reproduce-tests .venv/bin/python -m unittest discover -s tests -v
.venv/bin/python tools/run_job.py G1_LOGIC reproduce-g1-plots .venv/bin/python tools/g1_figures.py
```

job label必须新建且唯一；旧结果和receipts不覆盖。65项测试无PDE，不能充当Maxwell资格。g1_figures重建自己的结果表和图，不打开GT/PDE；它读取本轮保存JSON。可直接审核G0的复场/复响应/差分/端点NPZ及完整结果表。

## 重新独立计算有限模型方法

旧bank在[W-SLAM Release](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-w-slam-20261025)的controlled_E2E3数据档案；恢复包含FINITE_GRID、FINITE_BANK、OBSERVATIONS的finite_bank/online目录。不能仅下载本轮摘要代替旧bank。

```
.venv/bin/python tools/run_job.py G1_LOGIC reproduce-seven-methods .venv/bin/python tools/g1_replay.py --old-bank-root /path/to/finite_bank/online --output results/g1_replay_repeat
.venv/bin/python tools/run_job.py G1_LOGIC reproduce-extra-scores .venv/bin/python tools/g1_offline_extra.py --old-bank-root /path/to/finite_bank/online --replay results/g1_replay_repeat --output results/g1_extra_repeat
```

原cache-map本机路径只作来源记录。可移植入口只接受显式old-bank-root下相同允许文件、prior ID和字段；不自动搜索main或locked/inference目录。offline truth只由extra scorer在全部在线输出seal之后读取，其路径是online同级offline/FINITE_PRIOR_TRUTHS_AND_NOISE.npz。没有数据则BLOCKED，不能以合成阵列替代。

旧已付费用不能免费转为新部署优势。逻辑重放仅记录cached queries，非新增PDE/RHS/FE速度证据。归一化工具可显式指定old-bank-root，默认只对原保存M4 trace做有限修订；原F64失败数组不改。

## 实际物理重跑

本campaign的物理阶段已封账，直接在原目录运行G0驱动会拒绝。完整物理复现应新建独立目录/scope，复制本目录的src、tools、tests、vendor与冻结配置；不要覆盖原PHYSICS_LEDGER或结果。新scope的机器资源由复现者明确登记AUTHORIZATION，G0上限4000RHS/500FE/7200process CPU秒，保留原科学门。新建空results；可复制不含物理调用的source provenance信息。

已执行的命令顺序：

```
python tools/run_job.py G0 fresh-qualification python tools/g0_qualify.py
python tools/run_job.py G0 fresh-derivatives python tools/g0_derivatives.py
python tools/run_job.py G0 fresh-operators python tools/g0_operators.py
python tools/run_job.py G0 fresh-stress python tools/g0_stress.py
python tools/run_job.py G0 fresh-decision python tools/g0_decision.py
```

先检查初始资格与实际预算再进入下一驱动。g0_stress内部声明了固定C4 closeout，达到90%仅完成这些已声明核验。无“all自动扩展”命令；超过硬上限拒绝，失败仍计费。g0_derivatives的合法单侧修复已进入源码。固定数组作图依赖results/preflight/FIXED_MESH.npz，可用零PDE mesh_benchmark生成，不能据此声称新增full-wave。

## 输出与限制

所有PDE调用先reservation、完成记结果；wrapper聚合process CPU只加一次，不再加kernel/子层费用。PHYSICS_LEDGER为物理/费用权威，COST_SUMMARY为快照；逐job receipt包含实际命令和失败。

PUBLICATION_READBACK仅验证可访问性，不升级科学门。网页生成用build_site，成本/表更新用finalize_results。post-publication accounting资产补入最终收据，明确其覆盖范围。无hash、无备份流程、无未授权NN/3D/新观测。
