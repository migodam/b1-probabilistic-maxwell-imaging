# 独立静态物理与费用审查

结论：指定II/III费用与账本逐项一致，未发现Full-forward重复收费、未消耗的forward credit或缺失的inverse setup。发现一处真实GPU增量序列化遗漏，父线程已修正新进程代码；历史累计cap一直正确。另有原始回执/附录/快照范围需要按下述方式读取。

本审查冻结在ledger sequence647516、917,318,397字节的完整提交前缀；审查95个CLOSED jobs。3个活跃job只记录名字及排除理由，不读其进行中结果或当成已闭账。无新物理、模型训练或推断启动，无GT值读取，无科学或novelty判定。源结果、授权与旧账本未改；只登记本审查自己的零物理analysis回执。

## 1. 指定RHS对账通过

| 项目 | 完整complex RHS | 核对方式 |
|---|---:|---|
| event reuse两job | 14,147 | 原full-grid8,003 + reduced补记6,144，后者只记一次 |
| decision star | 62,770 | inverse4,608 + Full audit8,192 + GMRES source8,192 +实际preconditioner41,778 |
| interval pair | 118,720 | Full65,536 + interval53,184 |
| III两轮 | 459,712 | Euler291,648 + Heun168,064 |

物理state/LU也一致：event pair768/768，star2048/1024（含1024 iterative material-state attempts），interval pair14264/14264，III16060/16060。reduced768个LU另列，不能伪装为576阶物理LU。

Star raw standalone GMRES汇总46,402漏了source-system8,192；`COST_ADDENDUM.json`已正确改为54,594。实际62,770早已收费，无需再加8,192。它是counterfactual standalone arm汇总修正，不能再次追加到actual费用。

Interval的53,184已包含8个original R0各576 RHS=4,608以及6,072次Full fallback×8=48,576；不是只报fallback而隐去warm setup。Full arm为8192×8=65,536。共同StageI anchor获取仍在原账本，不能将此warm对比当从零获取观测/context的cold费用。

### 历史event arm表格必须使用附录

原report的1859-vs4096只表示full-grid/setup列，不能作为当前all-complex RHS对比。第一轮Full arm=4096；reuse arm=1859 setup+4096 reduced=5955，总10051。第二轮pCN审计Full=2048、reuse reduced=2048，总4096。两轮合计14147。故第一轮reuse在此口径下RHS比Full多约45.4%，原“54.61%减少”只能保留为历史full-grid-only范围，不可复用为all-complex节省声明。当前cap对账已正确，原始历史文件保留，使用 `RESOURCE_UNIT_CLARIFICATION.json` 及本补充表即可，不能再收费。

## 2. 真正发现的GPU行序列化遗漏：已修正未来代码

审查启动时的Book `_update` 对 `after.stages[*].GPU_seconds` 和cap做了正确增加，却未在每行顶层写 `GPU_seconds`。28个历史事件的顶层delta求和会读为0，而完整after差分为4202.726937542秒，全部已进入cap。

父线程随后将新代码改为顶层显式 `GPU_seconds=gpu`。只读重建通过：`DERIVED_GPU_DELTAS.jsonl`逐事件列出历史差分；`GPU_RECONSTRUCTION.json`记录原cutoff和总额，没有修改旧行或再次导入GPU费用。旧源证据在本审查 `source_snapshot/core.py:142`，未来修正见当前 `src/p3r4/core.py:142`。

已运行进程仍可能持有旧imported Book，所以它在补丁时间之后关闭也可能缺该字段。读取器应按“字段是否存在”决定使用after差分，而非只按时间分界；不要将旧进程迟到的一行误报为新代码回归。

## 3. 回执只读汇总的一个范围例外

`external_neural_audit_profile_reconciliation/CLOSED.json`的.001335625秒是登记动作自身。它的独立导入事件sequence87105还有1.885279750秒，因此该job总归属wall=1.886615375秒。外部CPU1.116593秒和RSS538,624,000 bytes也在事件metadata，原CLOSED profile只有登记CPU/RSS。

这不是未收费：ledger after正确包含导入工作。需要event-aware读取，而不是只sum CLOSED字段。父线程已在 `GPU_LEDGER_SERIALIZATION_CLARIFICATION.md`说明范围；原回执不能改写，也不能把外部GPU/CPU再次导入。

本审查自身也分开登记主扫描和后续只读处理：主扫描约17.77秒，后续基于已观察CLI耗时保守记20秒（仅本审查CPU墙钟、GPU0），精确范围见 `AUDIT_OWN_COST_SCOPE.json`。不把登记回执自身冒充导入工作总时长，不编造未测CPU core time。

## 4. physical wall与controlled wall不同

| 已闭arm | result记录wall | CLOSED controlled wall | physical-operation wall |
|---|---:|---:|---:|
| interval Full | 229.782692s | 229.817106s | 58.898236s |
| interval evaluator | 355.794177s | 355.824848s | 45.774842s |

结果字段在close前写出，CLOSED还含最后登记/持久化尾部。physical只覆盖被显式包裹的物理解算操作，不含全部cached矩阵作用、界计算、循环和持久化。完整端到端比较应使用CLOSED controlled或明确范围的总wall；不能把45.8-vs58.9秒的physical子段换称整条链耗时。所有95个已闭job均未出现physical wall大于controlled wall的异常。

## 5. source_snapshot范围不能倒推补齐

95个闭账job中29个没有SOURCE_REVISION；41个旧job虽保存p3r4模块，但没有其Python entrypoint；8个是非文件入口/其他窄范围，17个已保存entrypoint。详见 `SOURCE_SNAPSHOT_INVENTORY.json`。

当前Book已增加tools快照，注释明确older jobs仍是narrower snapshots。这个未来改进不能证明旧runner的未提交工作树版本。发布时保留逐job覆盖情况；若事后归档现有源码，应标记post hoc，不伪装成运行前快照。Git revision、argv与模块快照各有范围，不能互相替代。没有为此重跑或修改历史。

## 6. 模型数据与物理owner

- 恢复的32768 TRAIN、2048 validation joint与512/128 anchor contexts来自已付R2路径，`DATA_RESTORE.json`和执行argv有精确来源；恢复本身0新Maxwell RHS，不代表原获取免费。cold总费未在本次重建，必须保留旧owner。
- 已闭`risk_anchorbin_phase16_s2026102501`等argv确实使用旧512/128 anchor文件。新1024/256 context计划和packing contract尚不能替代其as-run数据资格；cutoff时validation contexts仍有活跃job，本审查不将其称完整。
- Operator pack的10351个paid J记录、10063去重配对来自III459712 RHS，J-state acquisition414040与其他Full端点45672一致。50,315个real≤66 GN参数solve不属于complex Maxwell RHS。新T1全部MISSING的限制已经标明，没有额外inverse/adjoint查询。
- CNN_CLASS_EQUIVALENCE metadata记录8对象max差1.20137e−7、0物理RHS/0优化updates及旧residual-scale参数映射。这里记录其来源与费用范围，不作最终架构/novelty判断。

## 7. Representative静态费用控制

当前仍是准备代码，未发现已闭真实representative job。静态接口与mock检查通过：每parent/arm16states/256complex RHS；source路径每trial forward8+gradient8；receiver-J路径forward8+J32。两次可行性restoration已预留在同cap内，所有failed/trial费用通过原Book预约保留；KKT使用最终缓存gradient，不新增越cap查询。raw/clipped mu只读取冻结online prediction keys。

final严格检查原Phi<=282.873792，选择原合格anchor、接受且可行的SLSQP点或新Full验证可行的restoration点。KKT只是数值一阶诊断。此处没有真实Maxwell representative验证，不能用4项decoder/mock PASS代替Output3物理资格。

## 机器入口与行动

`FINDINGS.json`为最终分类；`LEDGER_RECONCILIATION.json`保留原始扫描差异，28项GPU after mismatch已由只读重建解释，不是28个cap错误。`DERIVED_CLOSED_JOB_COSTS.json`提供scope-aware GPU成本视图。所有指定RHS合计、forward credit、Full事件数量和solver counter核对通过。

剩余行动是范围正确的汇总/发布：GPU缺字段回退、reconciliation事件只计一次、旧full-grid表格附上all-complex口径、旧source覆盖和数据owner如实列出。未发现需要再次补扣物理RHS的项目；科学支持、效率与最终claims仍由父线程判断。
