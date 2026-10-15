# 阶段裁决

|阶段／主张|状态|证据与边界|
|---|---|---|
|E0 prior / acquisition / noise 可复原|PASS|64 TRAIN parents、512 材料条件原生成器重放误差0；零新增物理|
|旧 T/G 可直接作为原先验竞争见证|FAIL|有 full prediction 的历史端点无 PHASE 通过；原 pool 的 L 不是 full likelihood|
|原缓存存在结构不同且 d_pred≤2 的原先验 pair|NOT_FOUND|8,064 比较，1,632 结构合格；最小结构合格距离44.8267177|
|E1 有至少两组独立原先验竞争见证|FAIL|新搜索最小49.2236228；四 pair 两端各32噪声重放均0竞争|
|E2 共同 corrector 四臂比较|NOT_RUN|E1 前提未满足，不能追加预算制造前提|
|新的 N3 先验保持修正|NOT_RUN|仅保留原N3及历史已付费诊断，没有借E2额度运行|
|P0 内生成策略改善／Best-of-4增益|NOT_RUN|没有四臂或合格多解集合，不能评分或用off-prior填K|
|GLOBAL_UNIQUE|NOT_ESTABLISHED|有限搜索不能证明不存在竞争解释|

**最终状态：NO_GO_FOR_GENERATOR_COMPARISON；AMBIGUITY_UNRESOLVED。**

现有 T/G 候选为 PRIOR_CONFLICT 的历史诊断。新的离线 benchmark anchor 是 P0_EXACT，但不是在线反演器发现一个正确单解的成绩。E1 输出中的 ONE_SUPPORTED_CANDIDATE_AMBIGUITY_UNRESOLVED 仅为校准台存在一个生成端点的保守状态，不应读成原 development 观测已获得 P0_EXACT 在线重建。

下一阶段仅建议：有新授权后，优先做原相位/原生成器保持的 N3 最小必要修正，并研究表示误差和噪声一致性；没有证据扩大生成模型、噪声或自由补偿。当前不启动该建议。

45 states / 680 RHS 后主动停止。E1剩余19 states /1368 RHS、整个E2额度均封存，不借用历史余额。
