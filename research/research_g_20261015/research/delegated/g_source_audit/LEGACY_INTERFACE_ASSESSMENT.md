# 完全相同历史冻结后端的有限接入核对

这是主Codex提出历史后端例外后的只读技术评估，不是T交付，也不是新物理执行记录。用户明确允许使用“T冻结physical handoff，或完全相同、已经核验的历史冻结后端”。本审计不增加corrector实现。

固定A30 commit `717d3d26a290248bd2b9aeb55f3d03694b407b84`及Research A commit `f1b674eb4893e39126701ed48f3d3dcef3098c3c`中的`src/research_a/e2.py`可明确定位既有物理adapter与普通材料GN函数。

`_gn_proposal(model,state,basis,observation,white,trust=.35,damping=1e-3)`使用当前材料field及receiver transfer形成material Jacobian。QP的actual material box为Reχ[-.5,4]/Imχ[0,.5]；trust=.35；SLSQP最多60步，ftol=1e-10；最后仅做合法材料射线收缩`_legal_step`，不是新增T的proposal-retentive corrector。`_bases(representations)[1]`是16个real Haar方向加独立imag方向，共32实坐标，source/receiver packing与A30一致。

需保持以下边界：

1. frozen field support仍是576格完整χ，但32坐标corrector不等于一般全field恢复能力；改变proposal搜索策略不应掩盖这一限制。
2. `_gn_proposal`会构建receiver transfer，实际32个adjoint RHS应计费；不能仅以8源full-forward RHS统计全部成本。每步新材料的receiver adjoint不可错误复用为旧state缓存。
3. 独立endpoint full validation和首次base重验均需在G账本支付；加载T/历史缓存只读复用不构成G物理加速胜利。
4. 原A E2比较material metric为RMS|Δχ|/1.2；G本轮固定dM采用real/imag物理范围归一化，两者不是相同去重计数。两组原始结果须并列保留，不能把新指标重算冒称旧实验成绩。
5. T-guided compensation-feasibility arm依然是WAITING_FOR_T，不能把普通GN/weak方向重新命名为T机制。

在任何新物理前，已登记并执行零物理saved-array核对（`LEGACY_BASE_AUDIT.json`，wall=.041423125s、CPU=.040732s，BLAS1）。选取规则是历史观测数据相容性，不读truth：

|parent/condition|原缓存σ|按原packing/white重算base residual|material box|原observation/background逐值一致|
|---|---:|---:|---|---|
|p250069 weak|.001270913712869431|20.929164702218827|合法|是|
|p250073 weak|.00119211798514519|20.891919121358494|合法|是|
|p250077 weak|.0011846108030513216|22.994967915752227|合法|是|

共同阈值23.785449004057032；field为576，prediction32×8。只白名单读取旧online数组`__base/__base_prediction/__observation/__background`及原pool的observation/sigma/background，没有打开legacy offline数组或真值。旧保存数组不含每base当前内部constitutive/solver defect，应由主Codex的新独立base full-validation回执补充，而不是填0。

因此静态source/packing/observed-target核对没有发现该legacy接入本身的冲突。是否满足正式E1条件仍由主Codex在获得独立full-base验证后裁决；小样本只能称3个已暴露weak parents的限定legacy pilot，不是4-parent/五臂完整T比较，也不是confirmation。
