# R6主D：紧凑数值交接

这两份NPZ只供下一轮GPT Pro离线阅读和再分析。DLOW与DHIGH各含256个不同parent，分别对应原sigma的0.5倍与2倍；两组不合并估计。本轮全部512个主D真值已经由根线程在全预测封存之后解锁，因此这个包**已经暴露GT**，不能再称作未见测试数据。它不含补充128例或CDEV。取样顺序来自冻结parent清单，包含全部parent和全部三seed，没有按真值挑选最好样本或seed。

入口是[数值清单](COMPACT_NUMERICAL_MANIFEST.json)、[DLOW无损包](compact/DLOW_COMPACT_NUMERICAL.npz)、[DHIGH无损包](compact/DHIGH_COMPACT_NUMERICAL.npz)。NPZ用`numpy.load(path, allow_pickle=False)`打开。`parent`给第一轴顺序，`seeds`给第二轴顺序，`action_names`列出六动作。每个`prediction_动作名`形状为`[256,3,576]`，是冻结材料风险动作的输出，**不是posterior样本**。三seed不是独立真值复制，六动作也不是可按GT选择的候选posterior。

`GT_material`保留原complex128材料对比度chi，576格对应24×24网格；`GT_native_q`保留66维原生先验坐标，前64维为Gaussian场，后2维为Uniform属性。`full_physics_native_q/material`与`low_physics_native_q/material`是Native96→Probit64、Native24→Probit64后的最终物理端点，原件键为`q0/chi0`。额外的`full/low_native_stage_q/material`是后续Probit之前的Native阶段端点。`*_gaussianized_x`单列保留Gaussianized坐标，不能与原生q混用。这些都是点估计或已生成真值，也不是posterior抽样。

`observation`保留32×8 complex128原始观测，`sigma`逐parent保留原数值；`GT_noiseless_response`、`noise_innovation`保留生成时已经存下的无噪响应及标准复创新。本次只核对`mean+sigma*innovation`与原观测逐字节相同，没有重跑Maxwell。`full/low_physics_response`是原件缓存的最终物理响应。`full_anchor_phase_code`来自Full最终anchor的decoder而非GT：-1为ROI外，0为host，1为transition，2为pure；`full_anchor_phase_fractions`三列依次host/transition/pure。

清单逐项记录原始路径、大小、已有封存身份及来源数组映射，并记录新GT原件身份。所有收录数组经过dtype、shape、C-order字节一致核对；每组5,397个源数组/切片全部通过。原件未改；新增物理RHS、训练和推断均为0。`build_compact.py --verify-only`可重做相同离线核对，不覆盖已完成文件。紧凑包没有J/H、current/field、全优化trace或历史训练权重；完整原始证据另见Release manifest，阅读包不能冒充全部证据。
