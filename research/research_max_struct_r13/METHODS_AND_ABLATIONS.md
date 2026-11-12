# R13 方法、对照与归因合同

本文件说明实际代码和冻结实验；恢复结论以封存后评分及中文报告为准。方法运行正常结束、完成原生迭代日程和材料恢复，是三件不同的事。

## 双场与事件

沿用 R12 的 `X=J/s`、`s=6`，状态范数为 `||X||F²+||Q||F²`。传播—数据集合采用完整背景 Green 与公开数据球；材料集合在每个单元合并所有照明，三维还合并三个场分量，共享一个实材料。反射参数 `omega=1`，无学习的投影或目标特定 mask。固定传播投影的 Cholesky、端口 Gram 分解以及每步矩阵动作均收费。

普通 DR：`c=PC(z); a=PA(2c-z); z+=a-c`。AP 则依次投影而不保留反射记忆。每份材料候选都需要另解该材料的原 Maxwell，不能用双场 gap 代替物理残差。

事件首先从一个公开状态 `a` 计算 `c=PC(a)`、`d=a-c` 和 `e=PA(c)-a`。精确最近点对的整数快进仅省略同一 DR 轨迹中的重复投影；近似对上的动作明确称有限状态 proposal。事件角度定律和适用条件见 `THEORY_DERIVATIONS.md`。

实际冻结候选包括：AP128 后整数/分数事件；同范数随机动作；法向 `t=2`；AP→DR；depth-one Anderson；事件后物理记忆重置。另一个开发得到的候选是 DR128→一次 A-shadow→首次分数事件→DR。它不需要 AP3000，且不被冒称精确 AP 快进。所有选择只使用公开状态。

在线事件额外三臂为相同准备下 `t=2`、无事件 A-shadow 重置、事件后物理记忆重置。它们在任何新材料评分之前登记，使用同一16对象和共同审计合同。公开旧 B 对象的区间中点反例被保留；没有据此加入中点或过冲参数扫描。

随机动作是单个冻结种子13013的对照，不代表所有随机多初值策略。核心的`random_integer`按整数事件的相同位移范数匹配；旧开发v1的`random_event`匹配分数事件，不能据它单独把整数事件差异归因于方向。倍增/随机控制仍沿用同一个可分辨gap及有限事件存在性检查；`t=2`只不使用阈值数值，并非无条件施加任何扰动。

## CIE、A+、G19 与完整 current

标量24²场景继承成熟 sourcewise CIE/Dual-CIE、global-PR Dual-CIE、R7 A+ 与 R11 G19。原生 Fourier masks 为 mf3/6/10/12，迭代数500/500/500/1500；共同材料576实自由度及相同 box/prior。

- CIE：sourcewise PR+，材料与 current beta 同为1/.2/.2/.2。
- Dual-CIE：sourcewise PR+，current beta 固定1，材料 beta 同上。
- global-PR Dual-CIE：上述双 beta，跨全部 source 使用一个 PR 步。历史标签 `global_dual` 专指它。
- A+：原 global current 更新，加 one-hop 材料残差与 Möbius/MM 条件读出。不是新物理信息，也不是仅改材料坐标。
- G19：原 global current 更新，加固定公开 Green 强子空间分母、完整分母导数与原材料内循环。rank由固定截断规则决定，名称不等于任何问题都严格 rank19。
- G19 full APCS：从起步即完整 Fourier APCS，DPCS 仍硬锁；beta、PR、材料反馈及先验不变。
- G19 full current：另外解除 DPCS 硬锁。其目标是表示消融，不是首次提出 NFFT-SOM。

mf12 在24²本来已覆盖完整 APCS。因此完整空间从起步的比较同时改变早期 continuation 和每步成本；不能把结果归结为末期容量增加。G19 继承物理 current 的合法 LS warm start；CIE/A+ 继承 cold coefficient start。共同材料起点不意味着内部辅助状态相同，转换缺陷单列。

`nfft_equations` 是原文 Eq3–6 的完整 current、实材料/box/prior、joint L-BFGS-B 适配，明确标 `NOT_ORIGINAL_NFFT_SOM`。它不能充当已完整复现成熟原 NFFT-SOM 的证据。

## 强原物理基线

direct L-BFGS-B 使用整段 CPU 预算，目标是同一个原 `Phi=D+R`。Full-J TRF 复用 FBG 的解析完整 Jacobian、同一 prior residual 和 `x_scale=jac`，不是弱化的一阶替代品。

**physical GlobalDual (beta=0)** 与 global-PR Dual-CIE 分开命名。前者每步解原 Maxwell，使用随材料变化的 total-field source Gram 和完整 Gram 导数；所有付费试探点最终仍按原 Phi 选候选。本轮只在兼容的标量共同端口场景运行，未冒称已推广到 source-specific Fresnel 或 vector DDA。

## NN、选择器与 corrector

只复用 R8 的24,737参数冻结网络及权重，无本轮训练。固定16对象按原公开输入产生 NN 材料，再分别送入 DR、G19、direct；不因 NN 初始 Phi 高而过滤。材料标签不进入网络输入或在线选择。旧训练总成本未知，不能给出训练摊销后的完整净加速。

所有方法主结果都是合法、预算内、已完成原物理审计的最小 Phi。secondary 为已付费球内材料中最小原 prior R，无合格材料则明确回退 primary。完整分母保留。

raw、corrected、primary 逐份保存。额外公开保留规则有两种：raw已球内则保留raw；以及已付费球内候选中离raw最小改变。后者在新 GT 打开前计算并封存索引。这些规则改变了材料选择偏好，不是同 Phi 搜索更强；原 corrector 全部成本照收，没有实测省去 corrector 的部署计时。

## 实际预算与审计

i9 主要比较为120秒整进程 CPU，proposal 的绝对 deadline72秒（含导入/setup），共同 corrector 最多48秒；direct/TRF 使用全段。付费前缀15/30/60/90/120秒。每个2秒 CPU bucket 后第一个完整材料更新允许物理审计，原生初始/阶段/最终审计仍保留并收费。

共同时间规则不强制造成相同候选数。一次不可中断线性代数可越过 deadline；该操作仍收费，但预算外候选不进入预算内表。

另有 M4 固定16对象360秒参考：CIE、Dual-CIE、G19、DR、direct、physical GlobalDual，以及单独 Full-J TRF。其硬件、NumPy/SciPy 与 i9 不同，分别比较，不能将跨机差异解释成纯预算作用。

条件救援对全48个 G19 primary 统一按原公开 residual/delta>1 触发；9对象触发，39对象原样保留。每个触发对象从相同父材料运行五臂。最终部署费用为全部48父 G19 加该策略实际救援，不能只报告触发子集或把历史状态当免费初值。

## 证据保存的限度

候选材料池和代表性完整物理场保存；部分正常结束的 DR/NFFT 辅助端点也保存。CIE/G19 并未保存全部 native PR 梯度、方向、stage state 与每个候选的辅助 current。不能据材料快照声称恢复了原优化器无缝历史，或给出所有辅助轨迹的因果证明。详见 `delegated/full_current/AUXILIARY_STATE_COVERAGE.md`。
