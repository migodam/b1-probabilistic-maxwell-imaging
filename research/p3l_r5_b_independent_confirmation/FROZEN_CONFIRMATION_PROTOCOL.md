# R5-B 独立确认：冻结合同 v1

冻结日：2026-10-27（Asia/Singapore）。本轮为新的明确用户授权，原 R5、R4 的封账不变。起点为 R5 `ad71a86ef3c28bf48049199f733080a0f2b95e33`，实际历史输入资产以 BLAKE2b 清单标识；不做架构、训练、容量或 head 搜索。

## 问题和主总体

唯一主总体是两批全新、互不重叠的原 P1 prior-predictive parents，各256。沿用此前只生成过 seed 清单的 A/B，原整数 latent_seed / noise_seed 不重选。另8个独立成本探针不进入测试，也不读取其 GT 评分。512份测试材料与观测均为新增 Full Maxwell 计算。

原 P1：Gaussian64 × Uniform2、原相关场、原 C² 饱和 decoder、被动复材料、24×24 TM/VIE、8源×32复接收器、256 ROI。sigma=0.00119211798514519。每源32实部后32虚部，Phi=||F-y||²/sigma²。原先验和 acquisition 不变。生成过程采用原 MaterialPrior.draw(default_rng(latent_seed))；噪声使用独立 default_rng(noise_seed) 的 sigma/sqrt(2) 实虚正态。

## 冻结动作和推断

主要动作：physics-only Probit、冻结 R4 Strong CNN、优化后 original Gaussian K32 fusion、R5 deterministic Q64 + matched head、R5 W128/N32768 distillation + matched head。三个训练槽固定 2026102601/02/03；前两种历史 head 映射至 2026102501/02/03。物理基线无训练随机性，其损失在三槽中重复只为数组对齐。额外 deterministic + frozen old head 仅解释性，不计入主 family。

Canonical 为单线程 CPU、batch64；物理 float64/complex128，网络保持原 float32。Q64 是原最终推断默认，不用成本图中 Q16 的数值快路替换它。原3e-6跨CPU/MPS逐像素门的 FAIL 不覆盖；CPU重放资格单列。K32 encoder 只执行一次/批，保留逐 parent、四个chunk8的原 Gaussian 随机数合同。

## 原物理锚点

直接导入历史 R4 的 `reconstruct(full_TRF,max_evaluations=96,polish=False)` 和 `probit_optimize(cap=64)`，不是重新实现或降低迭代上限。前者从零场、属性中点开始，native objective 为 Phi+||z64||²/2；后者以 native 最优有限点初始化，优化原 probit-coordinate objective Phi+||x66||²/2。ftol=xtol=1e-9、gtol=1e-7、exact TRF、x_scale=jac 原样。原有边界初始化 inward clip、有限最佳点保留及停止规则原样。不把 cap 和 scipy success 混为驻点资格。

每个候选状态在计算前计费并立刻无损保存 q、材料、预测、backward residual；最终另存 current、field、J。仅不重复落盘大 LU，不改变 live-cache、矩阵、优化步或求解容差。worker串行，receiver cache有界，parent结束后释放。未收敛但有原函数返回的最佳有限端点照常进入总体；无端点或接口错误先停机诊断，不能看 GT 重启或以廉价锚点补数。

## GT 隔离和封存

数据生成、online anchors、模型预测、offline evaluator 为不同进程。online程序只打开观测、同parent的原物理context与冻结模型。全部512份主要预测及seed输出封存、内容身份登记后，才运行可打开 offline_only GT 的评分进程。A批评分不用于改变B批任何合同；实际安排统一延迟至两批预测全部完成。GT不参与初始化、排序、停止、候选选择或异常重启。

## 风险和统计

主风险：raw full-grid complex MSE，分母576。先每parent对三个冻结seed的loss取平均，再按parent平均；图像ensemble风险另列。complex MAE、ROI MSE、passive-clipped、各anchor相区和单seed均为次要描述。相区只使用anchor field/halfwidth的±1，不能按GT分类。报告相区比例、每单元风险与对full-grid增量的贡献。属性解释受相区污染限制，不将风险图称为 posterior 或合法代表。

主 family 4比较：det matched vsK32、distill vsK32、distill vsStrong、retention≥80%。 pooled512为主，A/B为独立复验描述而非挑选显著的子集。20,000次paired、按A/B分层parent bootstrap；同时重采样所有方法和指标。Bonferroni4，双侧98.75%区间；retention主检验为线性contrast `(Strong-Distill)-.8*(Strong-K32)` 的单侧98.75%下界>0，且K32增量分母稳定为正。比率保留全部共同重采样分子分母，仅分母正且调整区间排零时解释。详见冻结统计JSON，不以点估计>=80%替代统计门。

## 预算与费用

独立总硬上限3,300,000 complex RHS，包含真值生成、8个成本探针、Native/Probit、J/adjoint、失败与fallback。旧预算不复用。每forward8 RHS、一次receiver-J额外32 RHS；缓存动作不虚构新RHS，但实际CPU/IO另计。每个新parent开始前保证完整配置最坏上限可支付；费用不足只停止完整parent，不降级配置。现有历史费用约1.93M/512，期望本轮足够，但以实际账本为准。

至少30次成本重复，分别测CPU batch1/64、provider-only、完整学习链、cold load/首推、physics anchor及全系统。报告p10/p50/p90；原K32为优化实现，确定性成本使用相同Q64路径。无法访问MPS时明确NOT_AVAILABLE。CPU/GPU时间、历史teacher支出、训练和部署摊销分开，不把provider约39倍称整个Maxwell系统39倍。

## 继续与停止

完整512在预算足够时必须完成，不因初期负结果停止。模型和协议不根据新测试修改。真实bug保留原失败和全部受影响预测，对所有动作公平修复并降级必要资格。只作确认；后续任何新模型必须另起任务。所有统计结论以实际数据为准。
