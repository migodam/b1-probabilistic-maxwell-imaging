# B1 × C1 Joint P2：零新增物理概率审计

日期：2026-10-13。七组共同概率检查实际执行完成，全部 PASS；warnings、errors 和非预期 physics/ledger imports 均为空。新增材料–频率评估 **0**。这里的 PASS 仅指有限概率恒等式、接口测度及负例检测成立，不是 Maxwell 精度、posterior mixing、现实 calibration 或后续科学裁决。

实际运行由共享 meter 的 `p2-e0-controls` 作业完成。脚本 body wall **0.348933 s**、CPU **0.343331 s**；共享作业 wall **1.315655 s**、child CPU **0.869773 s**，后者包括解释器及依赖加载。原始 stdout、stderr、child CPU receipt 和共享 job ledger 保留在 `results/jobs/p2-e0-controls/` 与 `results/JOB_LEDGER.jsonl`。本文件根据这次实际运行写成，没有重复运行概率控制或历史物理实验。

| 必要核对 | 实际结果 | 数值证据 |
|---|---|---|
| 共同材料 prior、decoder 与事件划分 | PASS | 与独立解析 prior 最大质量差 2.22e−16；分段 coefficient 积分差 8.33e−17；decoder 几何差为 0 |
| Proper prior-only posterior | PASS | 合成常数 likelihood 下事件 posterior 回到 prior，最大差 1.39e−16；实心原子质量 0.25 保留 |
| Truncated gain 边缘化与共享两频 gain | PASS | 独立一维 quadrature / dense Cholesky loglik 最大差 3.55e−15；条件 gain mean/variance 差 4.44e−16；known-gain 差为 0 |
| 同测度坐标变换与 prior sensitivity | PASS | coefficient 与正确 radius-Jacobian 计算的事件质量均为 0.4876646253；改成另一种 radius-uniform prior 后为 0.7506984792 |
| 窄/宽双峰连续 CDF、moments、quantiles 一致性 | PASS | 独立分区积分与解析 CDF 最大差 9.99e−16；归一化、均值、SD、quantiles 与 event mass 属于同一密度 |
| 复用高 ESS 漏模式与历史 coherence 负例 | PASS（负例被识别） | 最大 ESS 仍漏掉约 0.500025 的真实质量；历史方差违反分位数必要下界 98.638 倍 |
| Top-K、modeled remainder 与模型外语义 | PASS | 三个展示事件保留质量 0.6476844490；未展示已建模质量 0.3523155510；模型外质量为 null/unquantified |

## 连续生成 prior 与几何原子

两个 bump 支持互不重叠。对中心系数或环形系数，`z∈[-1,1]` 产生相同的无空腔材料部分，且原 coefficient-uniform prior 给它 **0.5** 的概率。两个部分独立，所以实心材料 prior mass 为 **0.25**；它不是额外给定的 topology 专家权重。

正几何尺度 (g) 与对应系数的关系为

\[
z=-[1-(g/a)^2]^{-2},\qquad
p_g(g)=\frac{g}{a^2}[1-(g/a)^2]^{-3},
\]

中心 (a=.45)，环形半宽 (a=.18)。该连续密度之外还必须保留零几何的 0.5 原子。积分 inactive coefficient fibers 是原 prior 的一部分，不能把它们丢掉，也不能将 coefficient-uniform 偷换成 radius-uniform。

事件 ID 按 `central_bin * 3 + annular_bin` 排列，中心 bins 为 0、(0,.10]、(.10,.20]、(.20,rmax]，环形 bins 为 0、(0,.06]、(.06,hmax]。独立解析 prior masses 为：

| 中心 bin | 环形 0 | 环形 (0,.06] | 环形 (.06,hmax] |
|---|---:|---:|---:|
| 0 | 0.2500000000 | 0.0332031250 | 0.2167968750 |
| (0,.10] | 0.0133243380 | 0.0017696386 | 0.0115546994 |
| (.10,.20] | 0.0557880880 | 0.0074093554 | 0.0483787326 |
| (.20,rmax] | 0.1808875740 | 0.0240241309 | 0.1568634430 |

事件互不重叠并覆盖声明的材料支持。边界 probes 包含 inactive 临界 `z=-1`、实心和两个最大几何，不读 truth，不预先输入 topology label。

## 噪声、gain 和 nuisance target

共同 likelihood 使用 Re 后接 Im 的实通道；各通道噪声方差为 σ²，而不是悄悄切换到另一种 complex variance 约定。gain prior 为 proper truncated (N(1,.15^2))，范围 [.7,1.3]，归一化常数为 0.9544997361。

独立核对包括两种方法：直接对该有界 gain 做一维 scaled quadrature；以及对 untruncated marginal 用 dense Cholesky 计算协方差与 logdet，再乘 posterior/prior 的截断区间概率比。五个固定 synthetic cases 使用边界与内部 gains。`conditional_gain` 返回截断之前 Gaussian 的 mean/variance，截断修正由 marginal likelihood 保留。

两频先合并通道，再对**同一个 gain** 边缘化。合成负控制得到 shared-gain loglik −5.6148552817，而分别边缘化两频的 loglik 和为 −5.6868384680，差 0.0719831863。这不表示哪一个更“好”；后者定义了两个独立 gain，是不同的概率 target，不能混用。

在同一个 toy likelihood `1+3*(r/.45)^2` 下，正确 coordinate transform 给 `r≤.20` 的 posterior mass 0.4876646253。保留原 0.5 原子而将 active radius 改成均匀分布，会变成 0.7506984792，差 **0.2630338539**。这是真实 prior 敏感性，不是坐标积分误差。

## 连续摘要与负例

窄/宽双峰例采用明确归一化密度：支持 [0,1] 上，0.55 的截断 (N(.22,.008^2)) 加 0.45 的截断 (N(.78,.045^2))。解析 CDF 与独立分区积分给：

- event `x≤.5` mass：0.5500000001103；
- mean：0.4719999479，variance：0.0785624064，SD：0.2802898614；
- q05/q50/q95：0.2093185781 / 0.2306814219 / 0.8349287086；
- 分位数必要 variance 下界：0.0100359394；实际 variance 满足它。

这次没有混合 node moments 与 prior-bin quantiles。另一个相关系数 0.995 的连续 Gaussian whitening 代数检查保持 covariance 和 Jacobian；它是一条 connected ridge，不被解释为两个分离模式。这些 toy examples 不代表 C1 已学会此类 posterior，密度训练与共同比较另列。

直接读取既有、已暴露的 B1 E0 `CONTROLS.json` 复用两项负证据，没有重建缺失 oracle、没有重复 physics：

1. 1,024 个均匀 proposal 点全部未命中 prior mass 0.0001、likelihood 10,000 的 rare region，ESS=1,024，却漏掉真实 posterior mass **0.5000250013**。高 ESS 只诊断已访问支持上的权重，不保证 mode coverage。
2. distinguishable/H0/Re(ε) 原生 mean、q05、q95 要求 variance 至少 0.0003185110915，而历史报告为 0.0000032290906。必要下界违反 **98.638 倍**，仍标为不同数值近似产生的不相容摘要，不误称 evidence 归一化错误。

## 输出与范围

Top-K 示例不重新归一化为 1。展示 mass、未展示 modeled mass、未探索状态及未知模型外材料世界分开：声明支持内的 prior-only 事件已全部积分；未指定的更大物理模型没有可量化 posterior mass，因此使用 `null/unquantified`。连续代表材料没有被赋予单点图像概率。

本脚本仅导入 shared material algebra 和 likelihood；没有 physical forward 或 ledger 写入接口，没有读取 locked truth，没有 NN 训练，没有 SHA-256、旧全套 campaign 或额外 smoke 检查。后续必须以实际 full-Mie reference、概率数值误差、mode discovery、fresh parents 与 compensation evidence 独立判断 Joint P2 的科学结论。
