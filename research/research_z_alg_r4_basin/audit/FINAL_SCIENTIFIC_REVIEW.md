# R4 独立科学审查

最新状态：`POST_V3_RESULTS_REVIEWED`。GT 恢复改善为 `NO`，六方法联合门均 0/3；有限数据目标可达性有受限正证据。首版记录原样保留，最新判决见末节。

审查状态：`HOLD — 首版静态审查，真实新实验未评阅`。日期：2026-10-30。文件名沿用交付合同，不表示已完成论文终审。主代理保留最终科学裁决。

## 当前裁决

理论核心恒等式和两个主要解析反例未发现致命数学错误；它们支持有限离散问题中“局部可辨识/校正”与“有限预算可达”必须分开。它们不推出真实 576 材料几何的失败原因，也不推出候选算法有效。历史诊断对局部非线性与方向旋转提供描述，未提供完整秩、全局 basin 或新鲜泛化证书。

新颖性与性能均保持 `HOLD`。有理扫描是标准 resolvent 工具；本轮最有价值且可证伪的算法假设是多个检出低谷分别校正后选择。代码目前默认校正器与文档 full-J 描述不一致，校正前近似去重也削弱了理想机制。

## 数学抽查

| 项目 | 审查结果 | 条件与具体理由 |
|---|---|---|
| `M_h Δj=D_hQ` | `PASS_CONDITIONAL` | 当前和端点状态矩阵可逆；材料与几何的有限离散范围明确 |
| `E=SM_h^{-1}D_hGM^{-1}D_hQ` | `PASS_CONDITIONAL` | 由 resolvent 差直接得到，非交换顺序正确 |
| 右侧 `M_h=M(I-M^{-1}D_hG)` 与左侧 `M_h=(I-D_hGM^{-1})M` | `PASS` | 两种因式分解在文中分别使用，没有把其范数和顺序混为同一对象 |
| 椭球证书 `sum h_i² ||(GM^{-1})_{i,:}||²` | `PASS_SUFFICIENT_ONLY` | Frobenius 精确量上界 operator norm；严格小于 1 给可逆充分条件，不是必要条件，准备的 N RHS 已明确 |
| 精确 ray `ΔF=t S(I-tK)^{-1}V` | `PASS_CONDITIONAL` | `K=M^{-1}D_pG,V=M^{-1}D_pQ`；只要求合法端点可逆，非 Neumann 小步假设 |
| 固定背景幅度式与导数 | `PASS_CONDITIONAL` | `T_a=D_aGM_b^{-1}` 次序正确；`Y'=L R² B`、`Y''=2LT R³B` 用到 T 与其 resolvent 对易，未假设 D 与 G 对易 |
| 余项/严格下降界 | `PASS_SUFFICIENT_ONLY` | 范数小于 1 及严格上界是条件；正则项实际变化单独计入，不能用经验估计替代证书 |
| 块二阶有理模型 | `PASS_CONDITIONAL` | 基包含 V、KV 时前两系数匹配；残差式含 `t²(I-Π)KU(I-tH)^{-1}C` 且 leading term 消失，故 O(t³) |
| 标量共同 Padé 限制 | `PASS` | 匹配二阶需 `b=βa`；输出二阶共线不能推出内部不变空间或全阶精确 |
| 小活跃集 Woodbury | `PASS_CONDITIONAL` | selector 与自由材料块范围清楚；全部材料 ray 没有被偷换成低秩 |

这些是解析核对，没有加公式精度实验。短基模型尚非当前 dense Schur 实现，不得把其 3s RHS 成本转贴给已运行的全阶 ray setup。

## 两类严格反例的有效范围

`MATHEMATICAL_THEORY.md` 的单材料 `G=1+0.1i,a=2` 例子：原损失导数符号确为 `(x-2)(-1+1.02x)`；下界 0.5 是严格单侧正残差局部极小。Möbius 响应全局单射且实 Jacobian 列非零。它反驳“局部满秩/全局唯一 ⇒ 连续下降可达”，但极小点在边界，不能被改述成高维内点失败定理。

`RESOLVENT_GEOMETRY.md` §7.1 的另一个单材料例子同样成立：`g=1/2+i/8`，导数因子 `(m-1)(1/2-15m/64)`，极大点 `32/15`，上界 20 为严格约束局部极小。两例参数不同且极小边界不同，但没有矛盾。均未证明系数来自本次 disk-pulse Hankel 自作用；“复系数 VIE 代数相容”是恰当范围。被动性的具体物理解释还依赖时间/Green 函数约定，不能从一个正虚部标量升级成 Maxwell 可实现性。

§7.2 tied-material 例子在 epsilon=0 的 `t0=2/3` 确有二阶导 486，非退化驻点允许小阻尼隐函数延续。其真值 t=0、共享材料约束与本次 [.5,20] 全部独立实材料域不同；原文已承认，不能搬成当前全空间内点反例。

`BRANCH_NAVIGATION.md` §5 光滑最小二乘反例也成立：初始 `(-1,0)` 的 Hessian `[[11,-2],[-2,8]]` 正定；区间线上唯一全局最小在 x=-1，x=1 是较高严格线低谷；在 `(1,0)` 全变量 GN 一步到 `(1,1)` 的根。它严格区分 multi-valley-correct 与 raw-argmin-correct，且只证明操作差别。它不是 VIE/Hankel 物理构造。

## 横向收缩的边界

`T_y=m-J†(f-y)` 在精确根、常秩光滑邻域中有 `DT=I-J†J`。满列秩时零导数可支持适当邻域内局部快速校正；欠秩时核方向为单位导数，只有在正则根流形、一致右逆/Lipschitz 等条件下才可讨论到局部根集合的法向距离收缩。根处投影谱本身不等于已证明一个邻域的统一收缩。

阻尼纯数据映射在零残差根的正奇异方向系数为 `mu/(sigma_i²+mu)`，核方向为 1。非零残差驻点必须保留 `sum r_k Hess r_k`；非驻点还有逆矩阵导数项。冻结弱坐标的单位模态来自算法约束；box 变 active set 时映射一般不光滑。原目标加 smooth/TV 后，GT 未必是此迭代固定点。

因此不得据此主张从失败点进好 basin、材料实际沿纤维走动、后验采样或 GT 收敛。

## 历史诊断的证据限度

`BASIN_GEOMETRY.md` 自述为已暴露真值的 offline diagnostic。36 个节点的 full-J SVD 可报告实际计算的完整奇异值数组；目前旋转只摘要 leading 32 右奇异子空间。97.4%–99.5% 能量集中和接近数值精度的尾值，既不证明精确核，也不证明 576 维 GT Jacobian 满列秩，更不证明局部 well-conditioned。**“GT 局部满秩且良态”不得出现在结论中。**

失败终点上的小负梯度探针仍可降低物理损失，因而这些终点未被证明是严格局部极小或 KKT。诊断支持局部预测可一致、下降可能偏离 GT 方向；不能据此确定原正则求解器停机原因。六点 GT 直线是 oracle 描述路径，不是优化器轨迹，也不给实际可用的在线逃逸方向。

Gram 条件、leading-subspace 旋转、resolvent 范数估计和有限步余项均不单独分离成功与失败；四个例子不足以确立通用阈值或因果机制。

## 代码与实际实验待核项目

1. `bridge` 的校正前 RMS=0.25 近似去重可删掉不同低谷；它需要修正或准确记录与收窄声明。
2. ray 校正默认 L-BFGS-B/adjoint；必须按实际方法解释，不能贴 full-J corrector 名称。
3. ray `raw_score` 与最终排序均含共同正则项；有限 scan 的 best-only 是 detected-top1，而非认证 global argmin。
4. `rational_rays.py` 仅取网格邻域局部最小并 bounded refine，无根隔离、窄谷覆盖或全 basin 保证。候选最终 `clip` 如超出仅浮点级修正，会破坏 ray-response 同一性；真实 forward 才是该材料的最终物理资格。
5. 校正返回有限时间内 best visited 原目标点，`corrected=True` 只意味着得到物理候选，不证明收敛或 KKT；候选材料去重也不是 basin 计数。
6. 每个 ray 的 N+s retained-LU RHS、全阶 Schur 与每节点三角 RHS 已在代码中登记；真实端到端收益和预算 overshoot 需看结果，不能从“不新 LU”推出便宜。
7. 分槽规则给后续候选分配剩余时间，前一候选提早结束会改变后续实际时长。这是确定的共享剩余预算策略，不能称严格相同实际校正时间。

## 当前证据状态

| 裁决 | 状态 |
|---|---|
| 有限离散核心代数 | `PASS_CONDITIONAL` |
| 可辨识不等于下降可达的逻辑分离 | `PASS_IN_STATED_ALGEBRAIC_SCOPE` |
| leading-32 诊断推断全576良态满秩 | `NOT_SUPPORTED` |
| 真实 Maxwell 多低谷校正收益 | `HOLD / RESULTS_NOT_REVIEWED` |
| 同预算优于普通多起点 | `HOLD / RESULTS_NOT_REVIEWED` |
| 全局 basin、全局收敛、连续 Maxwell 推广 | `NOT_ESTABLISHED` |
| 文献最近邻与全球原创性 | `HOLD / COVERAGE_INCOMPLETE` |

后续只评主代理实际交付的冻结结果与必要追踪，不扩大预算来挽救主张，不把历史开发成功改名为新鲜 holdout 成功。

## V3 新鲜结果独立追加复审

本次定向追加读取最终 `SCORES.json`、fresh contract、预测封存记录和最终算法说明；没有新运行、额外搜索、hash 或重复一致性校验。封存记录声明 18 个预测全部关闭后才读新鲜真值；评分 JSON 同样明确 truth-first-read-after-predictions-closed。这支持按所记录协议作为冻结后 fresh evaluation 解释，但本审查没有审计生成器隔离进程，因此不把独立文本评阅称为完整泄漏审计。样本数只有三个，同网格合成模型，且最终设置已吸收 V1/V2 开发负证据。

### 冻结结果

| 新样本 digit | detected-best-only 残差 | uniform_bridge 残差 | ordinary uniform_multistart 残差 | bridge top1 RMSE | bridge oracle best-K RMSE |
|---|---:|---:|---:|---:|---:|
| 7 | .762235 | .015425 | .054634 | 10.8908 | 5.0503 |
| 8 | .780809 | .111464 | .137025 | 10.5933 | 6.2721 |
| 2 | .911384 | .045760 | .041534 | 11.3274 | 6.3376 |

全部六方法的 joint gate 均 0/3；bridge 的残差仍全部超过预声明 .001 数据门，RMSE 全部超过 1.5，IoU 为 .1553/.1804/.2405，全部低于 .7。数据目标显著降低是连续量改善，不能被写成达到数据恢复门或成像成功。

bridge top1 RMSE 比保留的共同 anchor 更差；各 archive 的 oracle best-K 只是保留 anchor，新增候选没有材料 RMSE 改善。digit2 的 bridge best-K 约 6.337567 还略差于 data_lbfgs 的 6.337132；不能把这类数值近似并列写成材料收益。oracle best-K 是评分后诊断，绝不是在线可用策略。

### 机制判断与实际边界

uniform_bridge 相比同扫描 detected-top1 消融，在三例中均产生更低数据损失的实际材料候选。因此首版 `MECHANISM_BENEFIT_HOLD` 可以被细化为“对有限数据目标的可达性，存在限定正证据”。这不证明真实 GT basin 可达，也没有证明严格局部极小之间发生跳跃，因为 anchor 停机不等于 KKT/严格极小证书。

ordinary uniform_multistart 在三例也大幅降低数据残差，并在 digit2 优于 bridge。bridge 在 digit7/8 残差更低是该冻结策略的有限结果，不能升级为总体优于多起点；单个候选的校正预算分流、数据相关 seed 选择与实际停机共同影响结果。没有统计泛化证据，也没有成本优势定理。

实际 bridge 在 digit7 用 22.11s、2139 LU、51876 full complex RHS；digit8 用 30.69s、2478 LU、60000 RHS；digit2 用 28.15s、2479 LU、60000 RHS。后两项明确为 `BUDGET_STOP`，只校正了三个候选，虽计划最多四个。ordinary multistart 三例均达到 2500 LU，亦为 `BUDGET_STOP`。0 LU/RHS overshoot 的记录不能被解释成所有预定 seed 均完成；未完成候选和停止原因必须保留。Schur 与三角传递工作仍在账本中，不能只凭 LU 次数给出完整计算效率排序。

V3 文档已修正校正前阈值至 `1e-7`，最多四个幅度分散低谷、wall/LU/RHS 三重 seed 配额，以及已付费 anchor 复用。L-BFGS 使用仿射变量缩放，GT 未参与该变换；这改变优化路径及绝对停止尺度而非物理模型。首版 `.25` pre-dedup 的批评只属于历史快照。仍有效的限制是：网格低谷不完备、截断候选不代表全部 basin、corrected 不代表收敛、正则 ranking 不等同纯数据 ranking、局部定理不等同实际正则校正映射。

### 最终独立裁决

| 主张 | V3 后独立裁决 |
|---|---|
| 相比单个 detected-top1，多 seed 校正可降低原数据目标 | `LIMITED_POSITIVE — 3 frozen synthetic cases` |
| 达到预声明数据恢复门 | `FAIL — 0/3` |
| 成像/GT 恢复改善 | `NO — top1 worsens; no new best-K material gain` |
| 联合成功 | `FAIL — all six methods 0/3` |
| 总体优于普通 multistart | `NOT_ESTABLISHED — digit2 counterexample within this batch` |
| 标准 resolvent 成为新理论 | `NOT_SUPPORTED` |
| 全局 basin、576 良态满秩、连续 Maxwell、概率解集 | `NOT_ESTABLISHED` |
| 全球新颖性 | `HOLD — literature coverage remains incomplete` |

可接受的研究结论应把这次工作定位为：标准 resolvent 结构组织有限预算候选的算法研究与失败成像证据；它揭示降低当前数据目标与恢复真实材料之间的明显分离。不能为挽救 GT 主张扩预算、改 gate 或把低残差错误材料解释成已验证等价物理解。
