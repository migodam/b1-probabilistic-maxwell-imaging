# 附加诊断：不依赖 MCMC 的材料后验误差下界

本页是在模型、两批评价和原科学账本封存后，对**已经计算的固定样本检验**作出的 post-hoc 数学解释。没有新增 Maxwell、训练、GT 数组或模型选择。使用的是标准数据处理、KL 链式分解和二项区间，不称新 Maxwell 定理。

## 为什么一维残差事件可以否定完整概率质量？

真实联合 law 记为 P(dy,dc)=pY(dy)πG(dc|y)，模型为 Q(dy,dc)=pY(dy)QG(dc|y)。c 是复材料场，允许 plateau 原子和退化支撑，无须假定像素 Lebesgue 密度。

定义事件 E={GammaCDF256(Phi(c;y))>0.975}。在原离散 P1、正确独立复 Gaussian 噪声和精确统计目标下，P(E)=a=0.025。这个性质是 prior-predictive 平均意义，不能拿来要求每一份观测的 posterior 残差都服从同一 Gamma。

记 q=Q(E)。因为两者的 pY 相同，

`E_Y TV(πG,QG) = TV(P,Q) >= |q-a|`。

对事件及其补集作 KL 链式分解，

`E_Y KL(πG||QG) = KL(P||Q) >= a log(a/q)+(1-a)log((1-a)/(1-q))`。

其余两项是事件内/外条件 KL，均非负；若绝对连续性失败，KL 可以无穷大，下界仍成立。全变差采用 sup_A|P(A)−Q(A)|、最大值1的约定。

## 有限 parent 的不确定性

每份新观测预定 index0，一份抽样跨独立 parent 估计 q；不能选择残差最小的生成候选。使用单侧 Clopper–Pearson 下端 L=Beta^-1(0.05/52;k,n−k+1)。52 来自全部预冻结26模型×2批对象，而非仅当前8个实际 Full-PPC 行，保守覆盖从原固定 roster 选择报告子集的影响。

当 L>a 时，平均条件材料 TV≥L−a；平均条件材料 forward-KL≥BernKL(a||L)。在原联合 law、固定模型与**理想独立随机流**前提下，52项同时置信度至少95%。不是条件于已冻结PRNG比特的严格频率学证明，不是浮点定向舍入/连续PDE认证。

|固定 seed1 生成 law|批次|尾事件 k/n|平均条件材料 TV 下界|forward-KL 下界 nat|
|---|---|---:|---:|---:|
|raw Spline|512|510/512|0.95318|3.61276|
|raw Spline|256|255/256|0.93934|3.13445|
|几何 CFM，实际 NFE32 law|512|510/512|0.95318|3.61276|
|几何 CFM，实际 NFE32 law|256|256/256|0.94823|3.41367|

完整8行在 `evidence_summary/CALIBRATION_ERROR_LOWER_BOUNDS.json`。这并不声称每份观测都具有同样误差，而是对原 prior-predictive Y 的平均条件误差。它也不要求 CFM 拥有已合格的显式数值 log-density；界针对实际固定 NFE32 生成 law，不针对尚未数值求准的连续 ODE law。

## 归因与版本边界

初稿只对8行作校正；独立审阅确认代数和数值，在无法仅从紧凑计数证明“子面板在新观测前冻结”的地方提出质疑。主线程采用更保守的52项原固定完整 roster 校正，不改物理阈值、模型、原评分或旧 gate。主 checkpoint 的 pre-observation freeze、两批 online-before-GT 和固定样本身份另由原 chronology 审计支持，不能只靠本页计数代替。

这是在查看结果后推导的解释，不冒称预注册主要检验。原始8项初稿和审阅留在附加审计内。该界与三seed平均 proper-NLL 改善5.568/6.171nat并不矛盾：相对 KL 变好，不表示未知的绝对 law 已接近；这里的下界又只对应指定 seed1。少数边际 coverage 看似接近标称，也不保证 Maxwell response coupling 正确。

若原 forward/noise 定义存在模型失配、parent 不独立或样本选择使用残差，该结论需要重新核验。有限精度、离散求解器验证和真实连续物理仍是不同层级。输出材料风险动作/均值不是 posterior draw，不能被这个 Gamma 事件机械判罚。

这项补充只加强“完整 posterior 未合格”的解释，不取消本轮相对概率进步或独立材料风险收益，也不把标准统计工具包装成 IEEE TAP 级原创定理。

## 有限 K 能否独自解释生成均值不如 CNN？

对原 P1 decoder 的每个 ROI 单元，实增量位于[0,1.8]，虚增量位于[0,.18]；其余320单元固定。分别应用有界实变量的标准方差界，

`trace_M Cov_Q(chi|y) <= (256/576)·(1.8²+.18²)/4 = .3636`。

任意模型 Q 都适用，不要求 Q=posterior，也允许 plateau 原子。K 份条件独立材料样本的直接均值，在总体平方风险中的额外 Monte Carlo 项因而至多 `.3636/K`。本轮 K128 上界为 **0.002840625**。

原恒等式为 `R(sample_mean_K)=R(exact_model_mean)+E trace_M Cov_Q/K`，要求采样辅助随机数在给定Y后独立于GT。它不把Q修成真实posterior。两批rawSpline直接均值MSE约.02278/.02655，相对好的风险动作差约.016/.019，明显大于这个最保守MC项尺度。因此不能把差距主要归因于K不够；模型law/均值本身的偏差仍须解决。

有限测试给出的是风险估计而非精确总体风险；此尺度比较不是每个parent的确定界，也不把事后减去上界的数字当新预注册风险成绩。该结论仅用于直接K128模型样本均值，不应用于K32特征经过非线性CNN的融合动作、短链posterior均值或选择过的样本。
