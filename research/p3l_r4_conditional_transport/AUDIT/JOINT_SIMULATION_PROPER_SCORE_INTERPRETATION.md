# 不需要MCMC真值，也能判断哪些概率或风险误差确实减少

主线程数学复核。以下是标准条件期望、KL及proper score性质，不称为本轮原创定理。本轮的新证据必须来自冻结算法之后生成的独立P1对象；历史/开发分数不能替代该检验。

## 1. 材料MSE比较的含义

以复材料的实Hilbert损失 `||chi||_M²=sum|chi_i|²/576` 定义真正条件均值 mu(Y)。若动作T只依赖Y、已冻结模型及与真值条件独立的计算随机数，则

`E||T-Chi||_M² = E||T-mu(Y)||_M² + E||Chi-mu(Y)||_M²`。

因此两个动作的总体风险差等于它们到真正Bayes条件均值的L2误差差；不可约条件方差相消。合成联合模型中，每个新parent的一份GT足以无偏估计这种总体风险，不需要先有该parent的MCMC teacher。

这不证明某个动作已经等于mu，也不保证它对每一份GT都更好。有限测试的不确定性、训练seed与统计单位必须保留。

## 2. Full66 proper NLL比较的含义

原P1在Gaussianized x坐标有普通dx66密度。对两个完整正密度q1、q0，在积分有限时：

`E_joint[-log q1(X|Y)+log q0(X|Y)] = E_Y[KL(pi(.|Y)||q1(.|Y))-KL(pi(.|Y)||q0(.|Y))]`。

所以冻结模型在独立joint samples上的proper NLL改进，的确估计了**相对的平均forward-KL改进**。没有完整posterior reference并不让这个比较失去概率意义。未知的真正条件熵在两模型之间相消。

绝对KL仍然未知。NLL较好的模型仍可能严重不校准、Maxwell预测不相容。因此不能把相对改进升级为FULL posterior合格。CFM的数值log-density若没有合格误差检查，不得套用此结论；本轮冻结NFE32的CFM主要使用样本proper scores。

## 3. 材料energy score比较的含义

令P_Y是真正材料posterior pushforward，Q_Y是某个完整生成器的材料pushforward。用576复材料实化并除24后的Euclidean norm，

`ES(Q,c)=E_Q||X-c||-.5 E_QQ||X-X'||`。

条件于Y，真实P_Y下的期望分数满足

`E_P ES(Q,C)-E_P ES(P,C)=.5 D_energy(P,Q)`，

其中 `D_energy(P,Q)=2E||P-Q||-E||P-P'||-E||Q-Q'||`。Euclidean energy distance对有限一阶矩的分布严格proper，包括材料plateau原子/退化支撑；不要求像素Lebesgue密度。

独立joint测试的分数差因此估计平均材料conditional energy distance之差。128样本的不同索引pair估计用于避免样本自对项偏差。只比较不同seed平均损失时，parent仍是统计单位。

分数降低16%不等于未知的posterior energy distance降低16%，因为真实posterior的基准分数未知。它也不判断latent fiber内的概率是否正确：两个不同latent law可能经decoder得到同一材料law。Full66 NLL、材料proper score与数据依赖的Maxwell校准必须分开报告。

## 4. 有限K材料均值的附加误差

若K份生成样本条件独立，样本均值的材料风险等于真正模型均值风险加 `E trace_M Cov_Q/K`。本轮主分数保留原来冻结的K=128计算动作，不能在测试后换K或选择样本；可以另外从已存样本估计MC贡献，解释它是否足以解释成像误差。该诊断不修复错误的模型law。

## 5. 本轮应采用的科学措辞

- 风险主线：独立P1测试下的材料风险改进；不是每物体恢复保证。
- 概率主线：完整66维密度的相对proper-NLL/KL改进，或材料生成law的proper-score改进；不是未核验的完整posterior质量。
- 校准主线：fixed-draw Gamma256、数据依赖SBC、相事件与连续材料投影分别核验。任一必要检验显著失败都必须披露，不能被平均NLL掩盖。
- 物理代表/几何archive：只有明确的数据、先验与材料分离资格，没有归一化的posterior mass时权重为null。
