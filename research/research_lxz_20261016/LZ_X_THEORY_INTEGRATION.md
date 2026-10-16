# L/Z/X统一对象、但不混淆结论

MaterialPrior定义真实参考测度、原decoder与材料约束；MaxwellTarget固定同一forward/noise/nuisance/likelihood并计费；LocalGlobalSampler仅改变合法proposal与接受kernel。P0、条件smoothP0、未来P1、radial proxy及fullMie均独立identity。

L提供prior-informed tangent与有限噪声局部几何，要求chart光滑、材料方向非零、normal rank条件，不能把max归一化冗余当解族。Z提供全局latent proposal；由y确定的physics context没有新增观测信息，只可能改善算法表示。X负责固定材料的合法likelihood及密度修正，corrector不进入posterior kernel。

本轮local使用x=a+wΦ(u)，uniform-x prior推回standardGaussian，冻结R对称且RRᵀ+BBᵀ=I，因此pCN接受只需likelihood比。Global在dx下用q=.9normalizedC1+.1prior，接受min(1,Lnew qold/(Lold qnew))。固定.25混合两合法kernel保持target；状态相关selector需要反向选择概率比。TemperedSMC依次πβ∝p0L^β，带权重、重采样、prior-reversible mutation。

条件smoothP0的normal repair只做candidate discovery；不把多对一修正后的材料频率叫posterior。A17低秩Schur/feedback identities是可选加速关系，本轮没有合格实际增量，不强制拼成GN+Flow+MH。
