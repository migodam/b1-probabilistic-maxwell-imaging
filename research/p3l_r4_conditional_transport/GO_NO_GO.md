# 最终科研裁决

|问题|裁决|实际支持与边界|
|---|---|---|
|GO-A 材料风险成像|**GO**|512/256新对象相对R2改善13.10%/11.36%，相对强CNN改善5.36%/4.90%；普通Gaussian fusion一样有效，限定监督风险动作|
|GO-B 相对完整条件概率|**GO（相对）**|同容量Spline proper NLL改善5.568/6.171nat；静态Maxwell三角优于随机frame；几何CFM样本energy score改善约21%|
|GO-B 已校准FULL posterior|**NOT QUALIFIED**|Gamma残差严重失败、IS中位ESS1、校正短链不合格、CFM数值logq不稳定；不能给完整可信区间或mode mass|
|GO-C Maxwell计算|**GO，限定固定核**|quiet新8对象重放同决定同路径，wall降66.17%，含setup；非连续PDE认证，非混合增量；普通GMRES也有效|
|GO-D 静态几何学习|**CONDITIONAL GO**|NLL、energy score及小幅风险消融有重复改善，但不同任务指标并不全面一致|
|GO-D 动态Lie bracket与当前learned-J|**NO-GO 当前实现/规模升级**|距离保持差、路径昂贵、无材料统计增量；surrogate不过误差门，动作排序弱于便宜普通更新|
|GO-E 强IEEE TAP核心创新|**尚未建立**|有实质算法与独立实验，但通用CNN融合、Flow、LIS、resolvent和retrospective bounds均有成熟先例；没有新的Maxwell专有核心定理与跨尺度有效性|

停止扩大：当前相阈值mass head、active-row exact更新加速、四段路径、当前完整J回归，以及以反向KL强迫response-fit的当前参数。保留完整负结果，不改目标/噪声救gate。

下一步最有辨别力的正面实验：在同一P1/新独立对象上，构造既保留proper-NLL覆盖又把固定样本Phi降到原噪声尺度的可逆full66条件输运，并与同容量普通模型比较；或者在更大二维系统中，以相同统计误差/决定正确性证明包括setup与refresh的端到端收益。它们需要新的冻结合同，不能把本轮未用额度当永久授权。

风险成像已经值得继续；它应定位为物理基础图上的amortized material-risk estimation，而不是已解决的Bayesian posterior。A17在denseTM24没有独立净优势，删除其强制依赖不影响其他路线。
