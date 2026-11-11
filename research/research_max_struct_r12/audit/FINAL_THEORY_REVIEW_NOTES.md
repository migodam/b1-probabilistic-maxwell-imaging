# R12 数学与实现范围：供root最终审查

本记录阅读了更新后的 [THEORY_DERIVATIONS.md](../THEORY_DERIVATIONS.md)、实际 [field_current.py](../src/field_current.py) 和已保存NN/high3D结果。不修改core、不跑新资格实验。以下是最终稿必须核对或明确的点；不是worker代替root的最终理论/novelty裁决。

## 1. conditional shadow-gap bound：代数核对通过，条件不能删

在当前符号下，\(K_m=I-D_mG\)、\(e=Q_c-E^{inc}-GJ_c\)。由 \(J_c=D_mQ_c\) 有

\[
K_m(J_c-J_{phys})=D_me.
\]

因此稿中的 \(J_c-J_{phys}=K_m^{-1}D_me\) 符号和次序正确。由 \(Q_a=E^{inc}+GJ_a\)，\(e=(Q_c-Q_a)-G(J_c-J_a)\)，block norm给出 \(\|e\|\le\sqrt{1+s^2\|G\|^2}\,g\)。再用 \(\|J_c-J_a\|\le sg\) 得到稿中的port bound；没有发现该条件推导的代数冲突。

root须最终保留：

- \(g\) 是真实 \(a-c\) 的scaled product Frobenius norm，不是 \(z-c\) 或未缩放current差。
- G、port及resolvent的范数使用与Frobenius/vectorization相容的induced 2-norm；源特定ports对应block-diagonal port norm。不同源共享同一G/材料响应是当前物理前提。
- \(K_m\) 可逆是条件，不是对整个legal box已证明的性质。H正定不推出材料Maxwell矩阵良态。
- 未独立控制 \(\|K_m^{-1}D_m\|\) 时，该界不是uniform或可直接计算的验收证书；不应为了包装结论新增重复norm/机器精度实验。
- \(g=0\) 时的“uniqueness”是给定m的有限物理场唯一性，不是材料唯一性，更不是连续实验有效性。没有GT距离界。

## 2. 必须注明iteration0 audit是baseline占位

`run_projection` 的 `record(0,...)` 将初始物理state同时放入A/C字段，`gap=0`，但没有调用A投影。若初始数据不在球内，那个A字段不属于A。

这不影响随后AP/DR迭代，却影响audit解释：iteration0的0 gap不能用于conditional bound、feasibility或convergence判断。方法文稿已将其标为baseline/占位；最终图表或下游分析同样须排除该占位gap。若root以后修正metadata，应保留原记录/版本，不把它当科学轨迹变化。

## 3. 第一更新恒等式与暗cell材料标签

完整物理场初态在C中，故最近场对零距离且 \(P_Cz_0=z_0\)；omega1时AP/DR首次更新均为PAz0。这个场状态恒等式成立，暗cell材料标签不唯一不会改变零距离场对。

但“不变”不可扩写成材料标签和prior必然不变：若一个cell全部源/分量的x=q=0，确定性box tie-break可返回另一材料标签。理论中错误球内物理解固定，应明确是场状态；相应材料读出是否同一m、该m的Maxwell可逆性和prior仍由物理审计处理。omega≠1时第一次DR是relaxed PA，不能继续宣称与AP完全相同。

## 4. C投影精确性的适用范围

有限实box、共享实标量、Euclidean scaled product metric下，消去q后的有理二次式和合法根/端点枚举足以给全局cell nearest pair。共享范围必须包含全部source和三个component，cell-major布局不可省略。

该“exact”不涵盖：复材料、各向异性tensor、额外非局部prior、不可微branch切换。对2×2 Gram退化、端点切换或非零DR法向附近，不能不加条件地套用光滑tangent-space局部速率。固定点shadow条件与全局迭代收敛是不同命题。

## 5. A投影：理论精确与浮点rank判定分开

为消除 `T=WsS` 的阅读歧义，建议最终稿统一写 \(T_j=sW_jS_j\)、\(y_j=W_jY_j\)，每源可有mask/white block但只共用一个lambda。跨源相关white未实现，不由当前代码暗示其支持。

数学上H=I+s²G*G正定；port Gram可有零特征值。球低于null残差floor不可达，等于floor通常需要lambda无穷的伪逆极限。这些条件已在代码处理，最终稿不宜写“总有唯一有限正根”。

浮点实现将小于 `32*eps*max(d,p)*maxeig` 的正特征模式视作数值零，并用 `128*eps*energy_scale` 判边界。它是带明确rank/容差的dense数值实现，不是无限精度的集合投影。理论精确式与实际数值rank近似必须并列；近floor结论尤其不能超出该口径。

grouped v2只对完全相同T分组、多RHS应用并聚合power；它没有改变数学投影。但求和顺序会变，局部output差约1e-15的一个benchmark不保证3000步非凸轨迹bitwise一致。已启动v1结果保留、新cohort预先固定v2的版本边界必须保留。

## 6. 噪声合同与两种原物理目标

raw iid complex的99%全球半径 \(\sigma\sqrt{\chi^2_{2n}(.99)/2}\) 与本轮合同一致。API任意W/mask不能自动沿用同一个rawδ；需要相应covariance/rank合同。真实mask/Fresnel也不能套用synthetic共享S。

scalar-TM使用 .5全球norm平方归一的数据项加原smooth/TV prior；R8 vector使用逐源norm平方均值（无.5）加quadratic prior。两者都叫原Phi，但数值不跨模型比较。AP/DR raw球与vector candidate source-normalized objective不同，这是明示设计，不应隐藏成同一metric。

## 7. 已完成high3D提供的是selection边界，不是梯度错误证明

matched DR raw物理CN99为true，随后同Phi corrector使D升、R降、Phi降、材料误差升且出球。它符合soft D+R优化，不与A投影约束矛盾；corrector不在A∩C上求解。

mixed AP还出现D/R同时降、材料误差升，说明“把数据改善当GT证书”同样不成立。cross-grid clean差异26.1567%、coarse cell phase大，必须保留finite-DDA/model mismatch限定；真coarse标签出球不证明不存在任何其他合法球内coarse解。

baseline raw exact physical current未存，但其已知source-normalized D给出

\[
P D\min_p\|Y_p\|^2\le\|r\|^2\le P D\max_p\|Y_p\|^2.
\]

六个baseline raw下界均大于delta，故球外结论可以由界认证；缺失的精确global residual仍不能补造。RAW_CORRECTOR_DECOMPOSITION保留此区分。所选corrector evaluation序号准确，accepted iteration未保存；不从总nit猜测。

## 8. secondary规则需root完成独立确认审查

规则是已付费审计池上的 `(R,Phi)` 词典序minimum，受legal/完整物理/CN99/cap资格限制。它不是对连续可行集求全局minR，也不凭空增加候选或可辨识信息。无合格点应unavailable，不能借primary兜底后混算secondary成功率。

注册时序应按真实文件：24观察生成后、确认inverse开始0及GTscore0时注册；primary/core更早在数据生成前冻结。两者不能合并称“全部生成前冻结”。动机来自高3D已暴露开发，当前不能叫secondary验证成功。

root最终还须核查：八方法hook一致、完整物理pool而非nonphysical shadow、phase/atomic cap资格、快照/检查费用计入、无候选率、失败与缺失阶段池覆盖。候选覆盖不同，潜在selector收益不能等同纯DR机制收益，也不能用后验GT最佳pool点当预测。

## 9. 最终novelty和机制归因边界

DR/TLS/固定Schur消元/resolvent传播都是已有工具，不把此处代数或单个工程速度称新定理。AP/DR同init同projectors才能隔离反射；NN组没有NN→AP，DR/CIE还改变current表示/辅助目标/费用。high3D共同corrector后相近端点和明显较高DR费用须保留。

NN四例为暴露开发；高3D是两个独立几何加一个相关配对，且coarse highindex未给continuum精度资格。root应结合独立24整组结果、closest work及访问限制决定贡献范围，不从这些完成开发结果直接推出普适盆地、普适NN prior鲁棒性或全局恢复。
