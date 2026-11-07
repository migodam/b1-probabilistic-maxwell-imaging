# 吸收型分母：由passivity得到的有界current变换

状态：根线程DERIVED / PROVED_UNDER_ASSUMPTIONS；候选成像效果尚未运行。该候选是在16确认修复失败、物理端点壁垒诊断之后提出，不改变48或50三维的冻结算法。只做曝光DEV，不将旧对象重新称盲测。

## 公式与证明

约定exp(−iωt)，出射有限Green满足Im_H G≥0；实材料m_i>0，D=D_m，S=D⁻¹−G。这里S是逆contrast减Green，绝不是接收矩阵。于是Im_H S≤0。对任意β>0，令

`A_abs=I−iβD−DK`，`r_abs=β A_abs⁻¹[J−D(Einc+GJ)]`。

材料仍是原实box内m，人工吸收只用于辅助残差分母；原物理F(m)、数据、prior、接收、源、最终评分完全不变。β外面的缩放仍为正实β，不把其平方误写成(−iβ)²。

先取K=G。此时A_abs=D(S−iβI)，current导数

`H_abs=β(S−iβI)⁻¹S`。

对任意v，

`||(S−iβI)v||² = ||Sv||² + β²||v||² − 2β Im(v*Sv) ≥ ||Sv||² + β²||v||²`。

因此S−iβI可逆，`||(S−iβI)⁻¹||≤1/β`，以及`||H_abs||≤β`。不要求m均匀，不要求G normal，也不需要与D同时对角化。这给出的是固定材料时完整current残差导数在物理current Euclidean空间的上界，不是native非正交Fourier系数的Hessian条件数界，更不是反演收敛或材料恢复定理。

等价地，`C=(S+iβI)(S−iβI)⁻¹`是收缩，`H_abs=β(I+C)/2`。这是经典耗散算子Cayley结构的有限矩阵应用，不宣称新的Cayley定理。相比实移位S+βI，吸收型移位利用了已知Im_H S的单向符号，而不依赖Re_H S为正。

若K=PGP且P为Hermitian正交投影，则Im_H K≥0。令S_K=D⁻¹−K，仍有上式resolvent界。实际完整G current导数为

`β(S_K−iβI)⁻¹(S_K+K−G)`，

从而`||H_abs,K||≤β+||G−K||`。这只是保守充分上界，未包含方向信息；不能把K截断后仍叫完整G保证，或宣称它在每个通道更好。对m≥0中零单元，可由消去零行及m→0的连续极限扩展有限模型current界；DEV m≥.02无需此扩展。未将该论证外推到任意非uniform-cell离散、未知μ或一般复材料。

## 正确材料导数与伴随接线

令c=−iβ，w=A_abs⁻¹e，a=A_abs⁻*(wΩ)，e=J−D_mQ。在固定J/Q时

`dw[h]=−A_abs⁻¹D_h[Q+(cI−K)w]`，

`g_i=−ρβ² Re Σ_s conj(a_is)[Q+(cI−K)w]_is + grad R_i`。

R10的Woodbury代码在adjoint末尾除以实a_diag=1+βm；扩展到复c时必须改为除以conj(1+c m)。直接把β换成复数而不修这个伴随是错误实现。current外缩放保持实β，因此原H-adjoint中的该标量不需另取复杂符号。

## 最小决定性DEV

固定shift_phase=−i，无相位/幅度扫描，沿用原β_j/β_m大小。两个新臂：absorptiveGreenBothBlocks（原19维Galerkin K）及absorptiveScalarBothBlocks（K=0）。后一臂用于区分吸收型标量收缩与非局部Green作用，不冒称成熟原论文CIE。两者使用同8步完整材料L-BFGS读出、同real-CIE seed候选（仅按新辅助目标接纳）、同global PR和原FFT/APCS、同prior/box/初值与共同原Maxwell corrector。

四个固定曝光对象：R10环0001、开缝壳0004、桥0022，以及刚揭示BothBlocks破坏的新16环0000（现在明确为曝光开发对象）。共8次，whole-child480 CPU上限，不占远端48确认或三维任务，不以GT选根。先完成新复伴随/真实材料梯度的一次必要接线核对，不重放已冻结恒等式。8个预测全部封存后再离线评分。

剩余问题是：受控的current增益是否能保留环收益并避免其它对象的错误box路径，还是只改善算子界而不改善材料。这是本试验唯一扩展理由。即使成功也只是DEV，不能改写当前确认结果；若失败，不扫描更多复移位来寻找正例。

## 先例边界

复移位作为波动方程预条件已有成熟先例：[Erlangga, Vuik, Oosterlee (2004)](https://doi.org/10.1016/j.apnum.2004.01.009)，其原问题是Helmholtz线性系统的shifted-Laplace预条件。本文候选在材料反演辅助残差中使用耗散Green符号，不因此宣称发明复移位，也不把线性正演求解效率直接当材料恢复优势。当前是定向先例检查，不是完成整个相关领域的新颖性排查。

更直接的限制是：已有[2020 CIE/WT-TSOM论文](https://doi.org/10.3390/electronics9111760)明确允许复β。因此K=0控制只是已知复CIE族中的固定选择及共享实材料适配，不能叫新算法；纯虚移位取其边界情形，资格由上面的当前物理符号证明而非直接照搬论文约定。若所有收益都被K=0解释，必须归给这个已知控制，不能归给非局部Green。

能量不等式与收缩预条件在[电磁FD正演](https://doi.org/10.1093/gji/ggw237)已有成熟研究；[convergent Born series](https://doi.org/10.1016/j.jcp.2016.06.034)也有人工吸收/预条件的先例。[Jakobsen的2026预印本](https://arxiv.org/abs/2604.20282)进一步给出CBS的Cayley/resolvent解释。这里仅核对了后者摘要与元数据，不声称审查全文证明；其正演结果不能直接移植为本轮实材料反演定理。本轮潜在增量只能是这类结构在非局部材料读出中的具体效果、成本与失败边界，而不是“首次发现Cayley/复移位”。
