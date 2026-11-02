# R2 最终独立理论与主动照明表述审查

2026-11-02。只读审核主线程 `THEORY_DERIVATIONS.md`、`experiments/active_family.py`，交叉核对 `src/profile_current.py`、`src/cie_global.py`、R1已资格DDA及既有文献证据。按主线程授权只修正显式归一化/符号、主来源链接和非数值结果措辞。**没有运行新实验、导数测试、反演、评分或隐藏GT读取。** 二维预测已封存，本审查不假设算法优势。

## 审查结论

**核心有限矩阵论证在所列假设下未发现数学错误。** 该结论是推导/实现接线审查，不是成像收益、统计显著性、连续PDE或全局恢复PASS。下列前提与限制必须保留。

### 1. 消元后可变度量与源归一化

`H=P(m)L` 可逆需要L可逆及P各对角非零/有限，故β与 `1+βm_i` 的合法无极点条件不能省略。在本轮total-energy目标下，明确

`Sbar=S/||Y||F, ybar=y/||Y||F, Hbar=P L/||E||F`，

且两分母非零。以这些归一化量写 `K=Sbar Hbar⁻¹` 才得到与 `src/profile_current.py` 一致的profile identity。原文“吸入固定权重”可理解但不够显式；已按授权补充该定义，随后S/y/H均约定为本节简写。

`z*=−(I+K*K)⁻¹K*r`、Woodbury代回及谱界正确。多源共同标量归一化时逐列相加/同一空间factorization成立；不同Ωd/Ωs时需generalized Sylvester或分源加权系统，理论§2已经正确警告。不可将此profile与原生source-normalized CIE当作同一个辅助目标。

`profile`是最小化潜状态后的值，非积分边缘似然；归一化/行列式缺口已正确说明。κ=||K||无统一上界时profile小不统一控制原物理residual；标量反例只是矩阵反例、非已构造Maxwell材料，现文明确保留这个边界。小gradient也不补足缺失的resolvent条件。

### 2. 仿射current容量差

固定m、固定共同目标，若 `A=S*S+H*H>0`，完整current最小点的一阶项为零，故加权平方分解成立。其对仿射集合取下确界即给profile差。current系数冗余不妨碍physical-current投影结论，但若进一步对系数伪逆求导，需要常秩条件；本节没有错误地假设native Fourier系数全独立。

嵌套current空间使profile不增；不能推出材料误差不增。不同J_D偏移、不同β、不同正则器或不同源权重不属于这个定理的同目标嵌套比较。文中限定正确。

### 3. 包络梯度

对实材料，`w=β/(1+βm)` 的导数 `w'=−w²`。固定J时 `Q=E+GJ`，`de_i/dm_i=−Q_i`，故

`d(we)_i/dm_i=−w_i² e_i−w_i Q_i`。

因此理论§3与 `src/profile_current.py` 的实部共轭求和、`||E||F²` 分母、共同材料regularizer梯度一致。该包络梯度以current确切极小化为前提；浮点Cholesky仍是数值实现，不把常规舍入提升为严格机器证书。理论原来未显式定义Q；已按授权补充。未执行新的方向差分。

### 4. 完整源酉等变

`Ω'=U*ΩU`、J'=JU使固定目标不变。单位源旋转保持Frobenius梯度/内积、global PR scalar、全局exact line、共享材料readout、初始化/zero-padding和全局停止比。`src/cie_global.py` 使用full right metrics、Cholesky右白化、global变化比；与声明一致。

只运输objective不能使逐列PR或逐列相对stop等变；R2将这些区别明确保留。非酉变换还需变量内积/梯度运输；当前理论仅宣称酉变换。等变不保证更好basin、耗时或GT恢复。

### 5. 有限库与预测覆盖κ

有限库内任意两个预测相距>2ε，直接由三角不等式得到同一观测的ε球内至多一候选。如果真响应只在某候选的**当前采集预测空间**κ邻域，则半径扩至ε+κ确保至少一个覆盖候选合格；所有pair>2(ε+κ)给至多一个。文中“预测κ邻域”明确，因此命题正确。

κ必须对应实际所选采集、相同白化与范数；它不是未经稳定性证明就可由材料距离替换的数。没有覆盖界，候选内分离不保证覆盖真值。连续t族的端点最优source，也不保证全区间minimax或单调性/注入性；主文已明确。

### 6. 已知轴块化与主动实验

`ε−I=t b zzᵀ`令j_x=j_y=0，完整场仍为Einc+G_{:,z}j_z；z分量给

`j_z=(I−D Gzz)⁻¹ D Einc_z`。

在所需可逆性下这是精确块消元，不要求Born。沿z传播的两个横向旧源Einc_z=0，故整族旧散射相同（零）。新源的四列线性组合由已付完整forward得到；生成代码确实另外调用该公开DDA生成观测并加确定性有界噪声。它是独立synthetic realization，非独立solver或实验室采集；所有六个计算RHS已在各阶段计费。

fit只读取 `data_public/active_observation.npz`，score是另阶段读取封存参数；源码未将score误差用于选择拟合candidate。声明的真实值采样落在公开[0,3]族内。coverage control增加合法xx分量，所有受限候选旧源预测均零，所以若observed norm超过旧噪声半径，整个受限族都被拒绝；这是覆盖失配而非优化失败。

一份bounded scalar fit及有限参数误差**不能**认证连续族全局唯一性或所有t都可恢复。原 `SCORE.interpretation='added illumination resolves restricted family...'`强于现有证据，尽管FIT已有非全局认证限制。已按授权将源码及既有SCORE.json措辞改为：

`one realization compatible scalar estimate; no uniqueness certificate; finite DDA only, not CIE advantage`。

THEORY§5也改为“一个公开族内实例新增观测后的相容标量估计”，明确不证明全局可辨识、拟合唯一性或噪声恢复保证。

## 已执行的非数值修改

1. `THEORY_DERIVATIONS.md`：补Sbar/ybar/Hbar与非零归一化分母；补Q定义；限制主动实验结论；添加原始DOI/arXiv链接。
2. `experiments/active_family.py`：仅改SCORE解释字符串，没有改变模型、随机种子、源选择、误差、停止、预算或运算。
3. `results/active/SCORE.json`：仅替换解释字段，新增 `editorial_correction`保存旧解释及理由；所有数值字段原样保留。没有调用score或任何solver重跑。

## 原文归属

profile可变度量：[van Leeuwen–Yang arXiv2403.15292v1，Cor1 Eq10](https://arxiv.org/html/2403.15292v1)；exact state profile/envelope/Schur：[van Leeuwen–Herrmann DOI10.1088/0266-5611/32/1/015007，§4](https://doi.org/10.1088/0266-5611/32/1/015007)；二次投影：[Osborne ETNA28(2007)，pp2–3](https://etna.ricam.oeaw.ac.at/vol.28.2007-2008/pp1-15.dir/pp1-15.pdf)；inexact inner：[Aravkin等arXiv1606.02395](https://arxiv.org/abs/1606.02395)；mature CIE：[Zhong–Xu DOI10.3390/jimaging5020027](https://doi.org/10.3390/jimaging5020027)；相邻测量域state/data失配：[Sun等DOI10.1109/TAP.2017.2673758](https://doi.org/10.1109/TAP.2017.2673758)。这些链接已放入theory正文，证据等级/身份/版本与限流见 `literature/PRIOR_ART.md`。

最终算法收益、有限对象覆盖、实验独立性强度、论文贡献及发布结论仍由主线程结合封存后结果判定。
