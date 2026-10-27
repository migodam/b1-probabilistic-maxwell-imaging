# THEORY_AND_COUNTEREXAMPLES §2 的独立条件性审查

只读取付费A已保存的 `COMPENSATION_VECTORS.npz` 与 `COMPENSATION_BASE_CACHE.npz`；无Target/solver import，无物理RHS。原向量获取40 RHS保留为旧费用，不计本审查新增。

## 核心论证

固定其余材料而使ROI全体材料沿同一个复标量v变化，有限VIE系数矩阵对v仿射。只要矩阵可逆、源/接收映射不随conjugate(v)变化，逆矩阵及F(v)在该开域全纯。写a=∂F/∂v，则实部/虚部参数导数分别是a、ia。任何保留完整复通道的复线性selection P都满足P(ia)=iPa；普通实化后的两列范数相等且内积Re〈Pa,iPa〉=0。因此JᵀJ=||Pa||²I₂。在Pa≠0处是局部rank2浸入；Pa=0时rank0。

这一论证支持报告所写的**条件局部结论**，不是全域单射、保距、有限扰动旋转对称或一般非均匀材料场各向同性。a及其范数可随v变化，离开非极点域也不能沿用。

## 付费向量机械复核

- 保存的Jmaterial与BASE_CACHE的原J64、J65/.15逐元素相同；q也相同。
- J第二列与第一列经复乘i的实化旋转相符至float64误差；c与−κ(ReD·Jre+ImD·Jim)复核一致。
- 原1/2/4/8 sources整组复通道、额外偶数receivers、sources{0,2}×receivers{1,7,19}的6种saved-vector子选均满足Gram比例单位阵、两属性补偿余项近零。
- 单实属性补偿的余项比均为 |ImD|/|D|≈0.0809831321；两属性自由时c落在其二维切空间内。该比例对这些完整复通道子集不变，是共形结构的代数结果，不能把6次子选算6份独立物理证据。
- 对已存向量施加相同real/imag通道权重仍共形；只把real分量乘2通常破坏共形。完整数值与Gram见JSON。

## 需要保留或澄清的条件

1. **各向同性是充分条件而非必要条件。** 任意复线性whitening，包括一般proper circular Gaussian的Hermitian协方差白化，仍与乘i交换。报告“改变噪声协方差时不沿用”可保留为谨慎边界，但更精确是“必须重新验证whitening与复结构相容；一般real/improper协方差不保证”。
2. **参数尺度要统一。** Gram式针对v的Re/Im材料增量，不能直接用于Gaussianized属性X64/X65或两个不同区间宽度的native属性。付费脚本除去虚属性宽度.15后才复核；如果换Gaussianized坐标，还须除CDF导数。
3. **补偿是局部材料等式。** D/Q路径需Q>0且新D仍在合法属性box内部；边界或额外约束可能截断此路径。它证明相应3D均匀子模型的gauge，不证明Full66相关field源分层拥有相同正质量交叠。
4. **host端不能照搬pure端切空间。** Q→0时属性导数随Q坍缩，需独立核验主物理系数不为零；本次只复用pure侧向量，未替host端作新物理证明。
5. “这不产生新的材料fiber”措辞可能被误解；建议明确写“这是原饱和decoder已有的同材料gauge，不是Maxwell传播制造的新不辨识性”。
6. c在付费脚本中由同一J线性组合构造；投影近零是公式一致性，不是额外独立验证了有限h下的Maxwell三次渐近。最终科学判断仍由主Codex综合A实际误差与范围。
