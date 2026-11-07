# GN–NN R3：给下一轮理论分析的交接

## 最值得解释的实验证据

不是“网络能不能输出一张漂亮图”，而是以下四个同时出现的事实：

1. clean/noisy的GN–GT电流误差平均98.09%/91.85%位于接收零空间。
2. 自然GT电流零空间投影的材料NRMSE为0.00846/0.25706；后者64个对象全部无法通过物理数据门。
3. J-null的noisy材料误差在三种子为0.2933/0.3215/0.2720，J-full为0.1569/0.1652/0.1266。硬零空间约束在这组噪声条件下没有独立收益。
4. 监督投影GT电流比材料损失更准确地预测电流，却没有更准确地恢复材料。J-full可以输出偏离真实电流很远的辅助变量，读出材料后反而较好。

这要求把**电流表示、材料读出、物理本构相容和学习优化几何**拆开。当前数据不支持“零空间维度很大，所以更容易学到GT”的推论。

## 条件数学结构

令 $L=I-D_mG$、$d_i=\|Q_i\|^2$，在真实锚点、未裁剪内点共享实材料读出：

$$
\mathcal R_m[V]_i=\operatorname{Re}\langle Q_i,(LV)_i\rangle/d_i,
\qquad \mathcal R_m^*h=L^*D_{h/d}Q.
$$

$$
C=\mathcal R_mP_0\mathcal R_m^*
=D_{1/d}\operatorname{Re}[(LP_0L^*)\odot\overline{QQ^*}]D_{1/d}.
$$

它描述**允许非物理辅助电流时的局部材料读出**。C不是材料Jacobian信息矩阵，也不是物理解纤维的维数证书。

若 $P_0E=P_0(G-G^*)=0$，则 $P_0L^*Q=0$，C至少缺失d方向。冻结矩阵的有限精度诊断在1/3/6/12源和两个公共锚点给出秩575。TEST的GN–GT材料误差沿d方向能量均值只有clean1.888e-7、noisy2.289e-5。主要差距不能由这一维缺失解释。应研究C的非零谱、条件数与材料损失训练的pullback，而不是回到一个任意Weak64截断。

有限偏差也可精确写出：$j_a=j_*+V$、$q_a=Q_*+GV$时，裁剪前

$$
m_a-m_* = \operatorname{Re}\langle q_a,L_*V\rangle_{row}/\|q_a\|_{row}^2.
$$

对自然oracle，$V=-V_{48}V_{48}^*(j_*-j_0)$。该式提示：电流L2投影不自动优化材料误差；带噪GN可观测分量的误差经过L和总场权重，会在材料读出中放大。尚未完成的是一个可计算、对实际学习有效的条件界或预条件设计。

## 三个应当优先裁决的假说

**H1：硬接收锁定过强。** 研究 $Z=P_0Z_0+V_{48}a$，给a一个与噪声/数据门匹配的有界或协方差加权限制。J-full的正面结果支持做软可观测释放消融，但没有证明其收益就是来自放开可观测部分；还可能来自读出条件与网络梯度。必须与同网络J-full、J-null配对比较，按实际自由度和成本报告。

**H2：current loss选错了商坐标。** 多种辅助电流可以给出相近材料。欧氏电流误差可能在与材料无关的方向花掉学习能力。可以尝试材料相关pullback或本构法向加权损失，但必须保持最终原物理校正，且不把辅助电流误差变小当作材料成功。

**H3：V-physics的瓶颈是训练而非表示。** seed7的VAL/TEST好，17/27在VAL就不好；同样预算下J-full稳定得多。需要看输入支路梯度、归一化与特征利用，做少量明确消融。不能在已披露TEST上挑最好种子。当前两级voxelCascade只证明seed7下noisy NRMSE0.1157、23/64联合成功。

## 普通反演对照仍需补哪一刀

低对比带噪纯数据GN材料NRMSE1.0164，继续GN变1.3477，而固定弱正则45秒为0.2958；A=2.5 noisy弱正则甚至优于所有主NN。这是典型的噪声/正则问题证据。下一轮应首先冻结discrepancy stopping或有限正则路径。当前NN相对纯数据GN的大收益不能全部归为新搜索机制。

## 数据与评价的已知不足

- 数字与幅度循环耦合：1/4/7=.1，2/5/8=.2，3/6/9=.4。新确认必须独立交叉；现有64母体不可继续当新鲜测试。
- Otsu对连续GT自身仅47/64通过既定IoU门。原门未改；下轮应事前定义更适合连续材料的结构指标，同时保留RMSE与数据门。
- 三种子是训练重复，不是192独立材料；clean/noisy/stress也是配对重复。
- 575是数值诊断，不是一般Maxwell定理；自然oracle不是最优表示下界。
- 当前样本为同网格合成scalar-TM，不是3D矢量或实验室数据。

## 可直接使用的证据

[主报告](RESEARCH_REPORT_ZH.md)、[完整数学](MATHEMATICAL_THEORY.md)、[逐案例](results/RESULTS_BY_CASE.csv)、[机制链](results/MECHANISM_SUMMARY.csv)、[按母体配对区间](results/PAIRED_PARENT_COMPARISONS.csv)、[所有图](figures/all_test_cases/GALLERY.html)、[Release](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/gn-nn-completion-r3-20261107)。

最终建议：**继续完整voxel与完整current材料损失两条主线；把硬接收零空间作为受控消融，暂不作为下一轮唯一NN表示。** 本轮不是概率解集研究，不输出后验，也没有证明新的Maxwell分支。
