# 复杂场景：代数结构应放在 latent，还是 characteristic subspace？

这是响应本轮执行过程中用户追问的理论接口。**本文件仅推导，未实现 Full66 或新的 latent 成像实验，不计作 R3 的已验证性能。**

代数关系来自材料、本构、传播与观测。latent 可以规定允许的材料结构；固定传播或观测子空间可以承载消元及计算复用。两者不是互斥选择，也不是把同一个小型商环随意换一套坐标。

## L1：最简单的反例——一个 latent 不等于一次材料更新

若一个 latent 控制全部单元，$m_i=z$，则

$$M(z)=I-zG,\qquad\det M(z)=\det(I-zG),$$

其次数一般为 $N$。原来关于每个独立 $m_i$ 的一次性，并没有变成关于 $z$ 的一次性。一个标量变化 $h$ 对应 $D_{m(z+h)}-D_{m(z)}=hI$，非零时秩是 $N$，不是1。

更一般地，若 $m=m_0+Bz$，则 $\mathrm{rank}\,\mathrm{diag}(Bh)$ 是 $Bh$ 中非零分量数，通常仍为 $N$。因此小 latent 维数不自动给出同阶 SMW 更新，也不允许把 Q6 的四个实际单元换成四个稠密 latent。非线性 decoder 还可能完全破坏有限多项式结构；若 decoder 是已知多项式，当然可以重新研究其次数和稀疏性，但必须重做分析。

**成立的特殊情况：** latent 直接对应少数独立材料单元时，原多重仿射结构保留。若一个参数共享于多个单元，次数至多为它控制的行数，而非一。固定 support 的精确 ROI 消元仍可用，收益由未知 support 大小决定。

## L2：latent 可以直接进入全观测兼容性求解

设固定网格材料 decoder 为 $m=\mathcal D(z)$，固定参考 $m_0$，$\Delta(z)=\mathrm{diag}(\mathcal D(z)-m_0)$。全部观测已通过

$$W=W_0+CZ,\qquad Wq=T$$

精确消元。直接在 latent 中求

$$\psi(z)=\min_C\|W-(L+WK)\Delta(z)\|_F^2$$

即可保证每个被提出的材料都来自 decoder，不需要先生成任意像素图再投影回 latent。

假设内层最小二乘满秩、decoder 在当前点可微，令最优残差 $E=W-(L+WK)\Delta$、$V=L+WK$。对实材料 decoder，

$$g_i=-2\Re\sum_r\overline{E_{ri}}V_{ri},\qquad\nabla_z\psi=J_{\mathcal D}(z)^Tg.$$

若 decoder 输出复材料而 latent 仍为实，直接用

$$\frac{\partial\psi}{\partial z_a}=-2\Re\sum_{r,i}\overline{E_{ri}}V_{ri}\frac{\partial\mathcal D_i}{\partial z_a}.$$

**证明：** 内层正交性消去 $dC_*$ 项，再对显式本构材料项用链式法则。它不需要 decoder 为多项式，也不需要构造一个 latent 商环。饱和 seam、不可微 decoder 或变化内层秩需分别处理；先验约束和正则项必须明确加入。

**计算含义：** 外层未知量从 $N$ 变成 latent 维数。只有事先保证 $\mathcal D(\mathcal Z_{admissible})\subset\mathcal M_{physical}$ 时，decoder 才同时保证 box、被动性等材料合法条件；仅仅属于某个 decoder 的像并不足够。内层矩阵仍由实际物理 support 决定，不能把 $N\times(N-s)$ QR 算作 latent 尺寸的求解。原始 Maxwell likelihood 仍是最终排序目标，兼容性目标不是新先验密度。

## L3：固定 characteristic subspace 适合消元，但截断会改变问题

如果这里的 characteristic subspace 是固定 Green 或接收算子的 SVD 空间，需分两种情况。

1. **精确仿射数据消元。** 完整 $Zq=0$ 保留所有满足数据的 $W$。这是等价表示，真正省掉观测约束的反复优化；几何和采集不变时可复用。
2. **截断候选空间。** 若只保留 $Z_k$ 的部分行、仍保证 $Z_kq=0$，则 $W=W_0+C_kZ_k$ 继续精确满足观测，但可表达的物理传递矩阵变少。

对嵌套行空间 $\mathrm{row}(Z_k)\subset\mathrm{row}(Z_{k+1})\subset\mathrm{row}(Z)$，固定材料处有

$$\psi_k(m)\ge\psi_{k+1}(m)\ge\psi(m).$$

证明是最小化可行集合逐步扩大。完整物理真根也可能在截断空间具有正兼容性误差；因此不能把截断后的失败解释成没有对应材料。零代数残差依然能在可逆条件下验证物理解，但完整性可能丢失。

同理，材料对角算子一般不保持固定传播子空间：$[D_m,G]_{ij}=(m_i-m_j)G_{ij}$。固定 SVD 有助于投影、提案和近似预条件，不能无条件把非线性 resolvent 永久压缩为同一个小系统。若做近似，需要 R9 缺陷校正与实际残差规则。

SOM/TSOM 的多重作用可以用投影和投影像表达；两个投影的复合一般不是交集上的正交投影。应保留具体作用顺序及其物理意义，而不是用一个难以计算的“空间交集”替代实际算法。

## L4：适合现有 P3L/Full66 的下一步接口

```text
固定网格与采集：G、S、Einc，按合法未知 support 作精确 ROI 消元
材料模型：沿用原 decoder D(z)、其 Jacobian 和合法 prior/box
数据消元：固定传播/测量子空间，形成 W0、Z、ZK
候选生成：在 z 中优化共同材料兼容性；每个候选已经属于 decoder 像
物理修正：对 F(D(z)) 使用原全观测 Jacobian/adjoint
选择：预声明 likelihood/结构目标，独立 full Maxwell 验收
评估：冻结以后才计算 GT top1/best@K
```

若形状 latent 只改变固定网格中的材料占据，固定 Green 假设可以保留；若 latent 实际移动网格、边界、源接收或改变背景 Green，则 $K,q,L,Z$ 也变化，上述梯度必须包含它们的导数，不能继续声称全固定复用。

本轮没有运行这个 Full66 接口。最短的下一轮实验是：在原冻结 decoder 和相同观测/初始化下，对照原 latent Full-J 与上述 latent 兼容性提案＋原物理修正，保留完整内层 QR 与 decoder 费用。只有这一步能回答复杂场景中的实用收益，当前12单元成功不能替代它。
