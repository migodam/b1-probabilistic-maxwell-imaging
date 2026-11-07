# 公式与代码复用范围

| 资产 | 本轮用途 | 条件与不允许的归因 |
|---|---|---|
| R4-Diversity `src/physics.py` 的Model/Regularizer/material_basis | G/E0/S、保留LU、所有源共享RHS、解析Full-J、8→12→24材料参数化 | 完整576阶scalar-TM VIE，非3D矢量Maxwell；粗材料阶段不是物理矩阵降秩 |
| GN–NN R3 `src/networks.py` / `train.py` | 同一小残差U-Net、公共物理通道、共享实材料LS、TRAIN统计与VAL选择 | 本轮新网络全部从头训练；没有加载R3已见4/7类别权重 |
| GN–NN R3 `src/evaluate.py` | 相同全材料有界TRF校正、候选槽、物理合格后贴近NN先验的在线top1 | 同一冻结离散VIE；不是独立物理模型或全局收敛证书 |
| GN–Green R1 `results/fixed_spectra.npz` | 固定G/S右奇异基、19/557及9/39/528分组 | 强阈值.1是声明的定义；GIN数值满秩，没有虚构内部零空间 |
| GN–Green R1 `MATHEMATICAL_DEFINITIONS.md` | Two-fold顺序投影、数值交集与未分类余项 | 顺序投影不等于交集；重叠投影不能归一化为共同互斥物理分支 |
| GN–Search R2 `src/search.py` | 成熟全像素校正器、预算与box | 不采用旧Weak64/128，不重跑GT拉回零核证明 |
| FBG-R3实验的成熟结构CIE/Fourier-TSOM | 独立传统参考，原源码复制到本轮vendor，保留来源commit/path | 2D实材料适配，显式Fourier列；不是原论文3D实验逐项复现，不声称FFT复杂度 |
| 现有CUDA算子目录及A9 SSH transport | 查看资源占用、复用有限远端任务和P核分工 | 4060被其他任务占用时使用MPS；不安装GPU服务、不抢占别的研究 |
| A17 `Gaussian/A17/CLOSURE_R1/manuscript/CORE_THEORY.tex` 式4–7 | 区分直接电流路径与 retained-space feedback，启发补偿后的物理反馈 | 本轮RNN没有实际调用A17 Schur块，不声称其条件交换证书已经移植 |
| SOM-Net，Liu等，IEEE TGRS2022，DOI10.1109/TGRS.2022.3217162 | 作为电流／材料交替学习展开的直接先例 | 不能将通用电流—材料RNN本身声称为原创；[原文](https://arxiv.org/abs/2209.03567) |

## 本轮真正新增的实现

1. 同一两级current Cascade的完整、528接收核、557内部弱、509联合暗和软释放层；保留完整电流头而非任意小型weak latent。
2. 固定锚点归一化的本构／内部场损失，及条件启动的精确VIE微调与相同更新次数对照。
3. 固定投影中真实Maxwell电流的before/raw/final绝对能量、更新能量与带符号交叉项。
4. 新对象、未见类别、未见字体/规则三类测试；逐图标记训练身份。
5. 用户追加的共享权重电流—材料—固定Green反馈补充，二次数据与噪声容许数据分开。实现和运行状态见补充合同与最终清单。

## 递归反馈公式索引

完整假设、证明与先例边界统一见 [A17递归反馈附录](theory/A17_RECURRENT_FEEDBACK.md)。以下都是有实际计算用途的公式，表中列明标准工具来源，不以重新编号宣称原创。

| 编号 | 已推导公式／条件 | 计算用途与实现 |
|---|---|---|
| F-RNN-01 | \(j(m_{read})-j_a=-L(m_{read})^{-1}e_a\)，要求L可逆 | 解释辅助接收零空间与真实材料场之间的差；真实物理验收及current_realization_gap表 |
| F-RNN-02 | 剖面梯度 \(S^*r/a+\rho L^*e/b_0\)，材料逐行LS唯一且归一化固定 | 以G/G*代替每层完整逆；标准包络定理，`feedback_step` |
| F-RNN-03 | 固定材料二次直线步长 \(\|g\|^2/(\|Sg\|^2/a+\rho\|Lg\|^2/b_0)\) | 每个样本一个标量步；固定材料下降、再次读出后的剖面下降；标准最速下降/majorization |
| F-RNN-04 | 物理解处 \(H_{prof}=S^*S/a+\rho\mathcal N^*\mathcal N/b_0\)，核为 \(\mathcal T\ker DF\) | 区分材料不可实现法向、可观测材料方向、中性方向；条件局部分析，不是全局GT收敛定理 |
| F-RNN-05 | 噪声球梯度因子 \(\theta=(1-\delta/s)_+\)，曲率上界仍用完整Sg项 | 避免纯二次数据反馈继续推动噪声内合格物理状态；标准平方距离梯度，`feedback_hinge` |
| F-RNN-06 | 内点物理锚点的恢复导数 \(D[j\circ\mathrm{Read}]=L^{-1}\Pi_Q^{\mathbb R}L\) | 幂等斜投影解释硬暗电流在物理恢复后重新进入可观测方向；链式法则推导，理论分析，不伪称网络实现了该完整逆 |
| F-RNN-07 | 固定材料、状态梯度迭代误差 \(j^K-j(m)=(I-\tau L^*L)^K(j^0-j(m))\) | 说明迭代可近似逆的作用及病态时的慢收缩；标准Landweber条件结果，不作为NN复合映射的收敛证明 |

这些是受控算法实验的新增组合与诊断，不等于发明投影、CNN、Cascade、CSI本构残差或伴随微分。论文归因必须由实际运行及强基线决定。

## 当前冻结源码与数据

实际来源base commit见`audit/PROTOCOL.json`，CIE逐文件复制来源见`vendor/PROVENANCE.json`。不覆盖冻结R2/R3/R4/R9源码；本目录只通过新适配器调用。来源记录采用Git版本与实际文件清单，不进行SHA-256检验。
