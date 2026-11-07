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

## 本轮真正新增的实现

1. 同一两级current Cascade的完整、528接收核、557内部弱、509联合暗和软释放层；保留完整电流头而非任意小型weak latent。
2. 固定锚点归一化的本构／内部场损失，及条件启动的精确VIE微调与相同更新次数对照。
3. 固定投影中真实Maxwell电流的before/raw/final绝对能量、更新能量与带符号交叉项。
4. 新对象、未见类别、未见字体/规则三类测试；逐图标记训练身份。

这些是受控算法实验的新增组合与诊断，不等于发明投影、CNN、Cascade、CSI本构残差或伴随微分。论文归因必须由实际运行及强基线决定。

## 当前冻结源码与数据

实际来源base commit见`audit/PROTOCOL.json`，CIE逐文件复制来源见`vendor/PROVENANCE.json`。不覆盖冻结R2/R3/R4/R9源码；本目录只通过新适配器调用。来源记录采用Git版本与实际文件清单，不进行SHA-256检验。
