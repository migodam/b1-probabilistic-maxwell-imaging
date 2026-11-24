# 投稿稿 v2 修订记录

在原生 96² 网格、全环、2% 噪声和材料上限 0.75 下，两件物体的四区单源半定松弛加入严格电流能量帽后，对任意单位相干照明仍有精确零残差解。多源常权下界分别为 0.16070 和 0.11061，高于同一噪声球半径 0.09097 和 0.08803，仍能排除材料峰值不超过 0.75。

| 位置 | 修订内容 |
|---|---|
| §3.4 | 逐像素的扇区矩阵不等式与光学等式等价于实际电流来自同一个实对角材料。证明覆盖总场全暗的像素和上限为零的像素。 |
| §3.4 | 用约束嵌套串联全局、四区、十六区与逐像素极限。逐像素等价适用于实际电流及其秩一矩，半定松弛中的高秩矩仍可满足逐像素约束。 |
| §5.4 | 全盒电流能量上界为 2,340,476，覆盖给定十束入射波的任一单位相干组合。使用读数与噪声半径的上界分别为 4,982,633 和 1,035,001。两组精确拟合读数的可行矩分别只需低于 8,500 和 7,000 的能量。同一材料盒内有一件材料的实际电流能量超过 8,886.62，因此任意对全盒及全部单位相干波都有效的统一能量帽，都允许这两组矩。全部指定扫描点的最优残差都是零。 |
| §6 | 能量约束下的分离限于四区单源半定松弛。采用冻结的十六区权重，单个相干方向也能排除同一材料上限。 |
| 图表与引文 | 六张存档图随稿保留。34 条参考文献的已有出版核查随稿附上。 |

[能量上界的推导、区间验证与扫描证据](ENERGY_EVIDENCE_ZH.md)

## 送审文件

[匿名英文稿](manuscript_ip.pdf) · [匿名排稿源码包](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/optical-inertia-ip-v2/IP_Anonymous_LaTeX.zip) · [英文投稿信](cover.html) · [逐条引文核查](references.html) · [上一版](../ip_final/index.html)

*Inverse Problems* 允许作者选择双匿名评审。匿名排稿包只含编译论文的文件。论文在致谢和实验方法中披露 AI 工具用途与验证方式。送审前请由作者本人核对引文原文和 AI 辅助内容。 [期刊评审方式](https://publishingsupport.iopscience.iop.org/journals/inverse-problems/about-inverse-problems/) · [匿名核对指引](https://publishingsupport.iopscience.iop.org/questions/checklist-for-anonymising-your-manuscript/) · [作者 AI 工具指引](https://publishingsupport.iopscience.iop.org/questions/generative-ai-tools/)

## 排版核对

终稿 26 页，保留 34 条参考文献和六张图。匿名源码包可独立编译，排版日志没有未定义引用或越界警告。
