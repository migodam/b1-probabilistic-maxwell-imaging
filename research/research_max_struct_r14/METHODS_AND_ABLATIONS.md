# 方法与有限消融

## 原生四组合

`src/dualmetric_core.py`给原FFT-TSOM/global PR循环加完整state重入和正确的P1 current/adjoint。所有四臂从相同公开材料和冷DPCS系数出发。stage为mf3/6/10/12，最多500/500/500/1500步；首三stage不是按GT选择，最终stage达到预算时允许截断。每个材料对象的四臂在同一已验证P物理核上顺序执行，单BLAS线程。

Material为P0时保留原ordinary readout；P1时保留R7 Möbius/Arc两sweep及原prior接受。两者都使用实际材料相关变换的完整导数。P1的Green只在原disjoint2×2四cell内一次传播，self分母仍单独存在。01是原global A+-Arc更新，不把换名当作新方法。

## Shared material solver

0002/0017/0032/0037×四pair，material统一当前χ起点、同L-BFGS-B60/ftol1e−12/gtol1e−8/maxls25、box/prior、原m非增保护、无seed/Arc/MM。Current、PR、stage与物理corrector仍同一代码对象。该控制排除material求解程序不同的解释；不是独立确认，也不提供Mac对Windows速度结论。

0017的01仍恢复，其他三组合未恢复；0032四臂均好，0002/0037四臂仍差。因此至少这个关键正例不是Möbius/Arc的专属效果，也不是统一非局部metric自然拥有的效果。

## Self、β与次序

固定曝光16对象加self-only01/self-only11以及原global-Dual。Self-only只令Goff=0，保留Gii与β。完整01恢复8/16，self01和global-Dual各6/16，完整11与self11各5/16。0017/0045对非对角传播有决定性依赖；0024的global-Dual同样好、甚至误差更低，不能把全部β/读出收益记成空间传播。

单一候选01→11的schedule为[01,01,01,11]：仅原满Fourier入口切一次，利用原stage的PR reset，没有新增触发阈值。反向[11,11,11,01]只作曝光次序控制。前者保留01的8个联合达标对象，后者与11保留同5个。它支持早期可达路径重要的解释；不证明早期01对所有物理分布都好。

## 完整state干预

固定0002@128/512、0017@128/512/2500、0037@128/512共7state，复制后运行四pair×保留/清空PR，观察1/8/32步。J/Q/m、Fourier系数、历史gradient/direction、stage与权重都保存。干预比较固定更新次数，不把其中新增加的P1动作当免费，也不宣称等CPU。

干预包含完整原Φ、D/R、辅助port/构成缺陷、活动box、材料step与前一负物理梯度的角度。168个非初始角度有真实来源；初始没有step的56个角度保持缺失。有限差分保留正负扰动、两个尺度、分支与加法缺陷，不把单个JVP当成已证明的全程序Jacobian。

## 根与非根诊断

8²物理根实验在无prior、内点、独立current空间计算理想闭环；另跑真实一次有限程序、材料单独扰动和完整兼容扰动。后者明确揭示原readout guard的有限尺度deadband，不能由32步直接宣称渐近稳定。

非根material梯度分解为构成力、metric导数力和prior，直接复用原full derivative。该诊断只用原7×8 first-step state，不跑新forward，不作为新在线控制或物理信息来源。

## 原物理终点探针

五严重失败对象和两个预定控制各四pair，共28终点；保存原physical gradient、box投影梯度、沿该方向的retained-LU tangent、GN/先验曲率以及三种固定小步的实际Φ。所有试验材料仅在诊断目录，不插回在线候选池。它区分已知更好可行参考、当前点可微小下降、优化器容差/预算终止；不能由这组小步证明或排除严格局部极小。

## 新独立比较与成本

24新材料的11方法固定为00/01/10/11/01→11/CIE/Dual-CIE/global-Dual/G19/DR/direct。全方法同对象在同一个P或E物理核顺序执行，最多12 lanes（8P+4E），每进程单BLAS线程、同优先级。没有新CUDA移植或更换精度；复用FFT和算子，使预算投入科学比较而非新后端开发。

预算为整进程120CPU（proposal至72，之后余额可correct）；direct整120。完整费用包括导入、DPCS/SVD、重建辅助表示、固定投影、物理审计、失败、超时和导出。逐候选native CPU和外部OS child全生命CPU分开，最终以实际费用报告。

冻结的公开portfolio是01→G19→DR的已保存uniform子运行回放：当前incumbent球内则停，否则付费调用下一臂，并仅在Φ更低时替换。实际费用为全部调用子运行之和；它不是120秒总算法、不是warm restart、也不是离线GT并集。新确认中任何收益都需同时看这个费用差异。
