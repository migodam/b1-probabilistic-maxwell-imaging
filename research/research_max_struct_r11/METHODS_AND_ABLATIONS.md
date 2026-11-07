# 方法、消融与名称边界

Full G19是本项目对R10 Green分母材料读出的实验命名，不是独立标准文献算法，也不是原版A/A+的同义词。它仍建立在CIE/FFT-TSOM电流表示和相同原物理corrector上；本轮分别记账原G19、冻结A、contractedT等分母规则，不能把这些版本的全部结果直接归给原A+。

**fullGreen/fullG19中的full指保留材料依赖分母的完整导数，不指采用全部Green自由度。** 二维K为固定19维Galerkin截断，三维K为公开source空间rank120的P G P。全文把K=G只用于明确标注的解析机制极限。

## 用公式读懂Full G19

在固定辅助电流J的材料读出中，Q=Einc+GJ是该电流产生的总场，e=J−D_mQ是本构缺陷。它只有在e=0时才与当前材料的物理解一致。普通CIE使用

`w0=(I+βD_m)⁻¹e`，

Full G19将其改为

`w=(I+βD_m−D_mK19)⁻¹e`。

普通分母逐cell缩放缺陷；新增K19使不同位置的缺陷经过Green耦合后共同影响材料读出。分母可逆时，w=0与e=0等价，因此没有替换Maxwell本构方程；改变发生在尚未满足物理关系时的辅助搜索路径。非零残差的大小、方向及其材料导数都会改变，这可以帮助避开一个坏区域，也可以把另一对象推向坏区域。

“Full”具体保留的是`dw[h]=−A⁻¹D_h[Q+(βI−K19)w]`中的全部分母导数。只看`−A⁻¹D_hQ`便是另一种冻结分母的更新，不再是这里的Full G19。19也不表示材料只剩19个未知数：二维依然反演全部576个实材料值，原FFT-TSOM current容量、完整G正演及公共corrector均按共同合同运行。

因此，图中G19的收益属于“在CIE/FFT-TSOM上改变材料读出”的结果。是否来自Green分母，要与使用相同global current的Global-Dual比较；是否值得用，还要与主要基线CIE和Dual-CIE比较实际费用。Global-Dual是定位原因的消融控制，不是必须取代这两个主要比较对象的第三套论文基线。

## 二维实际比较

|方法|current块|材料块|角色|
|---|---|---|---|
|CIE-single|成熟sourcewise PR，βj=βm|普通CIE，原prior|主要成熟基线|
|Dual-CIE|成熟sourcewise PR，βj与βm分开|普通CIE，原prior|主要成熟基线|
|Global-Dual|global PR，同双β|普通CIE，原prior|隔离global/sourcewise与材料块|
|fullG19|与Global-Dual相同|A(m)逆残差，完整材料梯度|R10已有新分母的继承控制|
|freezeA|与Global-Dual相同|每次readout固定整个A|删去分母变化作用的机制控制|
|contractedT|与Global-Dual相同|保留z收缩、固定T(zref)|只冻结非局部feedback的控制|
|fullGreenBothBlocks|global PR采用同Green残差族|fullG19完整导数|R11修复候选|

48对象覆盖前六臂；另16对象覆盖三个基线、fullG19与BothBlocks。算法没有因新确认GT变化。BothBlocks没有增大Fourier/APCS容量，也没有改变最终Maxwell数据目标、材料空间或prior。它有额外small LU/三角解费用，绝不视作免费。

普通材料块沿用项目成熟_baseline_readout：含prior时保持该实现的求解，而不是为了“解析读出”删掉prior。新矩阵材料full/freezeA/contractedT保持冻结8步inner预算；不能将更准确的材料子问题直接等同更好的全波材料恢复。两个β不同意味着这些算法一般没有共同交替下降目标。

Green-current/ordinary-material这一缺失的2×2因子消融已在三个公开对象完成：环/壳/桥NRMSE=.4220/.5151/.5546，均未恢复。两块Green在旧桥的修复因此不能归因于current单独改善；完整四格表见delegated/diagnostics/current_only_dev/REVIEW_ZH.md。这是开发诊断，不改变已冻结确认组。

## 三维实际比较

sourceCIE、sourceDual、GlobalDual直接复用成熟三分量current/共享各向同性材料实现。新增full、lagged、contractedLagged在DEV中比较；新确认保留full与contractedLagged两臂。三维current仍是普通CIE，不把二维BothBlocks修复称为已经三维复现。

3D材料m为N实变量，场/电流为3N复分量；梯度先按三个分量与所有源求和，再投影同一实box。每次readout8次投影方向步；full以方向Gauss–Newton曲率初始化Armijo线搜索，lagged与contracted的固定二次数据项作可行方向精确一维最小化。R=0，因此不存在漏掉二维prior的伪推广；当前不声称适用于复介电或未知张量材料。

所有候选最后交给相同原物理corrector。三维原始final选最低已付费corrector Phi；成本前缀另在全部已付费initial/stage/corrector状态中按同一Phi选择。两者标签不同，不能偷偷用GT选择哪个状态更好。

## 为什么没有宣称全局凸化

R9已证明可逆Möbius坐标与严格递增loss变换保留局部极小，相关结论直接继承，不重复验证。R11固定A的单次材料块在凸prior下为凸；固定T的数据部分对z凸，但原R(m(z))未必凸。两种更新都不是moving目标的MM，也没有让Maxwell inverse problem成为凸问题。

Green修改确实可以改变辅助优化路径；它同时可能帮助和损害恢复。R11的科学问题是这个修改是否在独立材料、实际成本和真实矢量物理下值得保留，不是为某个变量变换冠以新的全局算法名称。

## 复用与新实现

- 普通基线、FFT/APCS、original-Maxwell corrector、2D评分：继承R8/R10；旧目录不编辑。
- R11新增material_readouts.py：full接口、freezeA、contractedT。
- R11新增cie_feedback_current.py：独立复制并只改变Green-current H/T、真伴随、方向二次线搜索与匹配日志。
- R11新增vector_denominator.py：真实3N共享N材料的低秩分母与三种readout；vector_cie.py保留ordinary global current。
- vector_cuda_dense.py只去掉每个新材料的empty_cache；库选择探针在公开材料上付费、只做一次，保持FP64和相同正演。
- NN、训练集、增广材料prior、新测量、扩大测试支持先验：本轮未使用。
