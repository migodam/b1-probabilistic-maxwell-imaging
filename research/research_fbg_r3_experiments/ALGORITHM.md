# 可运行算法与实际权限

本轮代码研究三种不同操作：成熟 CIE/TSOM 的联合电流优化、N1 的受控 current release、N2 的精确小支持多根。三者共享冻结 R4 的原始 scalar-TM 物理及同一原数据评分；它们不是同一个算法的三个名字，也不是后验采样。

## 输入和输出

公开输入：复 G/E0/S、匿名 Y、网格与采集几何、合法材料 box、公开初值、共同 smooth/TV先验及预算。全部材料为共享所有源的实数。4×4、6×6物理重新组装；576问题全部像素未知，不使用数字类别、支撑掩模或模板。

输出：一个原物理分数最低的在线 m，已付费且经原物理求解的材料候选、辅助/物理轨迹、停止原因和成本。候选可能是非零残差材料、不同初值端点或噪声相容点；不是默认的不同精确根。离线best@K、NRMSE与IoU只在预测封存后由评分器计算。

## 1. 原物理与普通 corrector

`src/physics.py` 是冻结 `research_z_alg_r4_diversity/src/physics.py` 的薄接口，保留 LU、多源 RHS、解析 receiver-transpose Jacobian、source adjoint和原材料正则。数据pack固定为source→receiver→real/imag。

```python
from src.physics import load_model
from src.full_maxwell_corrector import corrector
model, public = load_model('data/e1_public/obs_....npz')
answer = corrector(model, public['Y'], m0,
                   {'method':'trf', 'box':(.02,20), 'max_nfev':60})
```

普通对照含原生/倒数 TRF、原数据伴随 L-BFGS-B、GlobalDual及固定β=.2反馈Gram。不是有限差分Jacobian。共同corrector是成熟有界SciPy TRF；没有为了新方法削弱它。两个早期QP主动面加法舍入导致的TRF初值越界异常保留，后续共同接口显式投影到同一box再交接，未把旧失败改成成功。

## 2. 成熟结构 CIE-I / Fourier-TSOM

入口 `src/cie_tsom_baseline.py::run_cie(model,Y,m0,config)`。

1. 对公开S做薄SVD，按固定相对阈值1e-3确定DPCS，构造数据决定电流j_D。
2. 通过接收可见空间正交补投影Fourier列；mf是频域四角的边长，非模式总数。保留原生Fourier系数坐标和阶段零填充运输。
3. 每次做源分列PR-CG电流步及精确二次线搜索，再做共享实材料块。无共同正则时为合法modified-contrast区间LS；有TV时用相同材料先验的解析梯度有界内层优化。
4. 单β采用共同CIE目标；双β把同一物理m分别映射到当前块的rβ，不能原样跨β携带r。
5. 继承阶段并减小β/扩充Fourier空间。每个阶段末端重新解原始Maxwell，再选择在线物理incumbent。

小网格固定日程：4网格mf1→2、β1→.2；6网格mf1→2→3、β1→.2→.2。严格主对照采用系数平方相对变化1e-8、每阶段至少50轮、最多700轮及共同wall预算。原论文型1e-3停止只作为独立校准对照：在公开DEV上22/23轮会提前停止，而相同日程严格停止305轮把完整物理残差降到约4.27e-4/2.66e-4。没有依据早停原型宣称CIE弱。

本实现是从2019原始全文核实的二维实材料适配，保留DPCS/APCS、Fourier空间、交替、反馈日程和双β。CPU用显式IDFT列及缓存矩阵作用，**没有宣称复现三维论文图像或CG-FFT的速度复杂度**。CIE2011/2016原全文本轮未独立重取、一般相关噪声未实现，见独立来源审计。

## 3. N1：相同容量的选模、替换与软释放

入口 `src/controlled_release.py::run_release(model,Y,m0,config)`。

```text
solve exact physical anchor j*,Q*,L*
construct common receiver/Fourier/Krylov/random dictionary
form exact realified joint residual b, material B and current A_U
profile current -> H,g,p
for each equal-real-DOF candidate group:
    residual QR -> D,v,eta
    compute hard/soft/replacement material step
    rank proposal; pay full Maxwell JVP for two finalists
    retain residual pullback + curvature + common-prior model terms
execute local joint trajectory or two signed finite-current trajectories
read back physical materials; run common full-material TRF
select only by original full-data score + common prior
```

完整材料导数保留CIE的DP(m)·e。QR坐标软惩罚为γ||R a_C||²，而非裸系数ridge或非零列缩放。局部QP考虑实际box；理论步界只适用其明确的无约束/固定活动面条件，生产的有界数值QP不冒称有精确KKT证书。

所有N1选择器使用同一候选池、同一实际实DOF与获模开销。p=6时一个空间复模式对应12个实系数。修正后的公开可见backbone使4网格core156→active168实DOF，6网格204→216。成熟CIE的自由APCS容量与之不同，故“同DOF”主张限于N1选择/释放对照，不能直接扩写为与所有CIE版本容量完全相等。

初始backbone只给3个空间方向，存在严重辅助数据表示地板；该实现失败完整保留。第一项结构修正用公开S的可靠谱补足receiver backbone。该修正没有用GT，不能把它自身说成A17选模收益。随后又识别到可见/暗方向的分配与原生TSOM不同，故另建下述native-affine合同；两项修正各有独立源码与结果，未覆盖旧记录。

N1先比较CSI状态权重，再用固定β=.2的CIE状态权重作已声明消融，其他参数/候选池不变。两种状态权重都有真正的独立结果；没有只测试CSI就否定CIE法向理论。

### Native-affine N1（第二种独立表示）

`affine_mode='native_tsom'`固定公开数据决定的DPCS电流jD，以15个Fourier-APCS core加1个候选群构造同秩替换；6×6、6源时实际自由度为192实系数，与固定mf2 CIE匹配。所有selectors共享pool、原生source-normalized CIE目标、β=.2、prior与corrector。较强的多阶段CIE保留为额外参考，不能只用单mf弱化它。

联合轨迹从非物理jD起步，真实物理anchor仅用于原物理JVP、normal参考、评分与incumbent。DP·e不可删除，物理anchor的0/2/4阶定理也不能直接套在jD。内部每源白化与最终原数据/||Y||评分分开。所有fixed/task对照允许短期辅助目标下降而原物理预测收益为负，保留原物理最佳点；否则只会把CIE探索权限给新方法。

660次受控开发的四种有限current控制均从6/12增加到7/12联合恢复，task/normal/soft本身没有优于普通Fourier。该结果是generic finite exploration而非A17独有收益；另行冻结新的确认集检验，不能以此直接放开576维优越性门。

## 4. 有限法向探索

选择残差条件方向v+Dp，退化时用规定的局部特征方向。由当前残差QR坐标恢复真实current动作；以请求RMS上限2和current/||E||、defect/||E||上限3约束正负seed。低初值下，不能以很小的||j_anchor||作为唯一尺度而人为拒绝所有非局部移动。

若请求尺度不可行，最多8次显式可行尺度试探，两种符号共用同一个|τ|，保留实际RMS、CAP_LIMITED及所有拒绝成本；未假定非线性可行性全局单调。先固定candidate物理动作3步，再软释放最多80步，最后返回材料corrector。Fourier、随机、Krylov与直接材料±seed对照得到相同候选上限；实际seed/校正次数另报，不把配额当实际执行数。

本轮没有运行原生A17 Galerkin exchange来免费预测全部Full-J误差。N1材料相关评分受A17双通道思想启发；名称不构成A17必要性或原创性证据。

## 5. N2：小支持上的条件Q6/Q20

入口 `src/algebraic_multiroot.py::run_multibranch`。按公开2×2/2×3支持、全数据梯度与固定循环配额选择小块；支付k个完整RHS建立精确有限核。从既有源/接收选择少量通道，调用冻结Z-ALG商代数。共轭lift、实材料、box、节点与极点过滤后，用全部观测评价。

强对照获得相同精确小核：block-TRF、同核多起点、随机材料提案及普通全材料多起点。配额最多3个完整corrector，另有匹配**实际**代数校正seed数的对照。outside只在提案阶段固定，完整corrector释放全部材料。

数值解析的有限图册不是严格全根证书；子观测根不是完整观测解。没有使用另一端GT或oracle support初始化。

## 6. 预算与成本

每次完整LU、source/receiver/adjoint/tangent RHS使用冻结物理计数。`src/physical_budget.py` 在原operation执行前预付LU/RHS，停止时保存已付费原物理incumbent；不是另一个物理后端。共同wall在操作边界检查，实际wall跨界单列，不隐藏为0。Green/current模式作用、QR/SVD、内层材料求解与小核成本另记；相同LU/RHS上限不等于全部FLOP相等。

两个协议：每起点4s，K=1…5；固定总名义4s分给5个起点（每个.8s，同时每起点120LU/2400完整复RHS上限）。新鲜确认每起点4s、同物理上限、K=1…3。这些是有限预算比较，不外推K8/K16或全局收敛。

## 7. 576维 CUDA 成熟基线诊断

N1的有限开发增量与一个新鲜clean材料恢复尚不足以建立跨对象、跨噪声的独立优越性，N2也没有全数据增量，所以本轮不把两者扩大到576后宣称成功。独立运行的576实验用于补齐成熟CIE与固定反馈的强基线差距，采用旧R2已披露12个高对比对象。它不是新鲜盲测。

`src/cie_tsom_cuda.py` 复用现有SharedOperators/ObjectiveLane和一个CUDA上下文；公共receiver/Fourier基与G/S作用一次缓存，电流PR-CG留在各自stream，材料块复用CPU同一有界读出。每次返回的等待、传输及CPU材料成本如实计费。两lane有限FIFO持续补任务，不安装服务、驱动或新调度器，也没有宣称CUDA stream能够固定分配SM。

当前CIE576日程mf3→6→10、β1→.2→.2，45s总包络前55%电流proposal、余下共同原数据L-BFGS corrector；两个阶段都在完整576阶物理中验收。Native和固定β沿原R2历史强实现。所有图保持同一材料尺度；旧PIL数字按原始raster row0在上的约定显示，数组、物理和评分不作翻转。
