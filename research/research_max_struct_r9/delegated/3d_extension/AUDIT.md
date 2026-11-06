# R9 真三维扩展：入口、物理资格与严格基线审计

**状态：只读审阅完成，未实现、未运行、未生成新数据、未读私有GT。** R8保持原样。本文提供root选择科学合同的依据，不表示已经从2个案例增加到4个可信成功案例；52个R8比较臂也不是52种物体。

源资料来自实际 sibling `r8_worktree`。下文R8相对路径均以该仓库为根；R9路径以当前 `r9_worktree` 为根。原始2019合同已在 `research/research_max_struct_r8/delegated/paper_baselines/PAPER_CONTRACTS.json::CIE2019` 核对，主来源为 [Zhong & Xu 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC8320915/)。

## 1. 首要结论：可复用组件存在，但没有现成的严格3D A-material runner

最快可信的组件入口是：

- `research/research_max_struct_r8/src/vector_cie3d.py::VectorCIE3D`：R1全dyadic DDA、60有限电偶极入射、全矢量接收、原物理值/梯度、材料CIE/A读出、prior、已付费值缓存。
- `research/research_max_struct_r8/src/vector_pr3d.py::VectorPR3D`：sourcewise/global PR系数步，以及实各向同性A+坐标arc读出。
- 原坐标使用的实际调用示例：`experiments/vector_pr240_fourier.py`；但它写R8输出且自动离线读truth，**不能直接作为R9公开求解器启动或导入实验脚本复跑**。只复用类/纯组件，在R9另写public-only编排器并把scorer分离。
- 大电尺寸的已有传播组件：`src/vector_fft3d.py::VectorFFT3D`，复用A18_C padded卷积形式和R1自项；目前只有direct物理反演，尚无该后端上的CIE/APCS/A/A+集成。

### 阻断性混淆：R8的“A”同时改了current和material

`vector_pr3d.py:56–60` 为A/A+选择 `metric='A_TSOM'`，并把同一个metric传入 `current_pr` 和 `material`。`current_pr:13–16` 不仅用 `beta−s`，还将局部 `(I+diag(z)H)` 作用于整个current方程。因此R8结果是 **A-current + A-material 的联合适配实验**，不是R6/R9“普通global CIE current固定，只替换材料读出”的A-material延伸。R8的 `CIE_global` current beta也为6→1，而非所需的6→3。

R8全部旧分数/费用必须保留；它们只能作为已明确的混合适配证据，不可重命名为严格A-material结果。R9需要新的编排器把 `current_metric` 和 `material_metric` 分开，至少冻结以下六臂：

|R9严格臂|current residual|PR聚合|current beta|material residual/readout|material beta|
|---|---|---|---|---|---|
|ordinary CIE|`CIE_TSOM`|sourcewise|6→1|普通共享材料CIE|6→1|
|Dual CIE|`CIE_TSOM`|sourcewise|6→3|普通共享材料CIE|6→1|
|Global Dual|`CIE_TSOM`|global|6→3|普通共享材料CIE|6→1|
|A-material|**同Global Dual的普通CIE current**|global|**6→3**|`A_TSOM`材料目标|6→1|
|A+-material|**同Global Dual的普通CIE current**|global|**6→3**|相同A材料目标的合法arc提案|6→1|
|direct physics|无辅助current|不适用|不适用|原全Maxwell目标|不适用|

普通CIE与LSIE不是同一臂。若后续要声称复现原论文“LSIE失败、CIE成功”，应另保留明确LSIE对照，而不是将single CIE称为LSIE。是否添加由root决定，本文不追加任务。

相同current须共用相同DPCS/APCS仿射族、原Fourier系数度量、初始化、重启/transport规则及source metrics；只对比末端DOF数量不够。R8 `current_pr:18–20` 在mask切换时把系数置零。原2019全文§2.3 Step5明确把得到的contrast与APCS作为下一轮初始值；因此这不是仅文档差别。R9应复用成熟2D的嵌套mask精确transport思路，按(kx,ky,kz,component)键零填充原γ，保留Jd与投影器，允许重置CG历史但不把已有APCS清零。

此外，R8的ρ=10、quadratic smooth/l2 prior、bounded χ-space L-BFGS材料读出和固定48/240循环均是适配设置；原2019 §2.3为modified-contrast LS读出、跨轮APCS继承及sqrt(mean_source系数相对变化平方) <1e−3停止。R9可为公平共同prior保留有声明的适配，但不能把它叫原论文成熟实现。

## 2. 必须保留的矢量物理与A+适用范围

`research/research_max_struct_r1/src/vector_backend.py::TensorDDA` 使用cell-major三分量全电dyadic Green，包括纵向近场项。约定为exp(−iωt)、λ₀=1、k=2π、μr=1，电流密度j=χE、偶极p=Vj；体素自项为

`s = −1/3 + i k³ V/(6π)`, `G = V Goff + s I`, `alpha(χ)=Vχ/(1−sχ)`。

这是一套CM+RR点偶极有限离散，不是H(curl)兼容离散或连续误差证书。R8物理选择目标为 `Phi3D=mean_source(||Sj−Y||²/||Y||²)+prior`，回执残差为该source-relative平方均值的根；它不是R9 scalar的`.5*global-relative-residual²+TV`。新3D合同须把同一个Phi、prior和metric用于所有臂，不能跨2D/3D直接混用Phi数值或成功阈值。三分量不能改为三份互不耦合的scalar-TM。60源是三条圆上的有限距离电偶极，轴向分别为xy圈z、yz圈x、xz圈y；接收60位置×3分量，Y维度为180×60。保留入射和S的相同单位/体积因子。

实各向同性A+的条件结构可复用：每个独立材料变量只控制一个DDA体素，物理自项为共同标量sI₃，且局部H的自身3×3块为零。令b=β−s、zᵢ=χᵢ/(1+bχᵢ)，固定其他材料时，叠加所有源与3分量仍有 `R=C+D zᵢ`，故state目标为 `q|zᵢ|²+2Re(h zᵢ)+const`。现有A+两轮arc驻点/端点/旧值提案由完整相同prior进行接受，并非精确求解整个正则化条件最小值。

- **可用范围：独立实标量各向同性χ、共同sI₃、零Hᵢᵢ。**
- 复χ具有两个实自由度，原实arc不能直接使用；已有bounded Re/Im A材料梯度可复用。原论文有损coated cube的A+应标NOT_AVAILABLE，除非root另行推导/实现合法复材料读出。
- 不得把标量arc直接推广到各向异性tensor ε。
- 不得随意用“一个粗材料变量控制多个细DDA点”降低未知量后继续声称同一arc证明：组内H自耦合一般不再为零，条件残差可能出现更高阶项。这需要新数学审阅。

### χ8的盒约束不能只改一个bounds数组

现有 `vector_pr3d.py:44` 的A+端点和根可行性写死0/4；`:51` 的所有材料优化/物理corrector也写死[0,4]，`vector_cie3d.py::optimize`同样写死。只生成χ8目标而不改这些地方将制造必然恢复失败；只放宽A+又会制造不公平优势。

R9必须由一个公开 `material_box=(lo,hi)` 同时驱动：所有材料L-BFGS bounds、direct/corrector bounds、A+端点、根筛选、旧值合法性、初始化检查和offline合法性；共同prior及其梯度/条件接受也必须使用同一冻结定义。χ8需要hi≥8；可供root考虑公共类别盒[0,20]，但这不是已冻结配置。不能暗改prior强度来挽救GT分数。现有vector prior是N/edge归一离散smooth+l2，不是2D旧smooth/TV，也不是h⁻²连续H1；跨模型、跨分辨率须声明。

## 3. DPCS硬锁与公共noise-ball规则

`vector_cie3d.py:53–67` 的rank硬锁为 `sum(singular >= .1*singular[0])`，不是2D历史的.001；还默认 `whiten_apcs=True`，会把projected-IDFT列换成SVD正交range坐标。后者保留仿射family但改变PR系数度量。R9公平主比较应显式关闭白化，保留原projected-IDFT γ。

### 现有R8噪声不能直接冒充新规则的iid CN噪声

`experiments/vector_benchmark.py`、`vector_hollow_cube.py`、`vector_original_size.py` 都把复正态样本再次乘以 `.1*||Yclean||/||noise||`，强制每次噪声范数恰好为10%。这样得到的是固定能量的方向分布，不再是独立CN样本。原observations还保存clean_Y/noise，且没有独立public sigma。不能从clean_Y或已实现noise读取/估计sigma后，宣称套用了新的chi-square事件。

若使用R9已审阅的同一noise rule，必须在新R9生成合同中：

1. 公开一个有限、实、标量sigma，明确E|n|²=sigma²；生成 `sigma/sqrt(2)*(N(0,1)+iN(0,1))`，**不再按实现噪声范数重标定**。
2. optimizer只拿Y、sigma和已知acquisition/operators；clean/noise/material/reference单独留给生成器/scorer。新public loader应拒绝这些私有键。旧R8含clean_Y/noise的NPZ不能原封不动作为严格public-only输入。
3. 对raw S、raw Y做SVD，以同一完整噪声球半径选最小supported contiguous rank：

   `R² = sigma² * chi2.ppf(.99, 2*180*60)/2`, 即复观测10800、实自由度21600。

   三分量不能遗漏，不能误把nr代成60。固定同一个R²，不做每个r单独显著性检验或按剩余rank改变χ²自由度；同一个完整噪声事件逐点控制所有正交projector。
4. 使用现有R9 `src/baseline_adapter.py::noise_ball_rank` 的规则/数值支持门槛作为定义，复用一次receiver SVD的U/Σ/V*。新增3D prepare hook应据它生成V、Jd及同一投影后的Fourier基。不要先调用硬锁的base.prepare，再只修改rank字段；那会留下不匹配的Jd/B。
5. 若直接调用当前helper，它自己会做SVD；再调用base.prepare会发生第二次SVD，必须如实收费。更干净的R9改动是单独的post-SVD rank-selector接口，并保留历史固定.1作为真正3D控制。.001若另加，是新对照，不能叫旧3D设置。
6. 对非一致噪声、sigma数组、covariance/whitening/观测归一化等未支持标记明确报UNSUPPORTED，不能默默按每源/每分量归一化sigma。没有满足噪声球的数值supported rank应保留FAIL，不用GT或随意rank兜底。

物理/辅助目标中的source normalization不等于噪声whitening；新rank选择仍使用raw观测。新sigma规则也不提供模型失配预算、唯一性或材料恢复保证。

## 4. 四个有意义的高contrast几何合同：供root冻结，尚未实现或运行

现有R8只有四方块、包覆立方体、空心体+两杆这三个主要几何族；旧cuboid纠正成真正cube不能重复计为第四个独立几何族。推荐第四个使用既有box primitive组成异形双体。建议如下，**χ8数值是R9新高contrast改编，不是原2019参数**：

|几何族|可复用的实际定义|建议高contrast版本（root待裁决）|原论文/适用范围|
|---|---|---|---|
|异质四方块|`vector_cie3d.py::objects`，centers(−.14,.14,−.14)、(.14,−.14,−.14)、(−.14,−.14,.14)、(.14,.14,.14)，每边.16λ|χ=[2,4,6,8]，即ε=[3,5,7,9]|原2019为ε=[2,3,4,5]；必须另列原参数，不伪称精确复现|
|同心包覆立方体|同函数outer side.4λ/inner.2λ，解析box fractions|主实材料版本χouter=4、χinner=8，即`4Fouter+4Finner`；可做全六臂A+比较|原2019/R8是εouter=1.5+.3i、εinner=2+.8i；原有损版本应保留为单独原论文控制，A+不可用|
|真正空心立方壳+两杆|`vector_hollow_cube.py::cubic`，outer.30³λ/void.18³λ，中心(0,.05,0)，杆中心(±.2,−.2,0)，半宽(.035,.035,.18)|各实体χ=8；保留空腔和两根分离杆|原2019为ε2.5，具体几何数值仅Fig11未在正文给出；此几何是显式adaptation。薄壁/细杆在8³格点上未充分解析，须以独立前向细化定量报告，不能靠图像断言成功|
|异形双体：cube+rod|复用`box_fraction`和已有杆primitive；例如cube center(−.12,.06,0)、half(.09,.09,.09)，rod center(.15,−.08,0)、half(.035,.035,.18)，两者不接触且在.6λ盒内|cube χ6、rod χ8，产生非对称双体及矢量交叉耦合|这是R9新组合，现无已运行数据或专用runner；不得说它已在R8测试过|

四个主实材料合同可共同使用相同未知体素空间、material box、prior、source/receiver、noise rule、六臂和原物理预算。若root要求四个主对象中包含原有损coated cube，应明确只有三个对象可比较实A+，而不能将复case上的A+缺项伪造成失败或补用不合法arc。是否提高壁厚/杆宽或改grid，由root在生成新数据前决定；不得看GT后调整几何使之易恢复。

这些都是真正3D对象和独立体素未知量，不用少量已知形状参数代替完整反演。所有旧R8形状/分数已经暴露，适合作DEV；新增GT封存与prediction先封存后评分必须另做，不能叫旧对象盲确认。

## 5. 物理资格现状：哪些成立，哪些没有成立

已存在且可直接作为历史证据复用，不需要重跑冻结恒等式：

- 复材料非根导数：`RESULTS/vector3d/DERIVATIVES.json`，physical相对差3.83e−9，LSIE/CIE/A材料≤5.62e−10。它验证新增复适配的导数接线，不证明原论文基线成熟。
- 独立Mie几何修正：`experiments/vector_mie_geometry.py`、`RESULTS/vector3d/mie_geometry/RECEIPT.json`。r=.145λ、ε3、.4λ箱、固定24³子体素fraction，静电极化率保持MG规则 `χeff=3fχ/[3+(1−f)χ]`；h=.08/.06667/.05λ的远矢量误差6.378%→4.882%→3.875%。没有拟合或体积重标定，但也不是严格Maxwell cut-cell界面积分。
- 真实成像独立组装cross-grid discrepancy：8³→6³四方块15.161%，cuboid+杆10.069%，有损coated4.456%；10³→8³四方块10.607%。这里“独立”指不同网格/装配/求解，不是另一个独立连续PDE求解器。
- 同一个LU的60源RHS与同源adjoint是正常复用；不能把它称为独立物理验证。已付费场的代数残差也不能代替Mie/网格误差。

这些ε3的Mie百分比不属于.6λ成像h=.1/.075，也不属于新的χ8/ε9范围。若root选χ8，最高材料对比下的独立Mie或实际几何细化仍是必要新物理资格，且应在配置冻结后只做所需验证，不能重新做无关formula/error-floor扫参。当前沒有χ8连续精度证书。

## 6. 实测尺寸与费用；不把小模型成本冒充论文尺寸

|已有真实运行|未知材料/矢量维数|实际预算|CPU与内存证据|
|---|---|---|---|
|6³四方块，原Fourier γ|216 real /648 complex field|240 PR + max300 physical corrector|六臂17.47–40.84 CPU秒；整个六臂138.43 CPU秒/181.07 wall秒，peak RSS165.4MB。A/A+含上述current混淆，只能作资源参考|
|8³四方块/10³合成|512 real /1536 field；fine1000/3000|48 PR + max300 corrector；direct192+300|LSIE101.28、Dual114.63、A+115.50、direct121.52 CPU秒；350/350/351/495实际LU，source RHS21000/21000/21060/29700，adjoint RHS18060/18060/18120/29640。全四臂453.73 CPU秒/594.37 wall秒，cached阶段peak428.3MB|
|6³真正空心cube+杆|216/648|48 PR +300；direct192+300|10.89–14.22 CPU秒/臂；全四臂50.12 CPU秒/75.16 wall秒，peak445.1MB|
|原3λ，inverse12³ /synthetic16³，matrix-free FFT|1728/5184；fine4096/12288|direct20步+3 accepted corrector后截断，30步allowance未完成|379.59 CPU秒/395.87 wall秒，peak356.0MB，粗细观测差36.113%。只做direct，没有原尺寸CIE/A/A+比较|

证据文件分别为 `RESULTS/vector3d/pr240_fourier_coordinates/PR_COMPLETE.json`、`refinement8/{COMPLETE,SIMULATION}.json`、`hollow_cubic_shell_rods/COMPLETE.json`、`original_size/{COSTS,CROSSGRID,TRUNCATION}.json`。fine10首次生成的旧进程live sample曾见约1.2GB峰值；cached四臂的428MB不包含再次生成fine10。上述均不是χ8、新strict-A编排或8³×240PR的实测费用，不能线性包装成承诺时间。

缓存应明确开启 `cache_physical_values=True`；原坐标应明确 `whiten_apcs=False`。历史8³重复callback LU从653降至350、CPU227.39→101.28，证明这项既有paid-state复用有实际价值，不要重新引入重复物理求解。旧失败/中断费用在R8保留。

### 为什么原电尺寸不能仅改width或放大dense grid

原2019为3λ盒、60偶极/60矢量接收、forward60³/inverse30³、10%AWGN，coefficient相对变化的sqrt(mean ratio)停止1e−3；MF6→10，beta6→1，relaxed current beta依对象为6或3，hollow另有MF4/8/10的artifact控制。R8 .6λ盒是线性尺寸.2倍的电尺寸改编，不是原强散射条件。

`VectorCIE3D.__init__` 的接收/源半径.6、offset和对象坐标被硬写；只传`width=3`会改变网格但不改变acquisition/geometry，可能把源/接收放进DoI，不能这样声称原尺寸。必须在R9共同显式参数化grid、object、source/receiver位置、入射/S组装及E2，且检查源/接收在物体外。

下表仅为单个complex128 dense矩阵的算术内存，不是实测总RSS：

|grid|材料N|3N|单个G/L矩阵|
|---|---:|---:|---:|
|8³|512|1536|37.75MB|
|12³|1728|5184|429.98MB|
|16³|4096|12288|2.416GB|
|30³|27000|81000|104.98GB|
|60³|216000|648000|6.718TB|

实际dense实现同时留Goff/G/H/K/LU及临时数组，不能按单矩阵估总RSS。原MF10的显式projected-IDFT基在30³还需要约31.1GB；原尺寸必须用真正native 3D Fourier actions和matrix-free dyadic传播，不能把显式B放大后叫FFT加速。`VectorFFT3D`可复用padded传播/adjoint/GMRES骨架，但还需要严格CIE/APCS/材料读出接线，不能称现成成熟3D CIE。

12³/16³原尺寸已有36.1%模型差异，即使某个coarse材料NRMSE看起来小也没有关闭物理资格缺口。30³/60³规模/费用尚未测量，本文不承诺时间，也不新建CUDA/服务/调度框架。

## 7. 建议root的下一步裁决顺序（本文不执行）

1. 确定四个几何族、材料值与公共盒，以及是否把原2019有损coated作为独立保留控制；区分原论文参数和R9 χ8改编。
2. 明确目标是“小网格全矢量有限DDA开发”还是“原3λ电尺寸论文复现”。不能用前者的快费用证明后者完成。8³是现有最快有实测细化/费用依据的非平面全矢量入口；薄壳/杆和χ8资格缺口必须保留。
3. 在R9只实现新的strict current/material编排、公共box、raw-noise DPCS hook及明确的public-only数据接口；不修改R8，不复用会自动读truth的实验main。Common original Phi、prior、全未知空间与source metrics保持一致，保留raw/joint/physical-corrected输出和已付费轨迹。
4. root审阅新接线与物理资格合同后再冻结预算/运行。若未来需要一项新derivative wiring check，仅针对新增严格split/box/selector接线，不重复旧误差floor扫参。物理验证、算法成熟、材料恢复与原论文复现分开报告。

本审阅未生成任何新数值结果，也没有把R8的A-current混合结果改写为严格A-material证据。
