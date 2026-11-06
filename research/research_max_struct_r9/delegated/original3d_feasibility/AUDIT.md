# 原2019四立方体3λ：只读可行性与复用审计

**实际后续状态已更新于[EXECUTION_STATUS_ZH.md](EXECUTION_STATUS_ZH.md)：原尺寸60³数据完成，发现30³/60³ DDA响应差53.5%；父线程要求先诊断物理，五臂暂缓但未取消，GPU已释放。以下保留最初实施前审计。**

结论：**30³反演/60³生成在内存上有可信路径；现有R8接口不能直接运行成熟原尺度CIE。必须复用原dyadic FFT核，增加小RHS分块、流式外部算子和隐式原Fourier系数作用。不能把旧dense VectorPR3D简单改width=3。** 当前没有运行探针或反演，避免与其他四个数值进程竞争。父线程随后已授权R9原尺度扩展，以下为实施前证据及边界，不是已完成30³/60³的声明。

## 原尺度与几何

已读原文本地全文 `../r7_worktree/research/research_max_struct_r7/delegated/public_literature/CIE.txt`、R8 `PAPER_CONTRACTS.json` 和 `VECTOR3D_REPORT.md`。2019合同为3λ立方ROI、300MHz/λ=1m、60个电偶极源、同位置60个接收点×3分量=180接收通道，合计10800复观测。论文将数据写作60×180，本实现可用其转置180×60。三个半径3λ的圈分别为xy/yz/xz，中心(.2,0,−.1)/(.1,0,−.15)/(−.05,.1,0)，对应z/x/y偶极轴。四块边长.8λ、中心(−.7,.7,−.7)/(.7,−.7,−.7)/(−.7,−.7,.7)/(.7,.7,.7)，epsilon2/3/4/5，即contrast1/2/3/4。原生成为60³，反演30³，MF6→10、material beta6→1，relaxed第二阶段current beta3，10%噪声。

两个现有类的width语义不同，必须显式修复接口，而不是缩放物理问题：

- `R8/src/vector_fft3d.py::VectorFFT3D(ns,width=3)`：ns是**每轴网格数**，不是source数；P固定60。width只改变h/网格，半径3、圈心和四块几何是固定原尺度常量。width=3时匹配原几何。
- `R8/src/vector_cie3d.py::VectorCIE3D(ns,width=.6)`及其VectorPR3D子类：源圈半径始终硬编码.6、offsets为原文的.2倍；改变width不更新源几何。直接传width=3会把源放进3λ ROI，且仍构造巨大dense G/B；这是不可用的原尺度入口。

原文说明DPCS强子空间按噪声选L；不能把R8探索代码的固定相对阈值.1冒称原文唯一指定值。父线程的新原尺度合同选固定.1作为保守、明确披露的诊断规则。

## 可实际复用的API

| 现有入口 | 可复用部分 | 不能直接复用部分 |
|---|---|---|
| `R8/src/vector_fft3d.py::VectorFFT3D.G(v,adj=False)` | R1 dyadic近/远项、CM+RR自项、2n零填充FFT缓存；conj+分量转置的adjoint | 仅单RHS，固定reshape；需新增G_block(M,b)及伴随，不构造dense G |
| 同类 `.solve(x,rhs,adj=False)` | SciPy GMRES、对角预条件、原/伴随算子次序、逐RHS真实residual检查 | 逐源循环，restart60×8=480步上限；没有大型dense fallback，不能假定高contrast细网格仍会收敛 |
| 同类 `.physical(x,Y,gradient)` | 原完整数据目标与解析adjoint结构 | 原prior非原文合同，原callback会重复forward；新paper-reference必须R=0、阶段末验证、保留已付费状态 |
| `R8/src/vector_pr3d.py::VectorPR3D.current_pr` | sourcewise/global PR+、non-descent reset、固定材料精确二次line step、源能量权重 | 构造dense K、B、SB、KB；不是原尺度可执行实现 |
| A18_C `route-a/shared_backend/fft_maxwell.py::backend_from_arrays` | 已有多RHS卷积布局、conj-transpose传播，可参考分块接线 | 参数为原a/kernel/B合同，不能把chi当其a；CUDA L每次同步回CPU，SciPy GMRES不驻留GPU。其加载路径还含hash，不调用该加载器 |

已读 `research/research_z_alg_cuda_pipeline/OPTIMIZED_OPERATOR_CATALOG.md`。其576-cell标量resident LU、A17压缩endpoint不适合81000/648000矢量未知量。当前可信首选是现有CPU SciPy FFT/GMRES扩展；不以未运行GPU路线承诺加速，不安装后端或服务。

## 隐式Fourier/APCS与成熟sourcewise PR接线

令T_mf为三分量、orthonormal 3D-IFFT的原始低频掩码作用，V_r来自实际180×3N接收算子，P=I−V_r V_r*。存储原系数gamma，定义

```text
U gamma = P T_mf gamma
U* w = T_mf* P w = gather_mask FFT(Pw)
J = Jd + U gamma
Jd = Vr diag(1/sigma_r) Ur* Y
H gamma = p [U gamma − diag(m) G(U gamma)]
H* w = U* [conj(p)w − G* (diag(m)conj(p)w)]
p = beta_current/(1+beta_current m), repeated over3components
```

源归一data/state权重分别由公开Y/E逐列能量固定。每源独立PR系数与line alpha（Global只在配对臂改为全源内积），不得通过APCS正交化/SVD whitening改变gamma欧氏度量。MF6→10通过相同频率索引零填充gamma，Jd和V_r不变；阶段转换重置PR history。原文Eq16的MathML包含msqrt：delta=sqrt(mean_s||delta_gamma_s||²/||gamma_old_s||²)<.001，等价平方量门1e−6。先前把平方量门写作.001是继承误读，已在任何新反演前更正，来源见PRIMARY_STOP_CORRECTION.json。R=0时普通共享实CIE材料读出为bounded real modified-contrast LS，不需通用材料L-BFGS。

反演S可存233MB，receiver Gram仅180²；可分块累计S S*，在小Hermitian Gram上取特征对，再计算Vr=S*Ur/sigma，避免大B/SB/KB或全宽APCS SVD。强模态截断下应检查Vr正交及新U/G/S adjoint接线（父线程只授权一次有内容的新接线检查）；不宣称Gram平方条件数与直接SVD在所有弱模态等价。60³生成不需要DPCS或receiver SVD。

## 内存算术（complex128，十进制GB；不是实测峰值）

| 项 | 30³ inverse | 60³ generator |
|---|---:|---:|
|3N|81000|648000|
|dense G（必须避免）|104.98GB|6718.46GB|
|dense S，180×3N|.233GB|1.866GB（可stream）|
|单份全60源field/current|.0778GB|.622GB（可stream源）|
|dyadic kernel，(2n)³×9|.0311GB|.2488GB|
|一个4RHS padded FFT buffer|.0415GB|.3318GB|
|一个全60RHS padded FFT buffer|.622GB|4.977GB（多个并存会吃尽24GB）|
|GMRES restart60的61个单RHS向量|.0791GB|.6324GB|
|显式MF6 B|6.718GB|不需要|
|显式MF10 B|31.104GB|不需要|
|MF10原gamma×60源|.0230GB|不需要|

G_block采用b=2/4并复用kernel，不能一次pad全部60源。fine60生成可逐源/小组构造incident、求解、流式receiver S@J，只留下180×60的Y；无需同时保留完整S、inc和全部J。S/inc的dyadic临时tensor按source/receiver组块构造，避免现有 `points×60×3×3`约1.866GB暂存和隐式reshape副本。

局部A+的2³块只需同一个24×24平移不变off-diagonal Green，cell自身3×3块置零；不构造全域dense H。按8种voxel颜色/分块数组进行实scalar-isotropic材料候选，全部3分量/60源共同形成条件q/c，使用原固定current loss接受。它不是3D完整forward的rank-one结论；复介质不套此实圆弧。

## 已付费成本，不杜撰原尺度时长

原3λ R8实际16³生成用了18.0222s wall、60 primal RHS、5104 GMRES迭代、5292 FFT matvec，最大线性残差约2e−6。12³过程累计232717 matvec、3100 primal/1560 adjoint RHS，整个尝试379.5928CPU s/395.8676wall s、356040704B RSS；只完成20初始＋3accepted corrector，**不是完整求解**。16→12相同四块forward差36.113%，精度不合格。

仅作条件规模情景：假设完全相同GMRES迭代数与FFT效率，n³log(n)使60/16单matvec约73倍，18s生成约22分钟；30/12约20倍。它们不是实测时间、上界或收敛保证；更细网格高contrast的迭代数、cache和批处理效率均未知。不能由旧.6λ小模型秒数推断3λ成功。第一批60³真实source工作应直接写入生成progress/receipt，由已付费源估算余量，若不可承受，先交实际时间/RSS/残差，再由父线程决定中间独立网格。

## 自项与论文忠实度

现有FFT核与R1 TensorDDA一致：offdiag为h³ dyadic_green，self为 `s I3`，`s=−1/3+i k³h³/(6pi)`。隔离单cell对应CM+RR polarizability `alpha=h³chi/[1+chi/3−i k³h³chi/(6pi)]`。这是point-dipole/CM+RR离散化，**没有证明等于2019作者的体积分MoM自项/动态cell积分**；不得标作者软件复现或连续精度认证。

30³ h=.1λ、epsilon5内部波长约.447λ，仅约4.47cell/内部波长；60³ h=.05λ约8.94cell。原四块的.8边长及±.7中心在30/60网格上边界对齐，可避免任意partial-cell混合；原R8 volume-fraction chi近似在不对齐中间网格仍是额外误差来源。不能用同离散成功代替独立60→30或独立参考验证。

## 推荐的一条实际路径

1. 新增R9扩展，复用VectorFFT3D kernel/原source圈及R1 dyadic，提供G_block/G_block*、chunked incident/receiver作用；先完成代码，不与当前四进程抢CPU。
2. 在一个数值进程可用时执行原3λ/60³四块生成，源组2/4、逐RHS成熟GMRES并真实记录所有残差/迭代/失败/RSS。保持尺寸、源60、通道180，10%IID圆复噪声及公开sigma为新明确实现，不冒称作者噪声实现。
3. 30³全部27000实chi自由、air初值、[0,4]、paper refs R=0；隐式MF6→10、阶段各最多500步和Eq16平方根.001门（平方量1e−6）。先运行sourcewise CIE及relaxedDual，再配GlobalDual与同ordinary current的A+global；阶段末完整物理验证，初期不额外加昂贵corrector。
4. 生成/反演首次实际费用若超出可承受范围，给父线程真实receipt后再定中间gen/inverse网格，仍保持原3λ及原物体/采集。优先争取30/60；不默认退回电尺寸缩小。

本审计未执行任何计时探针、数值生成或反演；内存为维度算术，时间只有上述已有R8实测与显式条件情景。新实现及科学合同由父线程最终核准。

## 实施前主文MathML更正

已直接重新取得PMC HTML并解析MathML：Eq15没有额外rho（data/state1:1），Eq16与Eq17均各有一个msqrt。新paper-reference及matchedA固定rho=1、R=0；不能继承R8小型vector的rho10。Eq17报告sqrt(mean_voxel |epsilon_est−epsilon_true|²/|epsilon_true|²)，不是平方误差，也不是chi-relative-L2。错误旧门、原MathML和更正前反演数0保存在PRIMARY_STOP_CORRECTION.json。原文链接：[2019 primary](https://pmc.ncbi.nlm.nih.gov/articles/PMC8320915/)。
