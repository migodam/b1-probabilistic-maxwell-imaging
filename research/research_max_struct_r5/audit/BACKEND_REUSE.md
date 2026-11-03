# R5 后端、九数字与继承公式库存

2026-11-03；本文件只读历史源码、合同和回执后新增。新 Maxwell LU/RHS=0；没有运行 R5 三维成像、没有打开 truth 材料/support 数组，没有在线选择。稀疏缺失文件用 `git show HEAD:<path>` 读取，没有改变 sparse checkout。机器可读路径、访问方式和九数字映射见 `BACKEND_EVIDENCE.json`。本库存不作最终原创性或科学方向判断。

## 1. 可以直接复用的真实三维矢量后端

实际实现：`research/research_max_struct_r1/src/vector_backend.py`，Green 核来自同目录 `a9_green.py`，明确保留 A9 来源。

```python
model = TensorDDA(points, volume, k=2*np.pi)
state = model.solve(epsilon, directions, polarizations, receipt_path=None)
S_E = model.external_electric_operator(receivers)
S_H = model.external_magnetic_operator(receivers)
near = state.near(receivers)
far = state.far(unit_directions)
```

`points:(N,3)`、等体积 `volume>0`，`epsilon:(N,3,3)`，每源方向为单位三维向量且极化横向；内部列场为 `(3N,P)`，cell-major/component-minor。`state` 保存 E/H、j=(ε−I)E、J=−ikj、p=volume*j、exciting E、matrix、LU factor 与回执。`near` 返回总/散射 E/H，接收器必须在占据体素外；`far` 返回满足 E_sc(rn)=exp(ikr)F_E(n)/r+O(r⁻²) 的远场幅度，不能与有限距离近场混用。

物理是 exp(−iωt)、λ0=Z0=1、μr=1 的全矢量 CM+RR point-dipole DDA。`G=volume*Goff+(-1/3+ik³volume/(6π))I`，解 `(I−Gχ)E=Einc`，全部照明列共享一次 dense complex128 CPU LU。H 是电偶极 Green 的解析 curl，内部对称 cell principal self 为零。不是标量 Helmholtz 冒充三维 Maxwell，也不是有限 Born 截断；但它同样不是兼容 de Rham/H(curl) 离散或连续误差证书。

材料类别：

| 类别 | 可调用方式/实际证据 | 边界 |
|---|---|---|
| 实各向同性 scalar ε | 将 ε_i I 传为 tensor；实际 Mie 控制使用 ε=3 | scalar 是矢量 Maxwell 的材料限制，不是二维 scalar-TM |
| 实轴对角 SPD tensor | R1 `tensor_families.py`：I+b zzᵀ 与 I+3b zzᵀ | 已披露解析/开发同场控制，不是反演 |
| 复对称严格被动 tensor | R1 T=实SPD混合张量+.25iI，qqᵀ方向最大cell算子范数 .12，损耗裕量至少 .13 | 有限DDA数值分支；尚无该张量的独立连续solver对照 |
| 任意 reciprocal complex symmetric ε | API检查有限值、shape、transpose symmetry | API没有自动强制被动/SPD/box；调用方须显式约束，不能把可输入当科学资格 |

## 2. 独立验证已经到哪里

`experiments/qualify_backend.py` 的 `mie` 是解析 Riccati–Bessel 球参考，采用 Hankel h_n^(1)，不是相同DDA的另一组装。实际球 r=.12λ、ε=3、μ=1；`results/backend_qualification/QUALIFICATION.json`：

| h/λ | N | far vector relative error | scattering cross-section relative error | linear residual |
|---|---:|---:|---:|---:|
| .04 | 123 | .1046792878 | .2180790022 | 1.40e−15 |
| .03 | 257 | .0365537162 | .0598918396 | 2.01e−15 |

两误差随该次细化下降；远场横向误差约1e−16、1000λ近远差约3.18e−4。历史 qualification 为 `BOUNDED_DDA_QUALIFICATION_ONLY_NO_CONTINUUM_CERTIFICATE`，`M_Hcurl_status=NOT_RUN`；没有任意人为合格阈值。这是一个无损各向同性球的有界控制，不能直接升级为任意高对比有损/各向异性成像或收敛定理。

R1 同场分支的每网格核重新生成，不代表固定连续 tensor 在网格加密下同场收敛。固定粗方向插值控制在 R1 另外记账。R2 `experiments/active_family.py` 的三维实验是盒边 .24λ、64cell（细化512）、已知 b 支撑半径 .1λ、ε=I+t b zzᵀ、t∈[0,3]，仅一未知参数。旧源沿z、x/y极化；新增源由不同方向/极化通道相干组合，单位系数范数不是实际天线功率标定。12外部xyz近场接收、确定性复球噪声合同与二维Gaussian噪声合同不同。该演示说明限定族中的新观测价值，不能当自由体素三维成像。

## 3. FBG-R3九数字真实输入和约束

九数字是历史 FBG-R2 Stage4 已披露数据，根为 `research/research_fbg_r2_adaptive_feedback/data_contract/stage4_public/`；逐个匿名观测仅含 `Y:(48,12)` complex128 与 observation_id。算子文件 `operators/uniform12.npz` / `random12.npz` 包含 G:(576,576)、E0:(576,12)、S:(48,576)、二维centers、外部sources/receivers、k/radius/grid_shape/material_bounds等。已用标准库读取公开NPZ头部核实这些shape；没有新增数值运算。无须从GT生成新的数字或读truth进行在线选题。

完整九数字 label→匿名ID在 JSON，仅供离线展示，不能作为solver输入。该映射读取 `data/e3_eval_metadata_portable.json`，其中已有披露background等元信息；未用其值选择参数，没有打开其 `truth_file` 指向的数组。

FBG-R3 CPU API `src/physics.py:MaxwellModel` 包装冻结 R4 `Model/build_operators/Regularizer`，有 `forward/state/jacobian/tangent/regularizer/physical_score`。这是 exp(−iωt)、k=2π 的二维 scalar-TM disk-pulse VIE：j=m(E0+Gj)，F=Sj。全部576实材料自由，完整resolvent保留多次散射。它不能直接换作三维vector tensor后端。

E3真实runner为 `experiments/e3_fullwave_imaging.py`，不带 execute 返回NOT_RUN；历史Windows/CUDA依赖和当前CPU API要区分。精确已执行源还在 `results/e3_limited_diagnostic_20261101/executed_source_snapshot/`。历史两份manifest合起来为12 clean/amplitude12对象×4方法=48端点；pilot digit2四任务不重复。九数字是其中九对象，不是Z-ALG-R3另一套匿名九数字。CPU重构运行方式若不同于此历史执行必须重新标明，不冒称完全同一路径。

共同约束：公共常数初始化m=.1、box[.02,20]；无数字模板、无support先验；smooth与平滑TV各1e−4、差分尺度19.5、TV epsilon=.01。在线最终选择为原数据 `0.5||F−Y||²/||Y||² + R(m)`，不能以辅助current目标、纯residual或GT替代。若改成单位白化目标，R也须作对应尺度变换。每任务45s/2000LU/60000RHS是共同上限，不代表实际费用相等。Native/feedback保留8→12→24材料日程，CIE全576材料与current空间、mf3→6→10及最终共同full-material data corrector；不同自由度/日程不能被称同容量对照。

历史成熟dual-CIE为10/12 joint，九数字全部通过、ring/curved_bar失败；这是成熟基线诊断，不能归给N1，N1_E3/N2_E3保持NOT_RUN。GT成功率在历史报告已经披露，因此这些对象不能再叫fresh blind。

## 4. A17-X/P3L/Z-ALG 已有的数学，不可重新包装

下列文件在HEAD存在但稀疏工作区缺省时，用git show读取：

- `research/p3l_r10_bilateral_structure/BILATERAL_OPERATOR_THEORY.md`：T(D)=(I−DC)⁻¹D，F−Fh=OTE；精确差 ΔF=O(I−D2C)⁻¹(D2−D1)(I−CD1)⁻¹E。不需要D可逆。当前材料导数 dF=O_D dD E_D，Gram=(O_D*O_D)∘conj(E_D E_D*)，已知Khatri–Rao/Schur/resolvent公式，不是新的固定材料弱空间。
- 同R10 B13：近似current真实状态缺陷 r=DE−(I−DC)hhat，有精确 ΔF=O(I−DC)⁻¹r；任务伴随将误差写成 dual×r。若白化forward误差≤ε_f，则目标差≤||rhat_w||ε_f+ε_f²/2。没有数值可用稳定性/dual控制时只是恒等式/条件界。复被动固定形状的coercivity估计另需b_delta=Imδ/|δ|²+λmin Im_H B>0；不能对实scalar随意套用。
- `research/p3l_r12_bilateral_fixedpoint/BILATERAL_MODE_DERIVATION.md`：完整双边pair模态与Full-J局部信息等价。实先验链式法则 Hdata=(2/σ²)Re(Dx*GmDx)；局部decoder零方向须剔除，材料RMS分母仍576；加入原Gaussian prior曲率不新增观测信息。32 pair截断是启发式，不保证prior可达空间最优。
- `research/a17x_f3_competition_audit/THEORY_AND_SOURCE_COMPATIBILITY.md`：精确支持ΔX=UDU*下T=D(B+HT)、δj=PT、δF=OT成立；DDA中ΔX必须是实际离散本构/极化率的变化，不能混换。Schur端口更新与同一增广空间普通Galerkin完全等价，full−new尾项未消失。状态缺陷ρ_D=R_U Q_V B，theta=||Q_V||eps_H<1只是尾界适用条件；不满足不是物理奇异或材料不合法。
- 同A17-X：固定多源/三分量必须共享一个scalar材料；逐源或逐component单独拟合不保证材料可实现。已检验某些缓存端点直线current插值不满足共同scalar/固定相位，不能把端点同响应升级为连续gauge。
- `research/research_z_alg_r3_digits/REPRODUCE.md`：该独立历史九数字有90s合同、public_data/holdout输入，不能混成FBG九数字。用总场系统与电流系统验收是同一离散VIE不同组装，不是跨离散/连续验证。
- `research/research_max_struct_r1/THEORY_DERIVATIONS.md` §F：VU0=0给有限同场族；精确finite response与resolvent条件射线界已存在。局部tensor本构核谱不含接收器，不是完整Maxwell Jacobian谱；弱方向也不自动是可积分的材料纤维。

这些已有公式可以组成R5的正确实现和误差审查，但其成立本身不提供新颖性、自由三维成像收益、GT恢复或posterior校准证据。原创性与科研裁决留主线程。

## 5. 本轮仍未做的事情

本库存未开启新物理、GPU队列、科学参数扫描、GT评分或发布；未读取truth数组。尚缺任意三维成像任务的冻结prior/box/初始化、同后端独立基线与预算合同；没有将R2一维family或Mie球验收替代这些工作。系统python缺NumPy的只读尝试失败已保留为环境事实，没有安装包或启动求解。
