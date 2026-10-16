# 给 Codex 的下一轮 Prompt
# Research F0 — Intrinsic Prior Geometry Audit
# 默认只读与缓存代数；不运行新的 Maxwell states，不训练网络

你在 migodam 的 b1-probabilistic-maxwell-imaging 项目中工作。
先读随附 Research_F/report/RESEARCH_REPORT_ZH.md，尤其第 6、7、9、10、16 节。
本次目标不是复现已经付费的 finite_feedback_fiber_diagnostic，也不是扩大多解实验；
目标是判断旧报告的巨大 decoder 曲率中，有多少只是参数坐标加速度，
有多少是真实材料 prior 的第二基本形式，以及数据余核中不可补偿的二阶响应。

## 执行权限与硬边界

RUN_MODE = "audit_only"
MAX_NEW_FULL_STATES = 0
MAX_NEW_COMPLEX_RHS = 0
MAX_NEW_MIE_CALLS = 0
MAX_GPU_TRAINING = 0
MAX_NEW_3D_STATES = 0

允许：读取已公开/已暴露证据；decoder-only 求导/有限差分；对已保存 J 与材料数组做矩阵运算；
已保存状态方程的残差代数检查；随报告的六维/两单元代数单元测试。
不允许：重新求解旧 Maxwell states、新建 full Jacobian、调用 R0/LU solve 生成缺失响应、
重新运行 LXZ/Joint P2、重训 N3/C1、消费旧已关闭预算或修改 RUN_CLOSED。
所有新增 wall time、缓存读取、decoder evaluations、matrix products 单独记录。
缺数组就写 MISSING，不能偷偷重算并称为零物理。

原始证据入口：
1. b1 commit c76733676a7746627dfbf710c5f31b6375ae4335
   research/finite_feedback_fiber_diagnostic_20261016/
   START_HERE.md, DIAGNOSTIC_REPORT_ZH.md, METRICS.json,
   LOCAL_PATH_AND_CURVATURE.md, REPRODUCE.md, POSTHOC_CONSISTENCY.json,
   THEORY_QUESTIONS_FOR_PRO.md，以及其声明的公开 release/raw evidence。
2. LXZ commit dde0237eaa0d8c5507fedb07bb704c66a817028d
   research/research_lxz_20261016/RESEARCH_REPORT_ZH.md
3. A26 commit 2603f3bd16c48e48587fedc4363975779def17f5
   src/a26/decoder.py, src/a26/training.py, configs/a26/pilot.json
4. A17 closure_final/manuscript/CORE_THEORY.tex（只作理论接口参照）

新建独立研究目录 research/research_f_intrinsic_audit_20261016/，
保留原证据和实际执行版本 d74f747 的身份，不覆盖历史报告或修正其原账本。
先按 REPRODUCE 找实际可访问的缓存，输出 SOURCE_INVENTORY。
不要根据一个文件名推测内容。无法取得原 NPZ 时完成源码/数学审计后停止。

## A. 冻结正确几何与数据定义

沿用原 24x24 scalar-TM、8 sources/32 receivers、已知 gain=1、
源内 Re32/Im32 packing、sqrt(2)/sigma 白化、
材料 sqrt(cell area)*[Re chi, Im chi] 质量坐标。
禁止改变 noise、material box、P0 相位、原 rasterizer 或 max-normalization。

两个已知 TRAIN anchors 250033_weak、250034_ring 仅作几何校准。
它们自己的 GT 就是起点，不得把离开起点写成 GT 恢复。
N3、physics_R_base 是另列 oracle/chord 诊断，不能当合法 P0 anchor。

复核原 8 参数 Gaussian-CDF map 和同材料 scale gauge。
对 Dg 做 SVD 得到有效材料 rank 与正交 tangent frame；
不能把完整像素 nullity 搬成原 prior 自由度。
Geometry 可使用 quotient frame，但不要把剩余 7 坐标随便定义成均匀 prior；
本次不改变概率目标。

## B. 把旧 decoder 曲率拆成三项

对原报告已经使用且有完整缓存的方向 v：
Qraw = Dg
w = Qraw @ v
P_prior = orth(Qraw) @ orth(Qraw).T
d2 = D2g[v,v]
d2_parallel = P_prior @ d2
d2_perp = d2 - d2_parallel

从已保存数据取得：
C_Maxwell = D2F[w,w]     # 缓存缺失就 MISSING
C_coordinate = J @ d2_parallel
C_prior = J @ d2_perp

报告：
1. 三项向量范数、向量和误差、两两 cosine、方向性 cancellation。
2. 原报告的 C_Maxwell + J@d2 是否复原。
3. 用材料速度 ||w|| 统一归一化；二阶量按 ||w||^2 归一化。
4. material metric、prior metric 两套数值并列，不混单位。
5. A = J @ orth(Qraw)，Q_out = I - A A^dagger；
   不必构造512x512投影矩阵，使用投影动作。
6. ||Q_out(C_Maxwell+C_prior)|| 及 ||Q_out C_coordinate||。
   后者理论上应是数值零；若不为零先查 metric/packing/rank tolerance。
7. 明确当前 A 若满列秩，则没有非零 exact kernel。
   对旧 approximate-weak v 的二阶响应不能称“已证实 fiber 障碍定理实例”；
   它是固定初始速度下的 finite-noise 几何诊断。

D2g 若用有限差分，预先固定三个递减步长，报告收敛与减法消去，
active-max cell 改变、CDF 饱和和 rank 变化。
不得只挑最漂亮的一步差分。
保留旧 sign-tracking bug 及其已经发布的重放说明；
新后处理不冒称重跑了原轨迹。

## C. 加两个 decoder-only 对照

1. 平坦 image 的弯曲坐标单元测试：
   g(u,v)=(u+c*v^2,v), F(x1,x2)=x1。
   raw curvature=2c，intrinsic prior curvature=0。
   geodesic coordinate acceleration 去掉假曲率。
2. 原 P0 同 image 的局部坐标重参数化：
   只在光滑、rank不变的邻域内做可逆局部坐标变换，
   验证 intrinsic material curvature 与余核响应在相同材料速度下不变，
   raw JD2g 可以变化。
   不生成新的 full physics endpoints。

## D. 对 N3 条件均值解释做代码级裁决

直接复核实际 delta_loss 和 ssim_delta：
- shared material RMS；
- MSE 是复增量 abs^2；
- SSIM 对实/虚通道分别用固定范围 4.5/.5；
- full-grid box，不是固定相位投影；
- graph 无 measurement loss/full solver。

写清：固定相位条件下的理想 conditional mean 应保留相位；
N3 phase violation 不可只归因于 manifold averaging。
不重训，不重新计算旧完整成像指标。
不要把必要 phase 检查通过当作完整 P0 membership。

## E. 对下一阶段给可证伪裁决

若主要是 coordinate term：
建议先 material-native retraction/geodesic chart，再谈高阶网络。
若 intrinsic/cokernel term 显著：
给固定材料尺度的有限步风险预测，并说明需要怎样的独立 anchor。
若 gauge/不光滑/数值误差主导：
先修坐标，不建立“多解不存在”结论。
若所有已有 prior tangents 数值 full-column-rank：
保留 finite-noise posterior，暂缓 exact-fiber NN。

只设计，不执行新的 E2：
fresh cap 建议64 states/2048 complex RHS/0 GPU；
背景独立反演、合法 smooth-P0 条件目标、同材料距离对照、
GT 严格只作离线评分。没有新授权不得启动这个阶段。

## 输出

必须产生：
- SOURCE_INVENTORY.json
- CURVATURE_DECOMPOSITION.csv
- CURVATURE_VECTORS.npz（只保存衍生数组，注明原来源）
- INTRINSIC_AUDIT_ZH.md（结论、公式、原始数字、失败与缺失）
- REPARAMETERIZATION_TESTS.json
- N3_OBJECTIVE_AUDIT.md
- COST_LEDGER.json（new states/RHS/Mie/GPU必须显式为0）
- NEXT_STAGE_SPEC.json（默认 NOT_AUTHORIZED_NOT_RUN）
- tests/ 与可运行的只读审计入口

不要只交计划。能从缓存完成的审计实际完成；
缺关键缓存时交付已经完成的源码/代数检查与明确缺失清单，
不要自动转为新物理仿真。
不提交或推送 GitHub，除非我另行明确要求。
