# R4 物理数组与有限反馈审计

状态：`ALGEBRA_IMPLEMENTED_PHYSICAL_PILOT_NOT_RUN`。新增 Maxwell state / physical complex RHS / physical LU 均为 0。没有导入 A30 后端、构造 PhysicsBook、执行历史 stage 或修改旧 ledger。

## 冻结合同和数学

P1 target 为 `P1-CORRELATED-HOST-MIXTURE-TM24-v1-20261018`，prior 为 `P1-Gaussian64-Uniform2-C2-volume-fraction-v1`。24×24、576 材料单元、半径 .38 ROI、8 source、32 receiver、k=2π、gain=1。宿主 .1+.04i；增量实部 [.8,1.8]、虚部 [.03,.18]；Gaussian64×Uniform2，C² quintic，halfwidth .5。宿主 pocket 不是空气。

数据为 `S @ current`；packing 为逐 source `Re32, Im32`。σ=.00119211798514519，复噪声实/虚各 σ²/2；白化 `sqrt(2)/sigma * pack`，因此 Phi=.5||r_white||²=||F-y||²/σ²。probit x=(z,η) 导数必须再乘 du/dη；`compose_jacobian` 按调用者提供的 D 所属坐标解释，不能将 native Jq 冒充 Jx。

实际方程 `Lχ=I-diag(χ)G`，`Lχ j=diag(χ)e_inc`，`E=e_inc+Gj`。有限反馈注入必须是 `Δχ E0`。实现 `c=(I-diag(d)K)^-1 diag(d) B`，`Δj=Wc`，`ΔF=Oc`，`ΔE=GWc`。这里 W=L0^-1 E_I 由唯一物理 owner 支付并供应，模块不求它。cache 初始化同时检查 base/current/field/receiver 和 W 的方程残差；失败拒绝复用。

`active` 是相对原 cache base 的累计 changed-set 的超集。严格用实际存储非零判断，不偷偷截断小值。候选跨事件后已饱和的格子仍可能与 base 不同；属性变化也可能影响全部非零 fraction。少于66维 latent 或4个弱方向均不能给材料增量 rank。A17X/F3 的历史反证：17个 Δχ 为575个存储非零；3个 weak 端点在相对1e−10阈值仍为571。旧575阶增量不是P1事件稀疏性的正证据。

Jacobian 演化实现三项精确望远镜分解：(T1−T0)E0 D0、T1(E1−E0)D0、T1 E1(D1−D0)。其中 T=S L^-1。接收 T0 不能直接提供 active inverse columns W；W 本身通常也不足以生成完整 T1，还需要额外的 inverse actions 或付32 receiver adjoints。不得声称这里只更新了 resolvent 就得到了完整 J transport。

## 可复用文件和成本

完整绝对路径、文件数、样例 schema 见同目录 `CACHE_INVENTORY.json`。历史 R2 目录：`../b1-p3l-r2-risk-posterior/research/p3l_r2_risk_posterior_exploration`。盘点为路径条目，不是独立状态/parent数量：9,403 cache 文件；8,251 receiver_J；1,408 RAW_MODE_CACHE；1,408 JOINT_MODE_GEOMETRY；32 A/anchors。同一状态可能在多处重复。

可直接读取的 EXPOSED DEV 例子：

- `A/anchors/DEV_010.npz`，也有 DEV_004、DEV_005；含 q/chi/current/field/prediction/LU/pivots/backward_residual/receiver_T。
- `C/jobs/C_exposed_eval_r0/trajectories/DEV_025/dynamic_curvature/receiver_J/*.npz`；含 q、512×66 J、576×66 D、32×576 receiver_T。与同轨迹 cache 按文件 stem 对齐。
- `inputs/background/weak.npz` 中 G、S、incident、几何、descriptor 为冻结数组；R4 新 checkout 没有此未跟踪文件，需只读引用旧 R2 文件。

清单 schema 也抽读了一个 L/locked 缓存；这只验证历史存储形状，仍为 EXPOSED，不是新 DEV 或新 locked。

一次 Full forward=8 complex RHS+1个576阶 LU；完整 J 在已有 valid state 时32 receiver-adjoint RHS；gradient 为8 source-adjoint，已有 J 后仅矩阵乘；JVP 通常8 RHS，已有 T 则纯乘法；每个 active inverse column=1 RHS。四个 JVP 已抵 Full J 的 RHS，wall 则须实测。

A30 `cached_anchor()` 自带8 RHS 的 LU审计，不能误称0 RHS加载；R2 `Target.load_cache()` 的方程审计是矩阵乘，手工重建 native state，0新 RHS。旧 P1/P2/R2 Book 有闭场与stage guard，不能用于本轮，也不能改名绕过；由父线程建立独立新 campaign Book。

`FiniteMaterialCache` 显式保存历史 inverse-column acquisition RHS/owner，所有 cache 数组及 K/O/B/GW 的实际 bytes，setup arithmetic wall 和 evaluate wall；每次 evaluate 记 reduced factorization / small solve RHS columns，不冒充物理 RHS。`cond(small)` 是监控而非稳定性证书。新 rebase 必须另付与新 base 匹配的 current/field/inverse actions，模块从不自动重基或调用 Full。break-even 粗略需显著多于 |I|/8 次候选复用，还必须计小系统因子分解、条件监控、矩阵乘、内存、失败和重基。

## 实现范围和资格

`src/p3r4/reuse.py`：FiniteMaterialCache；J composition/evolution；good-Broyden；有限 secant 的 decoder/Maxwell 向量余项与 cross inner product；原 prior gradient 的 Gaussianized GN action；外供 Full 数组的观测误差对照；直线和 ellipse 事件根。模块不生成候选物理响应，也不独立求大系统。

复用既有 `p3l.decision_bounds` 的 saved-inverse Neumann 下界、primal/dual output bound、同一个 uniform 的 MH 区间判断。它们明确为 `CONDITIONAL_FLOATING_POINT_BOUND_NOT_INTERVAL_CERTIFIED`，有限 normal float64、标准gamma_k算术模型等条件没有经 directed rounding 验证。普通 condition number、随机 probe、训练分位数不构成确定性 β。没有可靠界时保留 `BOUND_UNAVAILABLE`/`REFINE`，由父线程 Full；本轮未产生确定性certificate。后验统一误差界还要求全支持或受控tail域，不能从几个端点推导。

`src/p3r4/geometry.py`：实坐标 `[X,Y]=DY X−DX Y`，声明 D metric 的材料绝对幅度/normal fraction；soft spectral mobility；带正 floor 的全维 covariance；含实际 determinant 的 forward/reverse Gaussian log density 与完整 MH ratio。禁用伪 determinant。场/metric必须是生产前冻结的可计算状态函数；历史cache依赖需要另行适应性或增广状态论证。模块未把四段 Euler loop 称为可逆 MH map，未声称 Lie full rank、连续路径可行、Hörmander、mixing 或 posterior 资格。

P1 seam只有 C²，固定光滑stratum内可作更高阶展开；hard projector在谱隙条件下通常只有C¹。弱frame符号与簇内gauge应通过projector处理，normal fraction不能替代绝对材料幅度和数据/prior成本。64场与66含属性几何须分开。

## 已测试与未运行

`tests/p3r4/test_physics_algebra.py` 共9项 PASS：9阶通用复矩阵有限恒等式/三项J演化；累计支撑和属性变化 fail closed；空支撑/输入所有权；stale inverse拒绝；Broyden secant；向量余项；解析事件根；有条件残差区间和同uniform判断；全秩Gaussian正反proposal以及Heisenberg bracket约定。测试只用通用合成矩阵，不是 Maxwell 物理状态。Full/reference矩阵解仅在此小合成fixture里出现。

物理pilot、收益/摊销、material MCSE、决策错误率、新parent统计、学习模型、Krylov和DINO对照均 `NOT_RUN`，不得据9项代数PASS升级科学结论。
