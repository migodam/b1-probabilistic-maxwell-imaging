# LXZ E0：二维既有物理资产与合法 anchor 审计

**原 development 观测中，现存 paid caches 尚未给出一个可直接用作原 P0 机制实验起点的、P0_EXACT 且 noise-compatible 的在线重建。** 本分工没有启动二维新物理。提供两个已知 TRAIN 材料的公开校准 anchor，必须与未知对象重建严格分开。

## 1. 数据、先验与接受条件

24×24 scalar TM/VIE，8 源、32 接收，complex128，gain=1。原复噪声 RMS 为 sigma，每实通道白化 `sqrt(2)/sigma`，512 实通道；`ell=.5||W(F-y)||²`，tau=23.785449004057032。各观测固定自己的既有 sigma，没有候选重新归一化、噪声放宽或 nuisance 变更。

原材料为 `chi=chi_b+v*rho`，`v=(1+.2i)/sqrt(1.04)`，a1=1.2；rho 非负且由原四类 TRAIN 生成器逐图峰值归一化。rho=0 是回到非零背景，不是空气。PHASE 必要条件、BOX、noise ball、点似然、P0_EXACT 和 posterior mass 是不同证据层级。只有实际原参数支持与端点精确重放才作为 P0_EXACT。

来源首先核对 `../b1-prior-constrained/research/prior_constrained_20261015/PRIOR_SUPPORT_AUDIT.md`、`SOURCE_MANIFEST.json`、`E0_E1_E2_DECISION.md` 与 delegated 原数组审计。前阶段已重放 TRAIN 512/512 条材料，误差0。原32 development 条件是 exposed；本次没有 locked 读取，也未改写历史结果。

## 2. N3：noise 相容与 prior 合格不能混用

从32个 A28 paid exact-state 文件直接读取 N3 prediction，对对应 pool 的 y/sigma 重算。**只有2/32进入95% noise ball，原 P0 的 PHASE 条件通过0/32。** 即使这两个数据相容端点也不是原 P0 合法 anchor：

|case|noise 相容|P0 必要 PHASE|离线 N3 与真 clean 的白化 prediction distance|
|---|---|---|---:|
|p250069_weak_a1_snr20|是|否|12.35540718523417|
|p250073_weak_a1_snr20|是|否|8.778176405627963|

全32条件的该 prediction distance 范围8.778176405627963–64.16729715102171。最后一列只在独立 `OFFLINE_N3_VS_CLEAN.csv` 用 development `clean` 计算，没有加载 development 材料真值、family 或 geometry 来指导在线候选。不能从 N3 落在 noise ball、或 N3 与真值共享同一个 y 推断二者具有相同似然/后验质量。

历史真 clean 对 y 的 noise norm 有4/32超过tau：250065 weak、250067 ring、250072 weak、250079 weak。保留其正常的噪声事件含义，不能称物理非法或改阈值消除它。

## 3. 可复用 P0 校准 anchor：两个，但不是 development 重建

主线程指定 genuine smooth chart 后，选择在算残差前固定为公开 TRAIN 的前两个 smooth parent250033/weak 与250034/ring、a1、snr20。用原 seed/family-local 规则重放原端点，两者逐元误差0，全部参数支持合法；未挑选最有利 noise seed。此前导出的250001/250017（shell/asymmetric）保留为历史准备文件，已从当前正式assets列表移除，不作为本轮E2起点。

|校准 anchor|TRAIN row index|既有 whitened noise residual|状态|
|---|---:|---:|---|
|TRAIN_250033_weak|263|22.313033017586633|P0_EXACT、noise-compatible、known smooth TRAIN calibration only|
|TRAIN_250034_ring|267|22.96957424921172|P0_EXACT、noise-compatible、known smooth TRAIN calibration only|

本目录保存它们的材料、paid clean prediction、既有 noisy observation、sigma、current 与公开 background；field只用缓存G乘已有current再加incident重构。两个同名JSON保存完整generator row及原先验证书，能够直接给原decoder，field/current/data的缓存方程诊断也保存于此。16个G source-column乘法及16个S source-column乘法是矩阵算术，不是full-solver RHS，不构建LU或Jacobian。它们来自监督 TRAIN 先验，没有 online discovery 含义，也不与原32 development y相容。不能用它们冒充已经解决 qualified base。两者没有已保存的自身 LU 或 exact anchor Jacobian，进一步付费准备必须由主线程单独决定。

E2若使用此两点，条件目标是原P0在family=smooth、amplitude=1.2及已知background下的分布；原四family中smooth权重.25被条件掉。参数保持原TRAIN独立uniform支持与实际raster-max归一化，inactive参数应积分掉，不变成field-uniform。smooth是前两个Gaussian和：shell/rectangle/第三Gaussian参数不活跃；两个Gaussian共同weight scale在max归一化后抵消，global/local坐标也可能冗余。这些坐标零方向不是物理材料弱方向。

另外，直接把每背景64个公开 TRAIN/a1 材料的 paid clean response 与32个 development y比较，作为不读 development truth 的有限 prior-reservoir 搜索。各 y 的最优白化残差仍为33.91839993065944–137.9364188410169，全部超tau。因此此已付费 reservoir 也未提供相容 development anchor；有限 bank 未找到不等于 P0 内不存在。

前阶段 `offline_only/E1_execution_2/f*_anchor.npz` 是已知生成材料的隐藏 benchmark 构造。可作清楚标识的离线物理校准，不能混成原 development 的在线反演或将其隐藏参数传给在线端。

## 4. 原 T/G corrector 的问题

已核对 T qualification/raw 与 G inventory 的原数组审计。T40个raw、G72个raw及所有保存 corrected/validation 记录，均不满足原增量相位/幅度/峰值必要条件。T 的绝对chi=0 core改变了非零背景且脱离原正增量生成分布；普通全材料修正不会自动把它恢复成P0_EXACT。

G采用32实参数 Haar GN/SLSQP、两步、BOX+trust，而不是原P0 retraction。72 validation记录只有36个精确不同材料，9个精确不同 noise-valid 材料；18条 noise-valid轨迹相对共同G base 的 Delta ell 最小13.370668、中位38.621476、最大64.453689，均不属于 log10 或 log100 的近点似然。G有效轨迹的 raw-parent 编辑保留投影中位数约0.006687，全部576像素被修正、低频/DC补偿显著；这提示corrector抹掉raw novelty，但该投影不能代替本轮持久结构验收。T虽保护严格空气core，也不能因此证明原P0支持或竞争 posterior region。

A29 pool 的预测字段为 L（背景线性响应）；历史2,100个候选缺 full prediction，因此 full ell 留空。不能从它们的线性残差挑选已验证的P0 anchor，亦不得为E0补跑full solves。

## 5. Weak/null 方向与 hard-raster 假零

原生成器有硬栅格化与按实际网格最大值归一化。几何参数改变而完整chi不变，只是同一离散材料的多重参数表示。此时任何response不变都不是非零材料weak/null方向，也不应以有限差分几何导数0宣称可辨识。

仅比较前阶段45条已付费P0端点记录：有4个同材料缓存配对（含重复运行保存及复用），没有发现不同材料而 prediction 精确相等的配对。排除同材料后，按每背景一个固定历史median sigma，最近既有 paid nonlinear secant 为d_pred=7.699171078917734、material RMS=.12247459814445766。这是既有offlinedata的有限pair诊断，未另验其结构门、不是局部Jacobian方向或竞争见证。它与前阶段“相对各固定anchor并过结构门的最小d=49.2236228”属于不同分母，不能替代原裁决。

本次没有新增材料采样、generator perturbation、Jacobian、JVP 或 adjoint；因此**原P0中的真正weak direction尚未建立**。要声称非平凡弱方向，至少需非零完整材料secant、合法P0端点、非平凡稳健结构变化及相应full response证据；原参数导数零、背景J的形式null space、off-prior全材料null space都不够。

## 6. Backend 与计费陷阱

可复用的原后端为 `../../A30/a30-maxwell-material-manifold/src/a30/physics.py`，其 `FrozenMaxwellAdapter` 保持同一 S/G、同一材料及照明。E0 **没有导入它**：

- `full_state` 先预约state再建立矩阵/LU/求解；不是缓存读操作。
- `cached_anchor` 读paid N3 LU后仍调用`solve_anchor(...role='cached_LU_audit')`，8 sources 会支付8个新RHS。因此零预算读取必须直接按NPZ key读取，不能调用该接口。
- `receiver_transfer` 初次构建需要32个complex adjoint RHS；之后的解析材料乘积便宜不等于初始Jacobian准备免费。
- `background/{weak,ring}.npz` 里有已付 `J(512,1152)`，但anchor是公开background chi_b，不是P0材料，也不是qualified development base。不得将它当作合法anchor的局部Fisher或weak-mode验证。
- T corrector/PC witness search在内存构建的实际endpoint Jacobian没有随state NPZ保存。本轮没有从LU重建它。A30缓存后端没有严格forward输出误差bound，只有离散数值验证范围。

## 7. 精确路径与复现

所有相对路径均以 `Gaussian/B1/b1-research-lxz` 为根：

- 原prior代码与证书：`src/prior_constrained/prior.py`；`../b1-prior-constrained/research/prior_constrained_20261015/delegated/prior/`
- paid TRAIN：`../../A30/a30-maxwell-material-manifold/inputs/a29_release/inputs/a25/train.npz`
- online y/N3/background：`../../A30/a30-maxwell-material-manifold/inputs/a29_release/results/a29/pools/<case>.npz`，只加载n3/observation/sigma/background
- N3 exact response：`../../A30/a30-maxwell-material-manifold/inputs/a29_release/inputs/a28/results/a28/exact_states/<case>.npz`，只加载chi/data
- 独立offline clean诊断：`../../A30/a30-maxwell-material-manifold/inputs/a29_release/inputs/a25/validation.npz`，只加载case_id/clean/data/sigma/条件metadata，未加载chi/current
- T：`../b1-research-t/research/research_t_20261015/qualification/`、`mechanism_v1/2/3/`
- G准确state字段路径：`../b1-prior-constrained/research/prior_constrained_20261015/delegated/g_cache/inventory.json`
- oldcorrector逐行审计：`../b1-prior-constrained/research/prior_constrained_20261015/delegated/cached_audit/E0_AUDIT.csv`、`CORRECTION_ATTRITION.csv`，G旧解释见同级`g_cache/G_CACHE_AUDIT.md`

仅数组运算的复现入口：

```bash
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 ../../.venv_nn/bin/python research/research_lxz_20261016/delegated/e0_2d/audit_cached.py
```

## 8. 本分工裁决边界

可向主线程交付两个**公开TRAIN已知材料**的P0合格校准anchor；尚无合格的原development在线重建anchor。不能将N3小MAE、noise-pass或T/G原BOX接受直接提升为原P0支持。E2是否只作明确的known-material机制台架，或应先处理合法基础重建，由主线程依据用户LXZ协议裁决。本审计不自动启动E2、不新建P1或sampler。

计算：三次0物理缓存重放累计计量wall0.8965635829372332秒；第二次增加已付P0端点pairwise secant表，第三次按主线程要求导出前两个smooth校准anchor。新states/RHS/LU/adjoints/JVPs/GPU均0，未新增material parent；TRAIN decoder重放不是新生成parent。另有缓存矩阵G/S各16source-column乘法，未计入solverRHS。文本阅读/编写时间另属orchestration，不伪装为solver wall。
