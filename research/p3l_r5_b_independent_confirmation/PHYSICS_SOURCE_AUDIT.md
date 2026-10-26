# Native96 + Probit64 真实源码审计

状态：**SOURCE_AUDITED; NEW_PHYSICS_NOT_RUN**。本任务只读真实源码、配置、历史context JSON与数组header；没有构造Target或Book，没有新Maxwell、训练、推理或材料评分。

## 冻结对象与入口

强anchor入口是 `tools/run_r4_contexts.py`，实际导入 `src/p3r4/inversion.py:reconstruct`、`src/p3r4/anchors.py:probit_optimize`、`src/p3r4/core.py:Target`。完整配置为 `research/p3l_r4_conditional_transport/configs/TARGET.json` 与 `PROTOCOL.json`。本地这些小源码/配置与R4 sibling直接逐字节比较结果见 [SOURCE_CONFIG_IDENTITIES.json](audit/independence/SOURCE_CONFIG_IDENTITIES.json)；实际历史执行仍应以每个job的SOURCE_REVISION与source_snapshot为准。

目标为原P1 scalar TM24，576 cells，8 source、32 receiver，wavenumber=2π；host χ=0.1+0.04i，ROI radius=.38，64 Gaussian DCT/offset变量与2 Uniform[0,1]属性，C2饱和混合decoder；噪声固定σ=0.00119211798514519，不按每parent SNR重标。观测32×8 complex，source-major每源Re32再Im32，white=√2/σ；Phi=||F(q)-y||²/σ²。background沿用已付费weak.npz；A30 CPU backend真实路径由 `B1_P3L_A30_ROOT` 覆盖或默认 `Gaussian/A30/a30-maxwell-material-manifold`。新实现须封存解析后实际路径，不以环境默认代替来源证据。

## Native96 的实际行为

`reconstruct(...method='full_TRF',max_evaluations=96)` 从 `[zeros64,.5,.5]` 开始，不读取truth，不以旧mode初始化。目标是 `Phi + ||z64||²/2`，两个Uniform属性只受[0,1] bounds，不加Gaussian惩罚。SciPy TRF / exact solver / x_scale='jac'；ftol=1e-9、xtol=1e-9、gtol=1e-7；96同时限制新forward状态尝试，包含失败与cache eviction后的重访。最多4个Point缓存并保留best。

此入口**没有传polish=True**，所以没有source-gradient L-BFGS polish；不能把P3L locked imaging其他入口的独立Native96+polish误当R4算法。final KKT是投影到属性box切锥的gradient，按sqrt(max(diag(JaugᵀJaug),1))缩放；它是诊断，不是额外收敛门。noise-ball 95% Chi²/2也是诊断，不是停止准则。达到cap与SciPy success分别保留，不可等同充分优化。

## Probit64 的实际行为

从Native best开始；仅初始化时将u夹至[1e-4,1-1e-4]再取ndtri。新坐标x=(z64,ndtri(u2))；目标 `Phi + ||x66||²/2`，是原Uniform先验精确变换后相对于dx的坐标mode，标识 `PROBIT_COORDINATE_MAP_NOT_NATIVE_MAP`。不能称原native q MAP，也不能把GN H=JxᵀJx+I当精确posterior协方差。

TRF/exact/x_scale='jac'，同ftol=xtol=1e-9、gtol=1e-7，max_nfev=64并独立限制fresh attempts≤64；保留current/best，完全相同native q可以复用。最终best的receiver J若缺失会补付32 RHS，无额外forward。返回full66 J、x、GN H、prediction/residual，context包括原native与probit诊断。

## 失败保留与调用方风险

Native出现普通异常或BudgetStop后，只要已有finite best便可能返回 `FAILED_BEST_RETAINED` / `BUDGET_STOP`；没有best才抛出。Probit普通优化异常也保留best，BudgetStop重新抛出；无best抛出。原 `run_r4_contexts.py` 不检查Native report status就进入Probit；若失败发生于预算边界，调用方必须在新ledger硬限额前拒绝后续动作。原runner对无finite best没有自动Native/prior fallback，整个job失败，不能默默丢掉该parent。

新确认流程若需要全部parent可评分，应在生成前封存“Probit有限best→Native有限best→原prior material mean”等明确策略、错误receipt与缺失字段处理；采用替代context会形成披露的fallback，不能伪装成同样Native96+Probit64成功。所有parent保留，失败不允许重抽truth。此处仅审计并提出要求，没有修改历史src。

## RHS 收费、历史费用与8-parent探针

`Book.reserve_state` 在A30 full_state构造矩阵/LU/solve前预约：每新状态8 forward RHS、1 LU；adapter callback消耗同一credit，不重复收费。receiver transfer解Sᴴ共32 adjoint RHS，Point内J复用不收费；source-adjoint Phi gradient若未有J则8 RHS，但此R4强入口没有polish。cache重用的历史冷启动成本必须披露，不能把40 RHS的单点forward+J称强优化anchor。

| 历史组 | parent数 | RHS中位 | RHS最大 | RHS总计 | Probit SciPy success |
|---|---:|---:|---:|---:|---:|
| TRAIN | 1024 | 3432 | 6272 | 3913008 | 673 |
| VAL | 256 | 3740 | 6176 | 1015096 | 161 |
| DEV | 128 | 3376 | 6136 | 479696 | 83 |
| LOCKED | 512 | 3376 | 6208 | 1934840 | 349 |
| CONFIRM | 256 | 3616 | 6208 | 987568 | 162 |

数据仅来自原context JSON，见 [HISTORICAL_CONTEXT_COSTS.json](audit/independence/HISTORICAL_CONTEXT_COSTS.json)；success为优化器标志，未success仍可有保留finite best。JSON中的 `native_failure_count=0` 只表示KKT诊断状态为COMPUTED，不能解读为所有Native优化器success。

保守上界：最多(96+64)×(8+32)=6,400 RHS/parent；truth生成另8 RHS/parent。512 confirmation合计最多3,280,896 RHS，未含额外机制实验/梯度/重试；实测context合计对照只能作为历史量级，不能提前报新运行实际成本。生成前独立8-parent成本probe上界51,264 RHS（8 truth + 8强context）；串行CPU，数值线程1，独立命名空间，不占用A/B confirmation；记录每parent fresh状态、LU、RHS、CPU/wall、peak RSS、Native/Probit状态及所有失败。probe只估成本，不做truth材料评分，不以结果重选模型；是否释放后续512由parent依据既有总预算决定。

旧R4 Book绑定已关闭campaign、绝对OUT、历史authorization/ledger，不能直接在新worktree运行并重开旧账。新runner必须使用新campaign ledger与明确预算；复用数学源码不等于复用关闭ledger。本次没有执行该release。

## 新确认适配与正式配置状态

对真实 sibling 原始文件直接比较 `git show ad71a86e:path` 与R4固定 `85b7bcbf39f0eeab7062ecce17531538151082bf:path`：core/anchors/inversion/prior、run_r4_contexts、TARGET、PROTOCOL全部字节一致。真实A30 `src/a30/physics.py` 与实际CPU vendor `a24/physics.py`、`a24/contracts.py`也全部等于A30固定717d3d26，见 [FIXED_REVISION_SOURCE_COMPARISON.json](audit/independence/FIXED_REVISION_SOURCE_COMPARISON.json)。仅小源码直比，没有文件hash。

新 `src/physics_confirm.py` 直接从真实R4 sibling导入原reconstruct/probit函数；Target继承原Target，仅keep_fields=False改变逐endpoint存储，persist_cache保存最终current/field/J/receiver_T而省略dense LU。live LU仍在原Point中，原优化器不会重载持久cache，因此这些存储变化本身不改变优化步。新Book是独立confirmation ledger，从未构造旧Book；原Target继承的TARGET/background仍读取原R4路径。`src/generate_observations.py` 实际分别使用manifest latent_seed与noise_seed的default_rng，符合分流要求。

已落盘 `configs/PHYSICS_RECONSTRUCTION.json` 明确Native96无polish、Probit64、原容差与失败保留；AUTHORIZATION给出3,300,000 RHS硬上限。FREEZE.json在首次源码查看时尚不存在，在随后状态snapshot时已存在（parent线程并发准备）；本审计没有执行release或物理调用。后续生成状态应以parent线程实际receipt为准。A/B本地parents文件的完整rows与原pending逐项完全一致；8个probe ID与confirmation/历史ID无交集，16个probe prior/noise seeds内部无重复且与A/B无交集。见 [PRE_RELEASE_CONFIG_STATUS.json](audit/independence/PRE_RELEASE_CONFIG_STATUS.json)。

520 parent（probe8+A256+B256）完整最坏上界3,332,160 RHS，比3.3M多32,160。历史实际成本较低不解除硬限额；runner必须保留RESOURCE_LIMITED提前停止，禁止降低保真、扩预算或删parent求完成。现有anchors入口在余额<6500时停止，是保守保护。新read_online校验keys/shape/sigma，但尚未比较target_id内容与原Target ID；建议release前补上该等值检查与finite检查。收费一致性的运行时验证不在本次只读审计范围内。
