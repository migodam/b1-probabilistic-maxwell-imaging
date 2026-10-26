# 独立 parent 与 RNG 元数据审计

状态：**PASS_METADATA_DISJOINT_WITH_SCOPE_LIMITS**。这是生成前的历史元数据审计；本任务新 Maxwell RHS=0、训练=0、真值材料评分=0，没有打开任何 q/chi/noise/observation 数组内容，也没有读取 checkpoint 权重或计算文件 hash。

## 保留原 confirmation A/B

沿用 sibling `b1-p3l-r5-parallel-ab/research/p3l_r5_parallel_ab/shared/pending_confirmation_manifests/confirmation_A_256.json` 与 `confirmation_B_256.json`，不重抽 seed、不替换 ID。两者分别 master_seed=202610260500、202610260501，各256 parent；原状态均为 `SEEDS_ONLY_NOT_GENERATED_NOT_RELEASED`。原 `audit/CONFIRMATION_PREPARATION.json` 声明新 truth=0、新 observation=0、新 RHS=0；本次对历史 R5 与当前 confirmation 输出目录扫描，未找到文件名包含这些 parent ID 的 NPZ。文件名扫描不能排除另名或打包资产，故不把它当内容级“从未生成”的证明。

本次枚举18个历史 cohort、37,240个唯一 parent ID，覆盖 P1 calibration/dev/material_check，P2 train/validation/heldout，P3L development/train/validation/locked，R2 train32768/validation2048/两批locked各256，以及 R4 DEV128/L1 512/L2 256。A/B与全部历史 ID 交集均为0，A与B交集0。A/B合计1,024个 latent_seed/noise_seed整数全部互异；由其生成的128位 SeedSequence状态，与历史74,100个不同 prior/noise状态无交集。这里的 SeedSequence状态是 RNG 元数据，不是文件摘要。

完整 roster、文件路径、字节数、mtime、检查结果见 [PARENT_SEED_AUDIT.json](audit/independence/PARENT_SEED_AUDIT.json)。历史 ID 并非每个都独立研究样本：R4强上下文训练有意复用R2 train/validation，是同一 parent 的新上下文；不可把它们重复计入样本数。

## prior 与观测噪声必须分流

| 历史阶段 | prior RNG | noise RNG |
|---|---|---|
| P1 | 每组 `default_rng(group_seed)`，批量 normal(n,64) 后 uniform(n,2) | 独立 noise_seed=2026101804，批量32次重放 |
| P2 | SeedSequence([2026101901,split_code,index,0]) | 同前三项，末项1 |
| P3L | SeedSequence([2026102201,group_code,index,0]) | 同前三项，末项1 |
| R2 | SeedSequence([2026102202,split_code,index,0]) | 同前三项，末项1 |
| R4 | SeedSequence([cohort_seed,index,0]) | 同前两项，末项1 |
| 待生成R5 A/B | 原manifest latent_seed | 原manifest noise_seed |

P1批量先normal后uniform，不能用逐parent重置同一个组seed来复刻。P2/P3L/R2/R4实际实现分别在 `tools/prepare_material_fiber.py`、`generate_p3l_observations.py`、`generate_p3r2_data.py`、`generate_r4_data.py`。SeedSequence输入列表的尾部0会受内部padding影响，单纯比较文本列表不同不足以证明RNG流不同，因此这里比较生成状态。新生成实现应完整保存原uint64整数，以 `default_rng(latent_seed)`、独立 `default_rng(noise_seed)` 分别绘制；不得把uint64转float、截成uint32、或从parent字符串重新造seed。

原pending seed整数的生成代码没有在历史R5可搜索Python源码中找到；manifest本身是可复用权威输入，但从master_seed到两个uint64的算法来源仍是缺口。本审计没有臆造 derivation。

## ID 派生 K32 流不是数据集 prior/noise

冻结 `src/p3r4/risk_fusion.py:parent_seed` 使用 `SeedSequence([20261025121,*UTF8(parent)]).generate_state(1)[0]`，输出uint32，仅控制原K32模型moment Monte Carlo。它稳定地依赖字节ID，没有使用Python随机化的 `hash(parent)`。本次对历史37,240 ID及待生成512 ID重算此小型元数据：历史内部、待生成内部、两者之间都未发现碰撞，见 [ID_DERIVED_MOMENT_RNG_AUDIT.json](audit/independence/ID_DERIVED_MOMENT_RNG_AUDIT.json)。uint32仍有有限碰撞风险；不能据此认定数据集材料不同。改名会改变K32输出，故保留原ID和冻结seed规则。

## 可用资产与范围限制

本次仅检查6个历史NPZ的ZIP目录、NPY shape/dtype header，其中4个文件只另读 `parent` 字符串字段：train32768、validation2048、旧anchor train512/validation128。所有ID落在已审历史roster，均与待生成A/B无交集。233个checkpoint共650,650,241 bytes，仅stat、没有反序列化。见 [ASSET_METADATA_INVENTORY.json](audit/independence/ASSET_METADATA_INVENTORY.json)。当前新worktree缺少的ignored inputs由sibling真实资产路径识别，没有复制巨大数组。

这证明所列历史元数据与待生成计划的ID/RNG分离；没有证明所有磁盘历史内容不重复，没有比较真值材料，没有审计未声明外部训练集，也没有认证后验或科学有效性。P1/P2/P3L/R2/R4旧test与confirmation已有暴露，应继续称历史EXPOSED，不可因新worktree而称fresh。8-parent成本probe须采用独立新命名空间与seed，不能从A/B各256 parent中抽取，也不能用probe结果改选已冻结模型。
