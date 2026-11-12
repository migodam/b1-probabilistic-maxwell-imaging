# R13 复现、来源与附件入口

本文说明已经存在的文件与执行入口，不代表尚在运行的队列已完成、独立真值已评分或论文已通过科学审查。方法结论由封存后报告给出。本地静态网页构建不等同于已经公开。

## 来源与实际执行目录

- 项目沿用既有仓库 `migodam/b1-probabilistic-maxwell-imaging`，R13位于 `research/research_max_struct_r13/`；不新建仓库。
- 历史基线checkout：`ec3a37da6f3ee90beb61e753067d4f639f4118f5`。
- 独立确认算法冻结：`3a91b35605634a18a731441693a46152eeb3363e`，见 [CONFIRMATION_FREEZE.json](audit/CONFIRMATION_FREEZE.json)。
- C2实际发运来源：`689850d725ddcbd6b65152595ac0ddf2b6eff1bd`；最终公开输入runtime包记录见 [CONFIRMATION_FINAL_BUNDLE.json](audit/CONFIRMATION_FINAL_BUNDLE.json)。Git提交标识是来源记录，不是重新执行SHA-256检查。
- 用户原ZIP：`/Volumes/migodam's-external-brain/archive/user-home/downloads/MAX_STRUCT_Research_Package_20261112.zip`；解包入口 [sources/pro/MAX_STRUCT_20261112/README.md](sources/pro/MAX_STRUCT_20261112/README.md)。该包是研究来源，不能把其中协议模板当作已经运行的结果。
- 当前完整任务/恢复入口：[TASK_RESUME.md](TASK_RESUME.md) 与 [原用户任务](sources/user/R13_TASK_FROM_REFERENCED_CHAT.md)。

实际Windows目录由父线程确认：

|用途|实际目录|
|---|---|
|继承的R11完整代码、vendor、既有物理实现|`D:/AI/MAX_STRUCT_R11_DEV_20261107`|
|C1开发与表示对照|`D:/AI/MAX_STRUCT_R13_20261112`|
|C2独立确认有限队列|`D:/AI/MAX_STRUCT_R13_CONFIRM_20261112`|

Mac已有Python环境为 `/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python`，已有NumPy/SciPy/PyTorch；此工作没有安装新包。大型输入、轨迹、NPZ及ZIP保存在外盘或专门远端目录。原运行不因网页构建而重启。

## 实际source/runtime包

下表为本地已经存在的包名；它们不是尚未创建的Release下载链接。

|包名|用途|
|---|---|
|`runtime_development_v1.zip`|首批开发runtime与公开输入|
|`runtime_refinement_v2_final.zip`|C1开发修订/表示对照runtime|
|`runtime_confirmation_C2_final.zip`|C2最终发运包，109文件、68,607,527字节；只含公开在线输入|

已有完成队列的完整结果归档在 `archives/`，逐包以对应 `*_ARCHIVE_RECEIPT.json` 为准。动态公开下载清单/Release链接只有在真实产生后才进入网站；目前不补造链接。被冻结的runtime包用于重现实际执行版本；当前 `src/` 的后续评分、绘图、文档改动不自动视为旧队列运行代码。

## 公开数据与固定NN初值

`experiments/generate_confirmation.py` 创建新48对象的公开观测和分离真值目录；`experiments/prepare_fresh_nn.py` 仅对预先指定的16个公开输入调用既有R8 `FeatureOperator` 与 `load_predictor`。它不训练、不按初始物理目标筛选，也不读取材料真值。

- 新48公开入口：`data_public/fresh48/MANIFEST.json`。
- 固定16个NN初值入口：`data_public/fresh16_nn/MANIFEST.json`；原观测/算子不变，原均匀初值以 `m_init_original_uniform` 保留。
- NN参数24,737，来源为 `research/research_max_struct_r8/RESULTS/training/physics/weights.pt`；[checkpoint及计时记录](delegated/evidence/nnprep/PREPROCESS_RECEIPT.json)。已有监督训练先验是额外信息，历史训练费用未知，不能计零。
- `sealed_truth/` 与各物理模块的封存真值不进入在线输入包；只有明确授权的离线scorer读取。

生成脚本会拒绝覆盖已有目标目录。复现生成应使用新的独立目标目录并保留原版本，而不是覆盖冻结输入：

```text
python experiments/generate_confirmation.py --root NEW_R13_ROOT
python experiments/prepare_fresh_nn.py --root NEW_R13_ROOT
```

两个入口仍需完整继承代码/权重的既有相对位置。仅复制单脚本不能替代完整依赖来源。

## C1、C2与本地控制的审计/计时合同

C1开发主要使用四个P核 `0,2,4,6`，mod50辅助步审计与各方法的原生阶段审计。具体以 `configs/development_v1.json`、`development_v2.json`、`representation_v1.json` 及对应runtime为准。

C2使用八个P核 `0,2,4,6,8,10,12,14`，每BLAS一线程。统一proposal审计规则是：跨过绝对进程CPU每2秒桶后，在下一次完整材料更新可返回的位置付费物理审计；保留原生初始、阶段、终点审计。不可中断更新/线搜索/LU可跨过桶，实际时间照实保存，不补造过去的审计。corrector和direct的每次已付物理目标评价保留。`audit_every=50`不能取代C2的共同CPU时钟解释。

实际启动入口：

```text
python src/start_confirmation.py --root D:/AI/MAX_STRUCT_R13_CONFIRM_20261112 --legacy-repo D:/AI/MAX_STRUCT_R11_DEV_20261107 --development-root D:/AI/MAX_STRUCT_R13_20261112
```

它等待有限的C1依赖释放资源，顺序调用核心48、固定16、公开残差触发的条件救援和物理合同队列，再归档。当前已运行的同一目录不会重复启动；新复现使用独立解包目录和明确资源分配。没有永久daemon、调度服务或新的CUDA后端。

单队列与单任务的真实入口为：

```text
python src/finite_queue.py --root R13_ROOT --legacy-repo LEGACY_ROOT --manifest R13_ROOT/configs/confirmation_core48.json --cores 0,2,4,6,8,10,12,14
python src/run_experiment.py --root R13_ROOT --legacy-repo LEGACY_ROOT --manifest R13_ROOT/configs/confirmation_core48.json --task-id r13_fresh_0000__dr
```

独立确认配置见 `configs/confirmation_core48.json`、`configs/confirmation_mechanism16.json`、`configs/conditional_rescue_policy.json` 和 `configs/physics_v1.json`。触发数来自公开物理残差，不能事先虚构；旁路对象原样保留并不代表已运行救援算法无损。

M1上的 `src/local_queue.py` 与 `configs/development_v3_local.json` 是另列的已暴露对象控制；它的CPU时间不与i9跨平台合并作速度排名。

## 封存、离线评分和图表

`src/seal_results.py` 只读取公开预测/元数据，记录全体预期任务、失败/缺失及文件mtime/字节数，**始终写入 `authorized_scoring=false`**。父线程检查终态队列与预期任务后才可授权。新确认遵循配置的联合封存约定，不能在另一确认队列仍可更改时提前读真值。

```text
python src/seal_results.py --root R13_ROOT --manifest R13_ROOT/configs/CAMPAIGN.json --output R13_ROOT/audit/CAMPAIGN_seal.json
python src/score_results.py --root R13_ROOT --seal R13_ROOT/audit/CAMPAIGN_seal.json --out R13_ROOT/analysis/CAMPAIGN --structure-gallery
```

第二条只有在父线程显式授权seal后才执行。它再次验证mtime/字节数和终态队列，随后才读真值；不运行优化器或新物理正演。分析输出放在冻结的 `RESULTS/CAMPAIGN/`、`queues/CAMPAIGN/` 外部，保护预测快照。

绝对色标图固定为共同χ范围；`--structure-gallery`另加逐行共享GT-derived物理色标、明确vmax和超范围饱和数量，仅作离线结构展示，不改变选择或原图。虚拟无GT集成夹具位于 `delegated/evidence/generation/scorer_fixture/`，永不进入研究图集或公开报告。

费用字段保留独立口径：外部child进程CPU、内部过程CPU、并行队列历时、physical LU/complex RHS与projector作用。累计history、辅助delta和外部/内部CPU不能相加重计。中断/失败尝试保留；未知时间仍为空。Fresnel/classic诊断δ不冒称已知CN噪声证书；跨网格模型误差与噪声半径分开。

## 本地网页构建

```text
python src/build_site.py --root R13_ROOT
```

该入口只枚举实际 `analysis/<campaign>/**/FIGURES` 与已有文档，生成 `index.html`、`GALLERY.md`、`GALLERY.html` 及每组图集。没有结果文件就不生成已完成占位卡片；没有下载清单就不显示虚构Release地址。父线程另行完成图像/科学审查与既有GitHub/Pages发布验收。
