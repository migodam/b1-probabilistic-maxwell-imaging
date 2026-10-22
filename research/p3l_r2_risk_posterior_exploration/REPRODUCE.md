# P3L-R2 复现：保存材料估计的零物理重放与完整重算边界

本页提供**已保存材料估计的重放入口**。重放读取冻结模型、线上上下文和融合动作，输出一个复数材料场；不生成新观测，不调用 Maxwell，不打开锁定真值，也不重新计算材料风险或密度评分。恢复证据、重放估计、重训模型和重新开展物理实验是不同范围。

本页编写时只核对源码、已有 JSON 回执、资产文件的存在性和本地运行环境，**没有重新执行重放、训练、fresh 评估或物理求解**。完整原始证据通过同一仓库的[本轮 Release](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/p3lr2-risk-posterior-20261023)交付；上传与匿名读取状态单独记录，不以本页声明代替实际发布回执。

## 1. 已保存的可复现证据

已有 [reproduction/REPLAY_00000.json](reproduction/REPLAY_00000.json) 记录：

- 父样本 `R2_locked_00000`，方法 `fused`。
- `physical_queries: 0`、`truth_read: false`。
- 相对保存估计的材料 RMS 差和最大绝对差均为 `0.0`。
- 冻结方法源码引用为 `ccf3987fd3910c656fa515fc076aa636fb89c47e`。

相应 [D/jobs/deterministic_material_replay/EXECUTION.json](D/jobs/deterministic_material_replay/EXECUTION.json) 保存了实际命令，[CLOSED.json](D/jobs/deterministic_material_replay/CLOSED.json) 记录该作业完成。这个单父样本回执不是全部 256 父样本或跨平台逐位一致性的验证，也不是后验校准证据。

**父样本计数：**独立训练 32,768、独立验证 2,048、首批 fresh locked 256、独立机制确认 fresh 256，总计 35,328。512 个训练锚点和 128 个验证锚点是前两组中的子集，不能另加为独立样本。DEV 历史样本属于已暴露历史资产，不能改名当作新的盲测数据。

## 2. 环境与安装范围

本次实际查询见 [ENVIRONMENT.json](ENVIRONMENT.json)：Python 3.13.13、NumPy 2.5.3、SciPy 1.18.1、Torch 2.14.0、Matplotlib 3.11.1；平台为 macOS 26.6.2 / arm64，10 个逻辑 CPU。NumPy、SciPy 和 Torch 的构建记录均指向 Apple Accelerate BLAS/LAPACK；具体 Accelerate 版本未报告。`threadpoolctl` 未安装，因此没有额外的运行时线程池枚举结果。

原运行解释器位于仓库上两级的 `.venv_nn/bin/python`；这是本机路径，不随证据包分发。仓库 `requirements.txt` 只列出 NumPy/SciPy/Matplotlib 的版本范围，**没有列出重放所需的 Torch，也没有提供完整依赖锁文件**。其他机器应自行准备 Python 3.13 环境及上述已记录版本；现有虚拟环境目录不能跨机器直接复制。依赖下载源及远程 wheel 可用性未在本次查询中验证。

以下命令假设 `PYTHON` 指向已经准备好的解释器。这里展示的安装示例没有在本次文档任务中执行：

```sh
python3.13 -m venv /absolute/new/path/p3r2-env
PYTHON=/absolute/new/path/p3r2-env/bin/python
"$PYTHON" -m pip install 'numpy==2.5.3' 'scipy==1.18.1' 'torch==2.14.0' 'matplotlib==3.11.1'
```

运行前统一设置数值线程；重放脚本还会显式将 Torch 的 intra-op / inter-op 线程都设为 1。本次查询时裸进程 Torch 默认值为 4 / 10，不能把裸进程默认值写成原实验限制：

```sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
```

全程使用 CPU。该工作流不需要 CUDA、GPU 容器、守护进程或任务调度系统。

## 3. 证据包的字节分片恢复

打包工具为 [`tools/package_p3r2_evidence.py`](../../tools/package_p3r2_evidence.py)。它用 tar 流加 gzip level 1 压缩，再将**同一 gzip 字节流**切分为 `evidence.tar.gz.part001` 等文件；分片不是独立 tar/gzip 档案。默认上限为 1,800 MiB，每片严格小于十进制 2 GB。

拿到完整分片、`receipt.json`、`INVENTORY.json` 和 `INVENTORY_ENTRIES.jsonl` 后，先确认 receipt 的 `status` 为 `COMPLETE`。失败目录与不完整分片必须保留，不能冒充完整证据。下面按 receipt 的顺序检查每片文件名和大小后流式恢复；不计算摘要值，也不额外生成一个多 GB 拼接文件：

```sh
# BUNDLE 是实际下载或本地生成的完整证据包目录。
BUNDLE=/absolute/path/to/evidence-bundle
# RESTORE 必须是新建的空目录；不要指向旧实验仓库。
RESTORE="$(mktemp -d "${TMPDIR:-/tmp}/p3r2-evidence.XXXXXX")"
set -o pipefail
python3 - "$BUNDLE" <<'PY' | tar -xzf - -C "$RESTORE"
from pathlib import Path
import json, shutil, sys
bundle = Path(sys.argv[1])
receipt = json.loads((bundle / 'receipt.json').read_text())
if receipt.get('status') != 'COMPLETE' or not receipt.get('parts'):
    raise SystemExit('Complete, nonempty receipt required')
for row in receipt['parts']:
    name = row['file']
    if Path(name).name != name:
        raise SystemExit('Unsafe part filename')
    path = bundle / name
    if path.stat().st_size != row['size']:
        raise SystemExit('Part size differs from receipt: ' + name)
for row in receipt['parts']:
    with (bundle / row['file']).open('rb') as source:
        shutil.copyfileobj(source, sys.stdout.buffer, length=1024*1024)
PY
printf '%s\n' "$RESTORE"
```

若管线失败，保留该恢复目录作诊断，换一个新的空目录修正后再恢复。正常恢复后，用打包目录中的逐行清单核对所有归档普通文件的大小；这只核对恢复路径与字节数，不证明内容同一性：

```sh
python3 - "$BUNDLE" "$RESTORE" <<'PY'
from pathlib import Path
import json, sys
bundle, restored = map(Path, sys.argv[1:])
count = 0
with (bundle / 'INVENTORY_ENTRIES.jsonl').open() as stream:
    for line in stream:
        row = json.loads(line)
        if row['type'] != 'file' or not row['archived']:
            continue
        rel = Path(row['path'])
        if rel.is_absolute() or '..' in rel.parts:
            raise SystemExit('Unsafe inventory path')
        path = restored / rel
        if not path.is_file() or path.stat().st_size != row['size']:
            raise SystemExit('Restored size/path mismatch: ' + str(rel))
        count += 1
print('Restored regular-file paths and sizes matched:', count)
PY
```

包内路径相对仓库根保存，恢复后会出现 `src/p3r2/`、`src/p3r2_mechanism/`、`src/p1_material/`、相关 `tools/*p3r2*.py` 及 `research/p3l_r2_risk_posterior_exploration/`。`receipt.source_revision` 表示打包时 Git HEAD；包内也保留当时工作区的源码字节，因此它不构成“工作区无未提交改动”的声明。冻结方法自己的 `source_revision` 仍以 `METHODS_FROZEN_FOR_LOCKED.json` 为准。

锁、临时文件、缓存和生成的压缩包不进入证据；研究失败记录及新产生的原始数组保留。符号链接只在清单中保留 `link_target` 来源，不跟随、不恢复链接目标。**`RESOURCE_LEDGER.jsonl` 与 `RESOURCE_LEDGER.json` 被明确排除**，主线程必须在记账打包作业关闭后单独交付最终账本；公开的 sequence/时间戳资源快照不能代替最终账本与封账记录。

## 4. 零物理材料重放入口

需要恢复出的源码树和同一证据包中的下列资产：

| 输入 | 用途 |
|---|---|
| `METHODS_FROZEN_FOR_LOCKED.json` | 固定父样本集合、模型路径和版本引用 |
| `configs/TARGET.json` | 原目标材料配置 |
| `L/contexts/R2_locked_00000.npz` | 已支付的线上 `observation/q0/chi0/J/prediction/residual`；不是锁定真值 |
| 冻结模型 `BEST.pt` | direct / residual 的已训练参数及模型内保存的归一化、几何信息 |
| `D/fusion/FROZEN.json` | fused 的冻结模型选择、投影和融合系数 |
| `L/estimates/R2_locked_00000.npz` | 可选的保存估计比对，不含本次重放所需的真值 |

脚本 [`tools/replay_p3r2_material.py`](../../tools/replay_p3r2_material.py) 支持 `fused`、`direct`、`residual`。当前实现即使选择 direct，也先加载同一份保存的线上上下文；不能仅凭一个在线观测文件替代全部上下文。模型通过 Torch 加载本研究保存的 checkpoint；无需外部 A30 求解器或历史 DEV 链接内容。模型、上下文或融合动作缺失时，现有入口不能凭空重建它们。

恢复后将源码根设为 `REPO`，将 `PYTHON` 设为已准备好的解释器。`--bundle-root` 是恢复后的**研究输出目录**，不是解压根本身：

```sh
REPO="$RESTORE"
CAMPAIGN="$REPO/research/p3l_r2_risk_posterior_exploration"
# 本机原环境可用：PYTHON=/absolute/path/to/Gaussian/.venv_nn/bin/python
REPLAY_OUT="$(mktemp -d "${TMPDIR:-/tmp}/p3r2-material-replay.XXXXXX")"
"$PYTHON" "$REPO/tools/replay_p3r2_material.py" \
  --bundle-root "$CAMPAIGN" \
  --parent R2_locked_00000 --method fused \
  --output "$REPLAY_OUT/fused_00000.npz" --compare-saved
```

输出为材料数组 NPZ 与同名 JSON 元数据。脚本拒绝覆盖已有 NPZ，因此每次使用新输出路径；JSON sidecar 也应保持同一新前缀。其他方法可分别换成 `--method direct` 或 `--method residual`，并使用不同输出名。去掉 `--compare-saved` 可以不读取保存的估计文件，但仍需要冻结模型与线上上下文。

重放脚本会将包含原研究输出标记的绝对模型路径映射到恢复后的 `CAMPAIGN` 内，并拒绝逃出该目录的依赖。不能将模型、上下文、融合系数或目标配置改成别的版本后仍称为同一次冻结方法重放。

这个入口不读取 `inputs/offline_only`，不会启动 `Target`、Book 或 Maxwell，不新增 LU/RHS。输出是监督材料风险估计，不能据此改称完整后验均值、后验样本、材料图像概率或新鲜盲测结果。旧的 locked 父样本在已有评分完成后已被暴露，重放不能恢复其新鲜性。

## 5. 历史资产和外部依赖边界

[源码来源](SOURCE_INVENTORY.json) 与 [历史复用记录](AUDIT/HISTORY_REUSE.json) 保存以下出处；版本号是已有记录中的声明，本次未重新验证外部仓库提交或远程下载地址：

| 资产 | 已记录来源 | 对复现的影响 |
|---|---|---|
| DEV_000–DEV_031 观测、历史锚点缓存 | 相邻 `b1-p3l-local-posterior/research/p3l_local_posterior_imaging`，历史提交 `c1171ab9c9006f5061f10293c7eadd5b977dd608` | 现有输出里是绝对符号链接；证据包不打入目标内容。重做涉及这些历史缓存的 B/C 审计，需要另取原资产并保留其来源 |
| `inputs/background/weak.npz` | 既有付费 P1 G/S/incident 背景数组 | 本研究中是普通文件，会随 OUT 打包；从零生成背景不属于免费复现 |
| A30 Maxwell adapter 源码 | 声明冻结提交 `717d3d26a290248bd2b9aeb55f3d03694b407b84` | 不在此窄范围证据包的源码目录中；新物理重算需要另备 |
| P1 decoder / geometry | `src/p1_material/` 与 `research/p1_material_target_20261018/inputs/geometry.json` | 窄范围证据包显式包含；训练和原目标重建不可换成不同网格/解码器 |
| Python 运行环境 | 本机 `.venv_nn`；历史 `configs/RUNTIME_ENVIRONMENT.json` 与本次 ENVIRONMENT.json | 虚拟环境不打包；需要另行安装兼容依赖 |

A30 的默认本地位置为 `Gaussian/A30/a30-maxwell-material-manifold`，代码也接受 `B1_P3L_A30_ROOT`。只有重新计算物理端点才需要它；保存材料重放不需要设置这个变量。历史记录的 8,374,472 RHS 属于旧封账，不能作为新实验的余额。

## 6. 哪些入口会开始新的工作

| 工具 | 实际操作与边界 |
|---|---|
| `prepare_p3r2.py` | 写授权和协议、建立历史链接。它不是获得新授权的替代品；不能在旧输出上重跑 |
| `generate_p3r2_data.py` | 调用付费物理生成独立观测、训练/离线标签；locked 生成要求先冻结方法 |
| `run_p3r2_anchors.py` | 支付 native/probit 重建和上下文 J；fresh 使用 native96 + probit64，matched 控制增加16个 fresh-state 上限 |
| `pack_p3r2_anchors.py` | 组合已经保存的训练/验证标签与锚点；只接受 train/validation，不能用来打开 locked 标签 |
| `train_p3r2_risk.py`、`train_p3r2_residual.py`、`train_p3r2_density.py` | 从训练标签重训并在独立验证集选择；无 Maxwell 也有 CPU/内存成本，不是保存模型的确定性重放 |
| `run_p3r2_fresh_evaluation.py` | 一次性顺序执行新 locked 生成、物理锚点、轨迹、前向诊断与最终离线评分；不是查看既有结果的入口 |
| `evaluate_p3r2_locked.py`、`evaluate_p3r2_fresh_trajectories.py` | 所有门控通过后才打开 locked 真值；本次复现文档任务没有运行它们 |
| `collect_p3r2_public_results.py`、`plot_p3r2_results.py` | 只整理已评分结果；要求两个 exact256 评估 COMPLETE，不以此替代评估 |
| `package_p3r2_evidence.py` | 只做字节证据打包；必须等其他作业关闭，拒绝覆盖旧证据 |

Book 的权威记录是 `RESOURCE_LEDGER.jsonl` 的最后一条 `after`，JSON snapshot 可恢复但不应替代权威流水。`active_wall` 是尚未关闭的预算预留；失败、中断和拒绝端点的成本及原始输出必须保留。相同已关闭 job 不允许重用，未核销的活动 job 不会自动重放。出现 `CAMPAIGN_CLOSED.json` 后，Book 拒绝在旧研究继续作业；`analysis_only=True` 的 Book 也不能构造物理 Target。

物理进程使用 `paid_worker.lock` 串行运行，receiver transfer 缓存最多 4 个条目。数据生成/锚点脚本在指定检查边界检查 RSS > 2,000 MiB；这不是操作系统内存硬隔离。`analysis_only` 包装器不持物理 worker 锁，主线程仍需防止多个分析进程同时积累大型数组。不要把这些源码策略误写成所有任务的实测峰值。

**真正重算数据、模型或物理证据，需要新的明确授权、预算与新的输出根／工作副本。**现有生成器、Book 和若干工具从 `p3r2.core.OUT` 固定定位输出，未统一提供可切换新实验根的 CLI；因此不能只加一个不存在的 `--root` 就称为独立重算。应由新的授权工作在隔离副本里显式配置新的 campaign ID、输出根和独立账本，复核目标配置与种子后执行。不要删除、清零、改名旧 ledger，不要移除 closure 来绕过门控，不要覆盖原模型或冻结配置。CPU-only 重训也应使用新的输出与独立记账。

对打包作业，主线程在其他活动作业关闭后，用 D 阶段分析包装器登记唯一打包 job，并将同一名字传给 `--accounted-job`。工具核验它是唯一活动的当前 D/analysis-only JOB_START，持有预算锁冻结新预留，并持物理 worker 锁；打包退出后由包装器记完 JOB_END，再另行保存最终账本。若旧 campaign 已封账，不应为打包删除封账标记；后续交付应由主线程按新的交付授权安排。

## 7. 本交付的依赖缺口

1. 零物理材料重放入口和一个已保存的成功回执存在；本次没有重新执行或扩展到全256。
2. ENVIRONMENT.json 保存当前解释器、版本、构建/BLAS 与线程默认值；完整依赖锁文件、跨平台误差保证未提供。
3. 完整复原历史 B/C 路径需要证据包外的 DEV 历史链接目标；物理从零复算还需要外部 A30 源码、其依赖及新的付费授权。
4. 证据包及最终账本是否已对外可下载，应以后续实际发布入口和主线程读回为准。这里没有执行打包、上传或 Release 检查。

## 8. 独立机制确认的附加源码

`src/p3r2_mechanism/` 保存第二批256对象的7家族21模型协议。9个新增训练运行与首批77个工件分开；随机decoder包络使用固定QR seed 2026102303，不根据测试结果重选。`MECHANISM_CONFIRMATION_FROZEN.json` 与源提交5d5332a1在第一份新观测前冻结。两批结果不可混成事前512人检验；第二批没有执行matched16。这里不提供绕过新授权和旧封账的重跑入口。

## 9. 已封账交付及可选依赖覆盖包

本轮已封账：5,956,328新增complex RHS、3.86296588小时受控CPU墙钟、GPU0。完整原始档为51个gzip字节分卷，receipt状态COMPLETE；最终ledger明确另交。恢复全部分卷后，覆盖FINALIZATION.zip，并把独立RESOURCE_LEDGER.jsonl.gz解压到本研究目录。不得移除CAMPAIGN_CLOSED标记来复用余额。网络上传与ledger压缩I/O有独立publication回执，不算作新Maxwell。

同一[Release](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/p3lr2-risk-posterior-20261023)另提供两个可选包，补齐主原始档故意未跟随的外部依赖：

- A30_FROZEN_SOURCE.zip：固定717d3d26a290248bd2b9aeb55f3d03694b407b84的既有src/a30及原许可文件；无需下载整个外部仓库。解压后新授权物理重算可设B1_P3L_A30_ROOT为A30_FROZEN。仅保存材料重放不需要它。
- HISTORICAL_DEV_PREREQUISITES.zip：32个已暴露DEV的80个原本符号链接目标，约176.6MB，作为同一新恢复目录的普通文件覆盖层。只复制既有字节、未解读数组、没有新增物理，也不含旧LOCKED文件。offline_only标签继续只供离线评价，不能用于online提案。这些是已付历史证据，不是新测试对象或新teacher。

pro_handoff.zip约22.25MB，是面向GPT Pro的精选理论与实际证据；不需要下载约95GB完整分卷便可先研究。其生成概要名为HANDOFF_README.md，另保留原仓库README.md，无重复成员。最初存在README同名的本地V1包已保留并标为交付替代，公开只采用V2。

原始分卷包含封账前报告快照；以本网页、最终账本、封账覆盖层及公开发布回执为最新交付状态。可用性与科学资格分开：FULL posterior始终未合格，不因下载或重放成功而改变。

## 10. A30源包的必要vendored依赖补充

A30适配层会导入`vendor/frozen_a29/src/a24/physics.py`，因此仅`A30_FROZEN_SOURCE.zip`不足以执行新物理重算。还需将同一固定commit的`A30_FROZEN_VENDOR.zip`解到相同新目录，合成`A30_FROZEN/src/a30`与`A30_FROZEN/vendor/frozen_a29/src`，保留来源和许可。使用新的Python进程避免canonical包版本碰撞。保存材料的零物理重放不受影响。原接口源包及第一次导出回执保留，新V2回执和[依赖补充](DEPENDENCY_ADDENDUM.md)明确纠正该缺口；没有修改solver、没有新LU/RHS，也没有重新运行历史实验。
