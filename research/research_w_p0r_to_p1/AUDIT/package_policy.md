# W 结束后证据打包政策

当前交付仅为 `tools/package_w_evidence.py` 与本文；没有执行打包器，没有创建 `CLOSED.json`、inventory、ZIP 或 Release，没有发布、提交或推送。实现只使用 Python 3.10+ 标准库；本轮仅做语法/静态检查，尚未执行真实打包验证。

## 必须先由父线程关闭研究

打包器的第一项实际操作是检查阶段根目录中已存在的普通文件 `CLOSED.json`。没有该文件，完整打包和 `--inventory-only` 都退出；工具不会创建该标记，不推断任何科学 GO，也不修改源文件、配置、输入、数据或账本。关闭标记由父线程在获授权的结束时写入，打包器不负责关闭研究。

默认阶段 ROOT 为本脚本上级目录，即：

`/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-w-p0r-to-p1/research/research_w_p0r_to_p1`

只在已有 `CLOSED.json` 后，在工作树目录运行以下命令。两次输出使用不同新目录，禁止覆盖前一份证据：

```sh
python3 research/research_w_p0r_to_p1/tools/package_w_evidence.py \
  --inventory-only \
  --output research/research_w_p0r_to_p1/release_evidence/final_inventory

python3 research/research_w_p0r_to_p1/tools/package_w_evidence.py \
  --output research/research_w_p0r_to_p1/release_evidence/final_package
```

默认固定历史代码目录为相邻 `research/research_w_continuation`。默认历史本地原始证据目录为 `/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/Research_W/research-w-continuation`。位置变化可分别指定 `--historical-root` 和 `--historical-local-root`；缺少历史源码/vendor 或历史本地 results 时明确失败，不默默省略。`--root` 可显式指定阶段目录；依旧必须有该目录自己的关闭标记。

## 收录规则：按用途，不按“小于某个字节数”

完整 inventory 对所有扫描到的普通文件给出 origin、相对路径、归档路径、字节数、角色、收录与否、原因及本地位置。运行环境、缓存、既有打包输出等目录不递归，记录被排除的目录项与原因。符号链接和特殊文件不跟随或复制。

白名单包括：

- 本轮 `src/tools/tests`、全部配置与公开 inputs；固定历史 `src/vendor/tools/tests`，包括 FEM、mesh、旧搜索和 vendored p0 几何源码。外部历史目录的可变源码不替代工作树中固定的可重放版本。
- 阶段和固定历史的根报告、HTML、JSON/JSONL、CSV/TSV、ledger、失败日志、图表、审计、文档与运行源码快照。项目和历史中的 requirements/版本声明/依赖锁定元数据一并保留；不包含系统 Python、安装包、`.deps`、虚拟环境、缓存或 node_modules。
- 所有白名单目录内的 response、latent、material、iteration、checkpoint、q、mean、prediction、Jacobian 和诊断数组。不设置小尺寸上限；重要数组即使大于单 asset 也收录，必要时以不解码的原始字节分块。
- P1A `rules` 和 P1B public full-bank/rule aggregate 的全部现存文件及伴随元数据，不以尺寸筛掉 aggregate。对协议列出的完整/截断 P1A 规则，另列 `missing_locally`；缺失不会被制造或写成已完成。此检查是完整性记录，不能把提前关闭的未运行阶段表述为运行成功。
- P1B inference 的整个输出子树，包括每 seed 的原始 trajectory、proposal library、finite reference、诊断和失败记录；offline/offline_scores 完整保留。其余 chain/particle 数组也由类型/角色收录，不只收最终图表。

已知完整原始场采用文件名规则区别：当前 `fN.npz/npy`、历史 `WORLD_*.npz`、明确的 `*_FIELD`/full-field 归档等。已明确为 mean/prediction/q/checkpoint/诊断的数组优先按证据用途保留。一般的 `ARRAYS.npz` 不能凭名字或大小排除，因为历史搜索和响应诊断依赖它。

存在一个由旧保存源码确认的例外：历史 `results/p0b/finer-endpoints/ARRAYS.npz` 混合保存 `*_data`、`*_points`、`*_epsilon` 和 `*_field`。原文件整体留在本地并列入 full-fields inventory。结束后打包时，只对这一明确命名的非 GT 归档，按 ZIP 成员名筛除 `*_field.npy`/`*_dofs.npy`，把其余 response/material/diagnostic NPY 成员流式复制到 `_NONFIELD_DERIVED/.../ARRAYS.npz`。派生件注明原始来源和每个成员的包含/排除原因；不读取 NPY 头或数组值，不解压被排除的场。原始文件不变。inventory-only 只读该归档的 ZIP 目录成员名/尺寸，不读其成员 payload，也不生成派生 NPZ，并明确把该派生件标为尚未创建。

## 公开 full-field subset 的严格边界

选择由 `configs/FIRST_JOINT_NOMINATION.json` 冻结的第一个 joint 见证，要求其独立 validation `RESULTS.json` 的两个既有标记 `numerical_pass=true` 与 `joint_physics_threshold=true`，并要求 validation MANIFEST 与 nomination 一致。不会按最好残差、文件排序或最新完成者另选见证。

当前实际目录和文件名/字节元数据已确认对应：

`results/B2/validation/joint_seed2_narrow_fullj/fields/`

其四个完整公开端点为：

- `cal_world0_fine/endpoint_0000/`
- `cal_world1_fine/endpoint_0000/`
- `joint8_world0_fine/endpoint_0000/`
- `joint8_world1_fine/endpoint_0000/`

每端点完整保留 `f0.npz` 至 `f4.npz`，共 **20 个频率归档**，以及 `PARAMETERS.json`、`FREQUENCIES.json`、`RESPONSE.npz`、`COMPLETE.json` 和端点中其余伴随 material/mesh 文件。频率归档整文件复制，不抽取、降采样或重算场、材料与网格。此次仅检查文件名和字节，未打开任何 NPZ；工具也不会把归档存在性表述为对内部数组的独立物理验证。

其余历史及本轮独立 raw complex-field 归档保留原本地路径，不删除、不搬走、不覆盖。`LOCAL_FULL_FIELDS_INVENTORY.json` 明示每项路径、字节、角色、是否公开及原因；`full_fields_all_uploaded=false`。公开页面/Release 应写“响应、检查点和诊断完整；原始场为上述见证 subset”，不得写“所有 full fields 都已上传”。历史混合归档的非场派生件仅保留响应/材料/诊断，不扩张原始场 subset。

## GT 与生成计划只在关闭后作为不透明文件证据复制

本轮开发没有读取 GT 数组或 GENERATION_PLAN 内容。打包器也不导入 numpy、不打开任何真值数组、不解析生成计划。结束后，仅允许 designated `results/P1A/data/` 与 `results/P1B/data/` 下的 data 证据按完整文件字节复制，包括已有 OFFLINE_TRUTH/OFFLINE_GENERATION_PLAN。复制过程是二进制流，不是数组加载、评分、推断或物理生成。相关路径的角色为 `POST_CLOSE_SEALED_DATA_OPAQUE_COPY`。

这些受保护文件若位于上述 data 目录以外则列入排除项，供父线程明确处理，不自动暴露。ordinary 报告中含“truth”的标题不等同于真值数组，仍按报告白名单收录。inventory-only 不打开任何 data/GT payload。

## ZIP、inventory 与恢复

单个 ZIP asset 上限为十进制 **1,800,000,000 字节**；`--max-asset-bytes` 只允许降低到不小于 1,000,000 字节。使用 `ZIP_STORED`，不反复压缩已有 NPZ，也不将多个 ZIP 再包成一个超大文件。每个 ZIP 都是独立可解压的标准归档，命名 `w-evidence-part001.zip` 等。

一个必要原始文件超过可用单 asset 容量时，按字节切为 `__chunks__/...partNNNNNN`；`REASSEMBLY.json` 记录目标相对路径、原始字节数、offset、长度、对应 asset 和成员名。恢复时解压全部 ZIP，然后按 offset 顺序二进制拼接；无分块文件已保留原相对路径。`_NONFIELD_DERIVED` 是明确标记的派生 archive，不冒充完整原始 mixed archive。

输出文件包括：

- `INVENTORY.json`、`INVENTORY.csv`：完整选择与排除记录。
- `LOCAL_FULL_FIELDS_INVENTORY.json`：本地场保留与公开 subset 边界。
- `PACKAGE_MANIFEST.json`、`REASSEMBLY.json`、`RESTORE.md`：版本、来源、分块恢复说明。
- `ASSETS.json`：inventory-only 时明确 `NO_ZIP_CREATED`；真实打包成功后列出各 ZIP 的实际字节数。
- 如复制中失败：`PACKAGE_FAILURE.json` 标为 `PACKAGE_INCOMPLETE_DO_NOT_PUBLISH`，保留 partial 供检查；不得将已写出的部分 ZIP 当作完整交付。

工具不生成或检查 hash/SHA，不附加内容指纹。复制前后核对源文件 size/mtime/inode/device；源文件变动时明确失败。该检查不等于内容哈希验证，不宣称密码学完整性。打包不触发 Maxwell、FEM、LU、RHS、Book、任务队列或推断。

全部 ZIP 与独立 metadata 文件应上传到**现有仓库的同一个 Release、多 asset**；没有创建新仓库的代码，也没有调用发布 API/CLI。实际上传、网页链接及最终科学状态仍由父线程处理，本工具只生成本地可检查证据。
