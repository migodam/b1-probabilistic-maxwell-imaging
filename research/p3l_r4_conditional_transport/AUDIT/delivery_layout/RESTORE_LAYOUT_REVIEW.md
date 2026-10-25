# R4 全量证据恢复布局审阅

本审阅只检查打包元数据和恢复约定，没有解包 175 GB 归档、读取科学数组、运行模型或增加科学结论。打包 `RECEIPT.json` 为 `COMPLETE`，源提交 `60066394392fe04fb6a5b30371ffa34c978ade9d`，封存研究账本序号 `1522869`。提交号是来源记录；归档包含当时实际工作树字节，不等于干净 Git checkout。上传是否完成应由最终上传回执和远端审计单独判断。

## 六个类别与恢复位置

所有类别必须解到同一个全新目录；它们共享归档相对路径，并不是六份可替换的数据集。类别只是储存分类，不改变历史复用、新样本或后验资格判断。

| 类别 | 编码 | 分片数 | 分片字节合计 | 精确范围 |
|---|---|---:|---:|---|
| compact_core | gzip1 | 1 | 1,323,321,353 | `compact_core.tar.gz.part0001` 至 `compact_core.tar.gz.part0001` |
| raw_physics_and_arrays | tar | 113 | 169,457,408,000 | `raw_physics_and_arrays.tar.part0001` 至 `raw_physics_and_arrays.tar.part0113` |
| raw_neural_and_scores | tar | 3 | 3,071,047,680 | `raw_neural_and_scores.tar.part0001` 至 `raw_neural_and_scores.tar.part0003` |
| reused_inputs | tar | 1 | 227,850,240 | `reused_inputs.tar.part0001` 至 `reused_inputs.tar.part0001` |
| derived_training_inputs | tar | 1 | 1,208,289,280 | `derived_training_inputs.tar.part0001` 至 `derived_training_inputs.tar.part0001` |
| new_R4_inputs | tar | 1 | 25,159,680 | `new_R4_inputs.tar.part0001` 至 `new_R4_inputs.tar.part0001` |

总计 **120 个分片，175,313,076,233 字节**。`compact_core` 是 gzip 压缩的一个连续 tar 流，其余类别是未压缩 tar 流。分片是连续流的切片，不能逐片单独解压。

恢复后的顶层布局为：

```text
NEW-r4-restore/
  workspace/
    src/
    tools/
    tests/
    research/p3l_r4_conditional_transport/
    research/p1_material_target_20261018/inputs/geometry.json
    requirements.txt
    ...当时选入的根目录文档
  dependencies/
    a30/src/
    a30/vendor/
    a30/requirements.txt
    a30/a30_config.json
    r2_baseline/local_geometry_cnn_s33/
  bundle_metadata/
    SUMMARY.json
    DEPENDENCIES.json
    EXCLUSIONS_AND_SYMLINKS.json
    entries/*.jsonl
```

`bundle_metadata/` 由 compact_core 流附带；完整 pack 最终回执本身在打包结束后生成，需要作为独立小附件保存。

## 恢复命令：拼接后直接流式解包

将 120 个分片和本目录的 `R4_FULL_EVIDENCE_PARTS.json` 放入同一下载目录。以下命令使用 Bash；先将两条示例路径改成实际绝对路径。`RESTORE` 必须尚不存在，避免覆盖现有工作。预检仅比较精确名称集合和字节数，不读取内容或计算散列。此文档只提供命令，本次没有执行这些解包命令。

```bash
set -euo pipefail
export LC_ALL=C
PARTS="/absolute/path/to/downloaded-r4-parts"
RESTORE="/absolute/path/to/NEW-r4-restore"

python3 - "$PARTS" <<'CHECK_PARTS'
import json, pathlib, sys
p = pathlib.Path(sys.argv[1])
m = json.loads((p / "R4_FULL_EVIDENCE_PARTS.json").read_text())
expected = {v["file"]: v["bytes"] for v in m["parts"]}
actual = {v.name for v in p.glob("*.tar*.part[0-9][0-9][0-9][0-9]")}
if actual != set(expected):
    raise SystemExit(f"Part names differ: missing={sorted(set(expected)-actual)}, extra={sorted(actual-set(expected))}")
for name, size in expected.items():
    f = p / name
    if f.is_symlink() or not f.is_file() or f.stat().st_size != size:
        raise SystemExit(f"Part type/size mismatch: {name}")
print(f"NAME_SIZE_MATCH: {len(expected)} parts / {sum(expected.values())} bytes")
CHECK_PARTS

mkdir "$RESTORE"
cat "$PARTS"/compact_core.tar.gz.part[0-9][0-9][0-9][0-9] | tar -xzpf - -C "$RESTORE"
cat "$PARTS"/raw_physics_and_arrays.tar.part[0-9][0-9][0-9][0-9] | tar -xpf - -C "$RESTORE"
cat "$PARTS"/raw_neural_and_scores.tar.part[0-9][0-9][0-9][0-9] | tar -xpf - -C "$RESTORE"
cat "$PARTS"/reused_inputs.tar.part[0-9][0-9][0-9][0-9] | tar -xpf - -C "$RESTORE"
cat "$PARTS"/derived_training_inputs.tar.part[0-9][0-9][0-9][0-9] | tar -xpf - -C "$RESTORE"
cat "$PARTS"/new_R4_inputs.tar.part[0-9][0-9][0-9][0-9] | tar -xpf - -C "$RESTORE"

python3 "$RESTORE/workspace/tools/package_r4_evidence.py" verify-restored \
  --inventory "$RESTORE/bundle_metadata" \
  --restored "$RESTORE" \
  --restore-times
```

`LC_ALL=C` 与四位数字名称给出正确的分片顺序。此流程不另写一份拼接后的巨型 tar。`verify-restored` 只走已存清单，检查条目类型、文件大小和纳秒 mtime；`--restore-times` 只修改新解包目录的时间，以兼容 tar/文件系统的时间精度行为。它不做内容散列或科学复算，成功也不能表述为科学结果复现。目标文件系统仍需支持相应时间精度。

## 依赖与环境设置

解包与清单检查需要 Python 3 标准库、支持 PAX 的 tar、gzip 和 Bash。打包机器为 `macOS-26.6.2-arm64-arm-64bit-Mach-O`，10 个逻辑 CPU、24 GiB 内存，数值线程数 1，MPS 可用，CUDA 未使用。该记录来自归档的 `workspace/research/p3l_r4_conditional_transport/ENVIRONMENT.json`，不是跨平台兼容性证明。

如需之后单独研究运行环境，可在解包目录之外建独立虚拟环境，并先审阅 `workspace/requirements.txt` 与 `dependencies/a30/requirements.txt`。根 requirements 声明 `numpy>=2.2,<3`、`scipy>=1.15,<2`、`matplotlib>=3.9,<4`；它不是完整环境锁文件，未固定 Python、PyTorch、OS 或 MPS 版本。Python/NumPy/SciPy/PyTorch/MPS 和操作系统均未打包。根据归档中的运行记录选择版本，并把任何重建环境记录成独立补充；本次不安装依赖，也不启动新的训练、推断或物理计算。

捆绑依赖包括实际 P1 runtime、`p3r2`/`p3r2_mechanism` runtime、A30 `src` 与完整 `vendor` 树、冻结的 R2 `local_geometry_cnn_s33` baseline 和 P1 geometry。A30 声明与实际提交均为 `717d3d26a290248bd2b9aeb55f3d03694b407b84`。R2 声明 pin 为 `6eadb4690c962aaa95821fc4056b677f0f924a6b`，打包时 R2 工作树提交为 `33471215c4e3005a9b810d5375e3ea2942d4d71b`，二者在依赖记录中分别保留，不应合并为同一身份。

## 最重要的三项可移植性限制

1. **历史绝对路径。** 冻结清单保留原始 `/Volumes/migodam's-external-brain/...` 路径与 mtime 身份。`workspace/` 与 `dependencies/` 是交付布局，不自动等于运行时代码预期的同级 A30/R2 工作树。需要恢复原挂载布局，或另写明确的路径迁移映射；不能静默改写封存清单。完整历史 R2 cache 未递归打包。本次清单报告 0 个 symlink，排除项仍需按 `EXCLUSIONS_AND_SYMLINKS.json` 审阅。
2. **打包后的补充另存。** 封存快照以源提交、账本序号与 pack COMPLETE 回执为边界；后来的交付脚本、恢复说明、上传重试/清理/最终远端审计和追加出版说明属于单独的 post-pack supplement。本说明也属于补充。保留原始分片，另附补充清单和独立来源/时间记录，不把后来的文件冒充成旧归档内已有内容。
3. **空间与条目数量。** 分片占 175.313 GB；原始逻辑文件内容为 178.237 GB，约 166 万文件及 74,675 个目录，库存 JSONL 元数据另约 1.114 GB。保留下载分片再流式解包，逻辑下限约 354.664 GB，尚未计入文件系统分配、目录/元数据与余量；建议约 400 GB 或更多可用空间，并确认 inode/文件数量能力。此建议是容量规划，不是本次测得的恢复峰值。避免额外落盘一份拼接 tar，否则空间还会显著增加。

## 最终 Release 需要的小附件

| 建议附件名 | 来源/状态 | 用途 |
|---|---|---|
| `R4_FULL_EVIDENCE_PARTS.json` | 本目录已生成 | 120 分片精确名称、字节数、类别、编码、pack 来源 |
| `R4_FULL_EVIDENCE_PARTS.tsv` | 本目录已生成 | 便于人读、表格与下载核对 |
| `RESTORE_LAYOUT_REVIEW.md` | 本文 | 新目录恢复命令、依赖与边界 |
| `R4_PACK_RECEIPT.json` | 原 pack `RECEIPT.json` 的原样副本，待主线程发布 | 最终 COMPLETE 状态与完整 parts map |
| `R4_DEPENDENCIES.json` | 原 inventory `DEPENDENCIES.json` 的原样副本，待主线程发布 | 来源提交、捆绑与外部依赖 |
| `R4_INVENTORY_SUMMARY.json` | 原 inventory `SUMMARY.json` 的原样副本，待主线程发布 | 类别与恢复条目规模；其内部状态是库存阶段，不替代最终 pack receipt |
| `R4_PACK_DELIVERY_RECEIPT.json` | 原 pack `DELIVERY_RECEIPT.json` 的原样副本，待主线程发布 | 打包成本、0 RHS、释放锁记录 |
| 最终上传 receipt 与远端 asset audit | 上传完成后由主线程生成/选定 | 证明完整远端名称、字节数及 uploaded 状态 |
| post-pack supplement manifest | 由主线程冻结后单独提供 | 明确后补文件的来源与旧归档边界 |

本审阅未上传任何附件；没有核验最终公开链接或匿名可访问性。小附件可以使用上表唯一名称，避免多个模糊的 `RECEIPT.json` 在 Release 页面上混淆。
