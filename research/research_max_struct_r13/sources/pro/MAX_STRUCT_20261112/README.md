# MAX-STRUCT 研究交付包

本包记录2026-11-12新加坡时间完成的独立理论研究、源码审查与真实Maxwell机制实验。入口是 `report.md`；19页可打印版本为包内 `MAX_STRUCT_Research_20261112.pdf`，并另行提供下载。所有新结果都是已暴露对象上的机制开发证据。

## 阅读顺序

1. `report.md`：中文综合报告、结论与边界。
2. `theory_audit.md`：完整命题、证明、反例和源—接收KKT联系。
3. `evidence_audit.md`：11起点/44运行/10候选池的独立审查。
4. `literature_review.md`、`bibliography.tsv`：35条主来源、阅读深度、原创性冲突。
5. `paper_experiment_protocol.md`、`CODEX_TASK.md`：待执行的论文级工作。

## 已完成计算

`work/normal_memory/` 包含三个实例的原始JSON、冻结场/电流/材料NPZ、普通/分数/整数分支、离线评分和完整轨迹。`work/evidence_analysis/` 包含原配对结果的新分解和候选池审查。`figures/` 提供出版尺寸PNG/PDF。所有GT评分与在线求解分离。

主结果不是新算法已经普遍胜出。阈值预测了三例首次大幅box切换；整数事件在C的AP停滞状态上改善了材料，但低于原始NN入口直接运行已有DR的结果；B持续失败。计时小计241.448473秒，624次物理LU/6240complex RHS；完整端到端秒数和峰值内存未测，不填零。

## 复现环境

本轮环境和版本在 `provenance/research_manifest.json`。数值实验使用单BLAS线程，Python、NumPy、SciPy；离线分析另外使用pandas/matplotlib/threadpoolctl。

历史算子、初值和原始投影实现来自提交：

`ec3a37da6f3ee90beb61e753067d4f639f4118f5`

请使用该提交的独立checkout，或读取等价文件，不要覆盖历史R4–R12。原论文全文和整个历史仓库未重复打包。

```bash
export MAX_STRUCT_REPO=/absolute/path/to/pinned/repository
export MAX_STRUCT_RUN_ROOT=/absolute/path/to/this/unpacked/package
python work/normal_memory_probe.py --case A__r10_twod_0004 --steps 3000 --follow 500
python work/normal_memory_followup.py --case A__r10_twod_0004 --steps 3000
python work/normal_memory_followup.py --case A__r10_twod_0004 --steps 3000 --integer-only
```

另两个case_id为 `B__r12_confirm24_0010` 和 `C__r8_test_0001`。这些命令会写同名结果；为保护本包的冻结结果，请先复制本包或将MAX_STRUCT_RUN_ROOT指向新目录，并创建work子目录。现有脚本是本轮机制脚本，不是生产求解器。

不重跑求解器，只做离线评分/作图：

```bash
python work/plot_normal_memory.py
python work/audit_saved_evidence.py
python work/plot_evidence_audit.py
```

后两项需要原仓库的保存材料、metadata和RECEIVER_DPCS。候选池审查还需要公开release中的10个CANDIDATES数组：

```bash
python work/fetch_paired_candidates.py
python work/analyze_candidate_pool.py
```

下载脚本仅使用公开GitHub release普通Range请求，约23.54MB；如果服务器不支持206 Range，脚本会说明失败，不会偷偷下载整个归档。原始检索与CRC回执已保留。不会把重新评分当成独立算法确认。

`provenance/original_scripts/` 是实际运行代码原件；`work/*.py` 是只修改路径配置的便携副本。路径编辑前后的hash在 `provenance/path_edits.json`。未修改任何算法公式或冻结数值输出。

## 后续实现

`CODEX_TASK.md` 是直接可交付任务。其接口是待实现设计，不是声称本包已有完整新算法、成熟NIE/CSI/WRI复现或三维实测结果。先复用仓库已经验证的Green/FFT/adjoint实现，再补因果对照和新的独立确认。
