# 复现与证据读取

## 读取完成结果

先读RESULTS_AND_METRICS.json和evidence_summary。完整Release包含模型、联合训练数组、原始posterior生成样本、物理端点/current/field、J/谱、chain拒绝重复、source snapshots和追加式账本。release分卷按EVIDENCE_PACKAGE说明恢复；不做SHA检查。记录文件size/mtime_ns用于溯源，不能把它们当内容密码学证明。

代码依赖Python3、NumPy、SciPy、PyTorch、matplotlib以及原A30后端；实际版本与CPU/MPS资源见Stage0源清单。原运行解释器为../../.venv_nn/bin/python；所有数值库单线程。已有本地MPS训练，未使用CUDA、云或新3D。

## 可执行入口

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 PYTHONPATH=src
python -m p3r4.train_density --help
python -m p3r4.train_anchor_density --help
python -m p3r4.train_geometry_cfm --help
python -m p3r4.train_risk_fusion --help
python tools/prepare_r4_fresh_risk.py --help
python tools/infer_r4_fresh_density.py --help
python tools/run_r4_interval_chain_locked.py --help
```

精确as-run argv、环境、source snapshot和预算收据在各stage/jobs/<job>/START.json、COMMAND.json及CLOSED.json；模型RESULTS/TRAINING日志包含训练seed、容量、选定update及normalizer。不得清空原Book重新运行；复现实验应建立新的独立输出及授权，或只读取已完成数组。

原封存manifest含当时绝对路径。下载证据后先按packaging清单把原根映射至新根，保留原封存文件，使用新副本/明确路径映射执行replay；不改原snapshot身份。源文件的新运行与旧mtime身份核验应分开，不用touch掩盖版本差异。复现关闭账本的正确含义是读取证据和生成相同算法的新结果，不是把新机器伪装成原进程。

输入身份：原joint x/q是TRAIN监督标签；anchor_conditions/context_x/context_J才是观测驱动context，二者严禁互换。新heldout GT只在offline_only目录，在线过程不导入其标签。模型元数据K修复见AUDIT/fusion/providers/CLOSED_K_METADATA_CLARIFICATION.json，只读compat保留原checkpoint字节。

## 报告再生

tools/build_r4_final_reports.py从closed评分和账本产生紧凑指标及中文报告，不重新求解、训练或选择模型。tools/plot_r4_final_risk.py及最终各figures目录的code/INDEX提供画图来源；图片均来自实际数组，不用生成图代替实验。

上述builder对应原封存报告草稿，不会自动生成后来追加的出版地址、训练解释或标准校准下界。最终版本默认阻止覆盖这些补充；若要还原草稿，请在独立checkout明确使用`--overwrite-final-documents`。只重算已公开材料风险可使用[可移植离线重放](PUBLIC_REPLAY.md)，不会运行训练、物理求解或改写原结果。静态网页由tools/render_r4_final_site.py读取最终Markdown生成。

原P1/A30/R2依赖与实际Git revision列入SOURCE_INVENTORY和证据包。发布不额外运行物理smoke campaign，读回网页或Release不当作科学gate。
