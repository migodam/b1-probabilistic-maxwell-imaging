# 复现

需要Python3.11+、NumPy/SciPy/Matplotlib，无需pytest、CUDA、后台服务。项目实际使用已有A25环境；新环境可安装requirements.txt。所有命令应从本目录运行。设置 `PYTHONPATH=src`；示例的 `python` 指你的科学Python。

运行前阅读w_config.json和PHYSICS_PROTOCOL.md。公开结果保留原样。缓存重画需要从[本轮Release](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-w-p0-20261022)下载Research_W_P0_ARRAYS.zip，按内部相对路径解压至本项目根目录。较大的墙体NPZ和逐网格field CSV放Release，其他摘要和四块小数组在Git。

若要独立从头重算，用 `python tools/new_run.py /absolute/new/empty/directory` 生成不含旧结果的新实验目录；它只复制协议/源码/只读输入，不启动求解，不是重置原campaign的账本。之后进入新目录执行下面命令。原公开项目及其负结果保持不变。

```sh
python tools/run_job.py local-analytic-001 python -m research_w.analytic --output results/p0a --source-results source_package/analytic_results/results.json
python tools/run_job.py local-tests-001 python tools/run_tests.py
python tools/run_job.py local-benchmark-001 python -m research_w.qualification benchmark
python tools/run_job.py local-cylinder-001 python -m research_w.qualification cylinder
python tools/run_job.py local-derivatives-001 python -m research_w.qualification derivatives
python tools/run_job.py local-wall-grids-001 python -m research_w.qualification walls
python tools/run_job.py local-parity-001 python -m research_w.qualification assembly-parity
python tools/run_job.py local-reciprocity-001 python -m research_w.qualification reciprocity
python tools/run_job.py local-wall-average-001 python -m research_w.qualification walls-averaged
python tools/run_job.py local-cylinder-average-001 python -m research_w.qualification cylinder-averaged
python tools/run_job.py local-derivatives-average-001 python -m research_w.qualification derivatives-averaged
python tools/run_job.py local-fourblocks-001 python -m research_w.four_blocks
python tools/run_job.py local-figures-001 python tools/plot_results.py
python tools/run_job.py local-summary-001 python tools/summarize.py
python tools/run_job.py local-report-001 python tools/report.py
```

仅重画已公开数组和指标不需要Maxwell：运行plot_results.py、summarize.py、report.py即可。公开 `results/jobs/*/RECEIPT.json` 保存真实执行命令及其scope，包含失败。

P0-A的独立入口：`python -m research_w.analytic --help`；上述命令独立保存在新results/p0a并比较原归档。原ZIP脚本若直接运行会写到它的旁边；只允许在自己的可写副本运行，不能写只读source_package。

本轮没有实现/执行完整四种见证搜索、未知目标joint正确器、P1或3D。它们不是隐藏的未公开成功代码。资格未通过时不能绕过gate去调用未来搜索器。
