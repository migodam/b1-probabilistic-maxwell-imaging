# 复现与读取

[完整复现包](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/gn-green-search-space-r2-20261107/full-reproduction.zip) · [阅读包](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/gn-green-search-space-r2-20261107/reading-pack.zip) · [Release](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/gn-green-search-space-r2-20261107)

完整包保留新研究代码、全部预测/基/候选数组、表格、日志与直接依赖的冻结R4物理、公开观测、原GN端点、已披露GT、上轮Green谱。解压后的顶层为`research/`，可以覆盖到单独的B1复现副本；不要把新运行写进原冻结目录。

## 环境

实际环境：Apple Silicon CPU，Python3.13.13、NumPy2.5.3、SciPy1.18.1；另需matplotlib、threadpoolctl。没有安装新库、没有使用CUDA。每个worker单线程BLAS，至多两个数值worker。运行费用会随机器改变，不能复制本机秒数作为自己的测量。

```bash
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
cd research/research_gn_search_space_r2
```

## 读取冻结结果与重绘（零物理求解）

```bash
python src/summarize_geometry.py
python src/final_summary.py
python src/plot_results.py
```

读取`results/RESULTS_SUMMARY.json`、`results/RESULTS_BY_CASE.csv`、`results/COST_LEDGER.csv`。只重绘不会产生新的成像证据。完整数组需从Release下载，Git只保留轻量结果和图片。

## 搜索入口与独立评分

```bash
python experiments/run_search.py --experiment e2 --workers 2
python experiments/run_search.py --experiment e3 --workers 2
python experiments/run_search.py --experiment e2 --stopping-profile tight --workers 2
python experiments/run_search.py --stationarity
python src/score_search.py seal
python src/score_search.py score
```

入口默认跳过已封存任务；在完整包中再次执行不会被算成新实验。重新运行时，应在独立副本保留旧结果后，将`results/e2`、`results/e3`改名为带日期的归档名；脚本会新建输出目录。E2必须先于E3，以便显式复用12个相同普通对照。正式评分器只在seal后读取`truth_after_freeze`，不会回调优化器。

单对象公共接口：`src.search.solve(public_arrays, anchor, task)`；公共输入可由`load_public(observation_id)`装入。`task`中的科学预算、权重和方法在清单中明确。GT、error_OFFLINE与评分器不进入这个接口。

## 空间与有限步分析

已有结果可以直接读`results/e1/*.npz`和CSV。若确实需要从头复算，在独立副本将整个`results`改名保留，并执行：

```bash
export GN_E1_ONECASE=1
for i in $(seq 1 52); do python experiments/run_geometry.py; done
unset GN_E1_ONECASE
python experiments/run_geometry.py
python experiments/run_geometry.py --fixed-reference-coverage
for i in $(seq 1 52); do python experiments/run_geometry.py --subsets; done
python experiments/finite_and_controls.py
python src/summarize_geometry.py
```

`GN_E1_ONECASE=1`每进程处理一个未完成端点，规避本机曾出现的长worker原生库损坏。最后无此变量的一次调用只完成已有CSV的收尾，不重做52个案例。重新运行的每进程公共reference准备也必须收费。`--subsets`每次只处理一个未完成anchor；12源结果明确复用，1/3/6源重新付费。

历史`--finalize`入口用于本次故障后的恢复审计，含原生失败说明；不要把它当新机器自动产生的失败报告。新运行应以自己的`process_cost_*.json`、`subset_cost_*.json`、E2/E3逐任务receipt为权威计费记录，不能复制本次COMPLETION/失败账本冒充新成本。

## 输出合同与限制

- E1：`cases.csv`、`dimensions.csv`、`current_subset_dimensions.csv`、`coverage.csv`、`drift.csv`和每案例基数组。
- E2：每臂预测、完整receipt、三权重与tight诊断分开；驻点诊断额外计费。
- E3：每臂anchor/top1/raw/filtered/before/after数组、原目标全轨迹、失败、在线排序。
- finite：432个GT定向离线探针与六幅度控制，不能作为公开数据搜索结果。
- 汇总：原始receipt为证据源；主报告已知成本下界不包括无法恢复的原生失败部分。

独立数据对象数是12种低对比材料，不是48端点或576条提案。没有NN权重，也没有隐藏GT种子；本轮没有训练网络。
