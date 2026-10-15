# Research B parent-level统计实现与已执行回执

脚本：`research/research_b/code/analyze_results.py`。已在原生 `Gaussian/.venv_nn/bin/python` 环境执行 `--stages E0 E1 E2 --job analysis-e0-e1-e2-parent-paired-v1`；这次调用在收到主线程通知E2完整后进行，不读取正在增长的E2作为完成结果。E3未混入。

## 统计定义

- E1/E2各32独立parent、每parent两个noise，先平均noise，再对32个parent做2000次配对bootstrap，固定seed26101569。CI为paired parent mean difference的percentile 2.5%/97.5%。
- 差值一律`policy−baseline`；Selected-MAE负数有利。多方法区间是描述性结果，没有多重比较校正或自动优越性裁决。
- 每stage/policy/mode输出initial、step1/2/3（限stage预算）与final。明确stop/none的最后合法行可向后保持；没有stop记录的残缺轨迹被拒绝，不能当提前停止。
- E0七个已暴露development cases仅作描述性汇总，不伪装成32parent，不给正式bootstrap CI。
- mode1是candidate profile选材；EIG、task risk、calibrate-then-risk、risk stop、bundle的选action仍标为reference-assisted，不能称整个pipeline candidate-only。
- ESS和finite-reference risk描述共享reference particles，不是candidate posterior；mode1 candidate regret是独立offline reference诊断。
- 源数据的“selected structure correct”实际使用12个事后central-radius/annular-width事件bin。它与posterior最大质量event决策分别输出，不外推为任意拓扑恢复。
- mode1保持同625固定候选池；mode2的刷新只是从现有2048参考支持中增加输出catalogue。Best-of-K变化不能称全材料新解发现。
- E1 fixed3.2是事后强fixed描述性baseline；E2 fixed3.2为E1之后冻结、同parent的暴露跟进。none比较同时披露测量block数差异。

## 产物

`results/SUMMARY_METRICS.csv`：400个stage/mode/policy/horizon组。

`results/PARENT_METRICS.csv`：11150个逐parent摘要行；E1/E2已经各自平均两noise。

`results/PAIRED_COMPARISONS.csv`：18324个配对比较/不可用声明。baseline为fixed2.4、fixed3.2、none；缺失指标保留null，不选择性删除parent。

`results/STATISTICAL_SUMMARY.json`：上述统计、availability、12行补偿机制表、forwardaudit与初始reference数值审计。这里不写最终科学结论。

## 关键已执行数值

以下均E2 final，delta为相对fixed3.2的Selected-MAE差，32parent bootstrap：

|选材模式|policy|delta|95% CI|
|---|---|---:|---|
|candidate profile|robust|+.002956737|[−.002879278,.010340840]|
|candidate profile|expected material risk|+.002889604|[−.000581489,.007753838]|
|candidate profile|event EIG|−.000778148|[−.005216503,.003527961]|
|reference selection|robust|+.004934768|[.000634832,.010134106]|
|reference selection|expected material risk|+.000291159|[−.001563784,.002216171]|
|reference selection|event EIG|+.002883724|[.000837980,.005081399]|

E2 fixed3.2的candidate/reference Selected-MAE为.090884949/.039812287；none为.162278670/.104045517。候选固定池Best-of-K=.062224210，Coverage@K=.40625始终不变。

## 成本范围

`COMPUTE_LEDGER.jsonl`中独立ANALYSIS receipt已append：wall2.243842416秒、process CPU2.19761秒、new Maxwell0、GPU0。`Cost.job`计费边界包含统计读取、bootstrap与四文件导出；不把它说成全部解释器启动、对话与交互式写代码费用。未调用`Cost.evaluate`、solver、训练、MC物理或A目录写入，也没有新增gate或smoke任务。
