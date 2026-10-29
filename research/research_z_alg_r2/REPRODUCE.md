# 复现与直接求解

[下载完整代码、观测、真值、收据与图像包](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/research_z_alg_r2/FULL_REPRODUCIBILITY.zip)。解压后进入`research_z_alg_r2/`；所有重新运行的结果应写入新目录。包是科学冻结快照，线上发布回读和最终交付时间单独见仓库根目录下本研究的`PUBLIC_READBACK.json`和`CLOSURE.json`。

## 最短入口：重新求一个已归档观测的全部候选

在 `research/research_z_alg_r2/` 中：

```bash
python -m pip install -r requirements.txt
python run_solver.py --geometry experiments/QUOTIENT_HOLDOUT_CONTRACT.json --observations experiments/quotient_holdout_public/QH4_03.npz --box 0.5 30 --method quotient --output /tmp/zalg_qh4_reconstruction.json
```

该真实入口已运行，返回2个完整验证合格材料，支付10个完整复RHS，原收据为 `results/CLI_REPRODUCTION.json`。若输出路径已存在，换一个新路径；程序拒绝覆盖证据。

本次机器使用：

```bash
/Volumes/migodam\'s-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python run_solver.py --geometry experiments/QUOTIENT_HOLDOUT_CONTRACT.json --observations experiments/quotient_holdout_public/QH4_03.npz --box 0.5 30 --method quotient --output /tmp/zalg_qh4_reconstruction.json
```

建议用自己环境的`python`，无需复现绝对路径。已用版本Python3.13.13、NumPy2.5.3、SciPy1.18.1；矩阵与图形采用CPU，脚本固定BLAS单线程。没有CUDA、服务或训练依赖。

## 48单元观测成像

```bash
python run_solver.py --geometry experiments/HIGHFREQ48_CONTRACT.json --observations experiments/highfreq48_public/HF48_00_relative_complex_1e-3.npz --box 0.5 4 --method hybrid --relative-noise 0.001 --output /tmp/zalg_hf48_reconstruction.json
```

该命令是可运行通用入口；归档论文计时来自原冻结`experiments/run_highfreq48.py`及对应source snapshot。通用CLI报告噪声尺度资格；归档实验原1e-8精确fit标志不被CLI更改。`--method svd`可重现未加校正的噪声失效，`--method fullj`使用固定中点。

## 四方法完整新保留集的原始执行命令

原执行：

```bash
python experiments/generate_quotient_holdout.py
python experiments/run_quotient_holdout.py --freeze
python experiments/run_quotient_holdout.py
python experiments/score_quotient_holdout.py
```

生成/执行脚本为保全原证据会拒绝在已完成目录重跑。实际一次driver guard预物理修复及`--freeze-repair`收据保留；它未改变核心求解器或参数。当前发布中的`private/`是**运行期间隔离、预测封存后才供评分/复现使用的真值**，现在随完整证据公开；重新读取这些值的重放应称exposed replay，不能重新称fresh blind。

重新计算相同四方法比较，使用独立输出的重放入口，无需修改原freeze或删除原结果：

```bash
python replay_holdout.py --output-dir /tmp/zalg_replay_new
```

该入口只读取公开观测和冻结的初始化，调用原四种求解函数，输出独立的逐例候选、物理费用与集合对比。它不读取GT；发布后的重放仍应标识为`EXPOSED_REPLAY`。源文件size/mtime freeze是原进程隔离的审计记录，在他机复制后会变化；不要伪造原freeze来声称新盲实验。要产生新的真正盲研究，应重新冻结合同并生成未暴露数据。

## 其它实验入口与对应结果

|入口|结果/合同|用途|
|---|---|---|
|experiments/run_branch_study.py|BRANCH_CONTRACT，branch_summary|24路径/六路径/FullJ8首次比较|
|experiments/run_adaptive_branch.py|ADAPTIVE_BRANCH_CONTRACT，adaptive_branch_summary|共享Schur核与可行切片|
|experiments/run_high4.py|HIGH4_CONTRACT，high4_*|独立新高对比度支持|
|experiments/run_quotient_benchmark.py|QUOTIENT_BENCHMARK_CONTRACT，quotient_bench_summary|101已暴露覆盖及5次6任务计时|
|experiments/run_quotient_control64.py|QUOTIENT_CONTROL64_CONTRACT|两例更强64起点控制|
|experiments/run_highfreq48.py|HIGHFREQ48_CONTRACT，highfreq48_*|全部48材料盲成像|

全部脚本的首次运行保护保留；不能在已保存结果上重复执行并覆盖原成本。部分历史初始化脚本引用同仓旧E1几何，核心`run_solver.py`、QH4与HF48输入是自包含的。

## 零物理重画与审计

```bash
python experiments/plot_final_results.py
python experiments/plot_highfreq48.py
python experiments/plot_blind.py
```

这些脚本只读已保存数组。主图/附图提供PNG、SVG、PDF。理论公式不依赖大量随机恒等式测试；严格推导、数值实施和完整物理验收分开。

## 费用与数据完整性

`COST_LEDGER.csv/.json/.md`列出原阶段345份规范收据（追加阶段后349份）及未保存完整原始端点的28开发RHS。历史anchor重复用于公平费用归属，但新增支付只计一次。数据生成也计费。没有新增SHA-256检查；本地路径size/mtime冻结记录不是密码学承诺或操作系统隔离。

## 追加六材料盲确认的可复现入口

```bash
python replay_general6.py --output-dir /tmp/zalg_general6_replay_new
```

它读取`experiments/general6_public/`和`GENERAL6_CONTRACT.json`，调用冻结的20维商/ReducedFullJ64求解方法。原始研究生成前冻结solver，预测封存后才评分，3份观测返回2/1/1对1/1/1根；原始完整成本48 RHS。重放命令已语法检查但未另行运行，避免对已完成科学结果重复付费。原始执行代码与完整收据均已实际运行。

重画新增材料根图（仅评分后可使用GT显示，不运行物理）：

```bash
python experiments/plot_general6.py
```

总账更新为53,951完整复RHS、349条规范收据；新增48个RHS的逐项分解保留在账本中。所有主文墙钟为同一Python进程内算法计时，包含物理准备与目标输出；公共模块导入及启动Python进程不计入各方法时间，也未另行测量，不能把毫秒表格理解成独立冷进程启动延迟。
