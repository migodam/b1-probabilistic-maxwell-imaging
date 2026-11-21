# 四个实验怎样接入，怎样用封存数组验收

在已有工作区，从项目根运行下列一条命令就能验收四个场景，并重算二维主半环24例。此次读取了原始观测、已选材料、必要的已曝光真值和历史对照表，新增反演和 Maxwell 场求解均为零。最能说明复算价值的是二维半环。30秒共同校正后，SRAP 在24例中的22例降低总材料误差，但内部误差中位数从0.621升到0.799。背景误差下降解释了总误差的改善。[逐物体复算 CSV](results/half24_pairs.csv)

```bash
Gaussian/.venv_nn/bin/python -B research/experiment_handbook/replay_saved.py all
```

公开仓库从仓库根运行同一脚本，将前面的解释器换成有 NumPy、SciPy 的 `python`。可用子命令为 `r11`、`austria`、`fresnel`、`vector`、`half24`。单独运行 `half24` 只需 NumPy、SRAP 的二维合成数组和历史表。新进程验收五项共 **1.260秒墙钟、0.564秒 CPU**，其中半环指标和历史表核对用0.099秒墙钟。进程费含导入、读取、计算和输出。历史反演费用见各节。[复算程序](replay_saved.py) · [完整进程计时](results/replay_process_cost.json)

### 数据放在哪里

|下文别名|本工作区从根开始的路径|公开仓库从根开始的路径或还原来源|
|---|---|---|
|`R12`|`research/delegated/max_struct_r1/r12_worktree/research/research_max_struct_r12`|`research/research_max_struct_r12`；补充实验位于 `supplement_paired/`|
|`R10`|`research/delegated/max_struct_r1/r10_worktree/research/research_max_struct_r10`|`research/research_max_struct_r10`；可额外比对真实源预测|
|`R8`|`research/delegated/max_struct_r1/r8_worktree/research/research_max_struct_r8`|`research/research_max_struct_r8`|
|`SRAP`|`research/srap_cie_array_study_20261121`|相同路径；数组附件解压到本研究目录|
|`FRESNEL_2005_DIR`|`research/delegated/a2_literature/data/2005_iop_21_6_S09`|取原始七列 `FoamTwinDielTM.exp` 所在目录|

复算器自动在公开目录和本地工作树之间查找。数组放在别处时，设置 `R8_DIR`、`R10_DIR`、`R12_DIR`、`SRAP_DIR`、`FRESNEL_2005_DIR`。需要读取兼容源码的 Fresnel 校准，`R8_DIR` 上两级应是检出仓库根，内含同版本 `research/` 依赖。默认写到当前目录的 `.artifacts/experiment_handbook_replay/`。可用 `--output` 指定本轮输出，保留本手册的冻结结果。以下各节的“历史入口”用于了解原求解合同，本次验收实际执行的是上面的只读命令。

## 1. 同一个真实坏终点，换哪一种后续算法

**入口是已保存的材料和同源测量。** 取 `R12/supplement_paired/inputs/manifest.json` 的 A 组三行，从 `inputs/A__r10_twod_*.npz` 读 `m_init`、`E0`、`G`、`S`、`Y`、材料上下界及 CN99 噪声尺度。三例各有24×24个实材料未知量，12源、48接收，矩阵形状依次为 `E0(576,12)`、`G(576,576)`、`S(48,576)`、`Y(48,12)`。历史工作树的绝对路径可以按 `parent_task_id` 还原为 `R10/RESULTS/twod/<parent_task_id>/prediction.npz`。

|R11中检查过的对象|原预测的真实方法|四个后续臂的实际起点|
|---|---|---|
|`r10_twod_0001` 环|`cie_single`|`prediction.npz/m`|
|`r10_twod_0004` 开缝壳|`matrix_green19`|同上|
|`r10_twod_0022` 桥|`matrix_green19`|同上；该终点缺乏严格局部极小资格|

**验收。** 本机比对3个真实父预测、3个公共输入和12个方法起点，材料数组逐元素完全相等。对12个已选 `primary.npz` 独立检查 `Q = E0 + G @ j` 和 `Y_prediction = S @ j`，最大相对差分别为 1.73×10⁻¹⁵以下和2.19×10⁻¹⁵以下；重算的复数据残差与该候选的 `primary.json` 一致。此验收不读材料真值。[12条验收记录](results/r11_replay.json)

**历史求解与评分。** `prepare_inputs.py` 形成公共输入，`src/run_paired.py --root <新实验副本> --legacy-repo <兼容工作树根> --task-id A__r10_twod_0001__dr --cpu 0` 是历史单任务入口。G19 重启会把父材料的物理电流投回原 DPCS/Fourier 空间、重置 PR；AP、DR 保持同一起始物理场，direct 从同一材料优化原物理目标。各法公用 box [0.02,20] 和同一 smooth/TV；A组限制为每臂新增240个整进程 CPU 秒，带候选的三法分配144秒候选、96秒共同校正，direct 可用全额。在线主选择取截止前完整候选中原物理 Phi 最小者。全补充共44个实际子任务，使用 `ALL44_PREDICTIONS_SEALED.json` 封存后才运行 `src/score_and_plot.py`；比较成本读取 `RESULTS/offline_v2`。这三例属于已曝光困难起点，历史材料评分中 G19 重启赢环，DR 赢开缝壳和桥。换算法时应带上G19重启这个强对照。

**还原。** Release 的 `runtime_inputs.zip` 含公共输入；按附带 `ASSET_MANIFEST.json` 顺序拼接 `paired_results_full.zip.part001/.part002`，再解压。工作区有父预测时程序会验证原始数组；缺少另行归档的 R10 父预测时，验收记录会标 `original_parent_array_checked=false`，仍可核对四臂相同起点和实际物理场。[R12补充报告](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_max_struct_r12/supplement_paired/) · [运行包和完整结果](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-max-struct-r12-paired-20261112)

## 2. Austria 的150²观测能否直接用于32²和64²反演

能。两组历史求解都从 `R8/RESULTS/paper_benchmarks/austria_eps4_400MHz/PUBLIC_DATA.npz` 读取**同一个** `Y(40,20)`。本例固定400 MHz、物体介电常数4、2 m正方域、20个平面波和半径3.75 m的40个接收位置。原始 `DATA_RECEIPT.json` 记录150²细网正演和10%复噪声。32²是计算适配，64²使用更细的反演网格。比较要保留这两种网格分别对应的材料真值和求解费用。

|已封存的逐源双β CIE|材料向量|输出复场|本次从同一Y重算数据残差|历史该臂CPU秒|
|---|---:|---:|---:|---:|
|32²，`<case>/cie_dual/PREDICTION.npz`|1024|40×20|0.0980128943|17.220|
|64²，`original64/<case>/cie_dual/PREDICTION.npz`|4096|40×20|0.1005967768|318.121|

**验收和历史入口。** 表中 `<case>` 是 `austria_eps4_400MHz`，两个相对路径均起于 `R8/RESULTS/paper_benchmarks/`。直接读取 `m`、`Y_prediction` 和封存回执，数据残差为 `norm(Y_prediction − Y)/norm(Y)`，与两个 `RECEIPT.json` 完全相符。32²原程序位于 `experiments/paper_benchmark.py` 的 `solve`；64²为 `experiments/paper_compare64.py --arm cie_dual`，使用 `acquisition(64,...,dense=False)` 和 `MatrixFreeModel`。复用观测时可以跳过 `paper_benchmark.py build` 的细网求解。历史32²评分入口为 `paper_benchmark.py score` 与 `paper_report.py`，64²离线评分入口为 `paper_followup_report.py`；这些程序会读私有材料并写旧账本，应在新副本中使用。公开历史评分中的相对材料 L2 分别为0.8257和0.8143；两行沿用各自物理迭代与成本合同，以同一Y衡量的残差和材料排序可以不同。[本次数据核验](results/austria_replay.json) · [R8实验说明](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_max_struct_r8/) · [R8归档的 paper_benchmarks 附件](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-max-struct-r8-20261106)

## 3. Fresnel 2005 Twin，先核对复数约定，再划训练与留出

**从原始七列文件开始。** `R8/src/fresnel_adapter.py::load_file(path,8.)` 在 `FoamTwinDielTM.exp` 中选出8 GHz的4338行。18个源，每源241个实际接收位置，用 `(360,18)` 矩阵和缺失值容纳；原始编号从1开始，NumPy 索引从0开始。列序为源号、接收号、GHz、带物体总场实部/虚部、空场实部/虚部。项目对原始差场和空场一起取共轭，接到 outgoing Hankel-1 的 exp(−iωt) 实现。本次直接对4338个复数逐元素核对差场与空场转换。

**空场校准只用移走物体的数据。** 发射与接收半径均1.67 m；每源在与发射位置相对的接收方向±30°内，用五个柱波模式 m = −2,…,2 拟合空场。按角序每三个候选留一个验证，其余用于拟合；保留奇异值超过最大值 10⁻⁴ 的模态。本次重新调用 `incident_fit`，18源空场验证误差中位数0.114904，生成的训练入射 `E(3600,9)` 与保存数组的相对差为0。这个校准仍包含接收喇叭方向图带来的模型误差。

**源依赖通道。** `SRAP/data2d.py::fresnel_selection` 和 `run_fresnel` 生成两个布局，各用9个源和733个有效复通道。保存的 `S(360,3600)` 保留全部接收槽，`Y(121,9)` 含填充；第 s 列只取 `valid_mask[:,s]` 为真的位置，并通过 `rx_ids[s,:]` 查原始接收器。连续训练源（零基）为0至8，交错训练源为0,2,…,16。两布局共用5个完整未见源9,11,13,15,17的1205个通道，并在5个共同已见源0,2,4,6,8上另留出390个接收通道。复算逐源核对训练值、两个留出数组及不相交关系。原始散射值、空场拟合误差、未见源误差分别有自己的分母。实物没有可逐像素打分的材料真值。[校准与遮罩验收](results/fresnel_replay.json)

**历史求解与评分。** 在60²、0.2 m域中，`SRAP/run_study.py --cohort measured_extra --arms dual,srap` 载入经过校准的公开包。按源归一化的训练场目标选择15/30秒端点，共同校正从80%时间前的公开种子开始。完整预测后 `score_study.py --cohort measured_extra` 重算未见源和另一接收侧的相对复场残差。历史30秒校正时，Twin 连续阵未见源残差 Dual/SRAP 为1.0564/1.0593，交错阵为0.9128/0.9357。这是原已评分表，当前回放验证了原始相位、校准与评分分割。

公开 SRAP 数组归档不含可替代原作者七列文件的完整测量输入包。有原件后，以下命令只生成 Twin 的公共包与留出包到新目录，不会生成合成训练集。按需设置 `FRESNEL_INPUTS` 让本手册复算器读取该目录。原件与代码的使用范围见原始数据来源。[原始实验与数据附件](https://doi.org/10.1088/0266-5611/21/6/S09) · [发布实验的复现说明](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/main/research/srap_cie_array_study_20261121/REPRODUCE.md)

```bash
export FRESNEL_INPUTS="$(mktemp -d)"
python -B research/srap_cie_array_study_20261121/data2d.py --cohort fresnel --output "$FRESNEL_INPUTS"
python -B research/experiment_handbook/replay_saved.py fresnel
```

## 4. 半球阵用的是三维矢量场，合成与反演各用哪张网格

在 `SRAP/vector/` 下，`generate.py` 以24³细网、每体元64个子采样点生成材料，`physics.py::Model` 构造全三维自由空间 dyadic Green、体元体积和辐射修正自项。8个不同发射点各有两支正交偶极，16个接收点各取两个正交分量，观测 `data/hemi_00.npz/Y` 为 `(32,16)`，含512个复数。12³反演网格有1728个实材料未知量，三个电流分量及全部发射共用每体素一个 χ。背景波长为单位，域宽0.6λ、半球半径0.9λ，观测为细网场加10%复噪声。正式另一布局是二维矩形面阵；早期共面圆环留在 `development_ring/`。

**历史求解与评分。** `invert.py::optimize("hemi",0,"Dual")` 或 `"SRAP"` 从同一 χ = 0.25 起步。两臂同用逐源双β及原物理目标，取60/180秒前缀；light 从48/144秒付费种子出发，用12/36秒做最多20个接受步。SRAP 的传播路由为5184维稠密谱分解，图投影用 FFT-CG，物理正演及伴随用 FFT-BiCGSTAB。冷启动费用包括每臂设置、路由和审计；两个预算前缀来自同一轨迹。三个分量的电流不是三种各自独立的材料。

**验收。** `score.py` 检查48臂完整后才打开 `sealed_truth/`。本次对已曝光 `hemi_00` 验证48个正式臂的封存清单、收发极化对正交性，以及24³真值按2×2×2求均值后与保存12³真值逐元素相等。直接从 `predictions/hemi_00_{Dual,SRAP}/materials.npz:light_180` 复得材料 L2 为 **0.555770/0.535641**，与 `objects.csv` 一致。该壳的非球形细化对照在半球阵从12³→24³的4.158%降至24³→36³的1.712%；独立 Mie 球复远场误差从7.946%降至6.152%。这些数衡量指定离散模型和尺度。[半球单例与网格验收](results/vector_replay.json) · [完整三维报告](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/srap_cie_array_study_20261121/vector_results.html)

## 由材料数组独立复算二维主半环结论

只取 `data/cases.json` 里 `cohort=confirmation` 且 `array=PW180_10x20` 的24个不同物体。固定10个平面波、20个接收，96²细网观测、32²反演，400 MHz、10%噪声。程序先读取全部48个 `results/predictions/<case>__{dual,srap}.npz`，确认 `m_30_raw` 与 `m_30_light` 都完整，再打开这24个已曝光 `data/truth/<case>.npz/m`。这一步从端点材料重新计算，历史表只在复算完成后对照。

误差向量取预测材料减真值，内部定义为真值大于零的体素，背景为真值等于零的体素。总误差、内部误差、背景误差均除以完整真值的 L2 范数；分子分别取全域、内部、背景上的误差 L2 范数。两种恢复均用真值峰值一半的共同阈值计算 IoU。每个物体的收益为 (L₂,Dual − L₂,SRAP) / L₂,Dual，以物体为单位对24个配对收益做1万次 bootstrap，随机种子20261121，取中位数的双侧百分位区间。以下数对均为 Dual / SRAP 的对象中位数。三维报告中的配对收益采用绝对 L2 差，不能直接与本表百分数拼接。

|30秒端点|总 L2|内部 L2|背景 L2|IoU|SRAP 胜/平/负|配对收益中位数及95%区间|
|---|---|---|---|---|---|---|
|raw|1.849 / 1.000|0.554 / 0.754|1.768 / 0.523|0.126 / 0.000|21 / 1 / 2|39.09% [27.18%, 47.35%]|
|light|1.882 / 1.120|0.621 / 0.799|1.779 / 0.838|0.136 / 0.084|22 / 0 / 2|31.90% [23.45%, 37.91%]|

空气初值的相对 L2 为1。raw 优于空气的对象分别为0/24与8/24，light 为1/24与8/24。减少背景伪影、超过对照与正确重建物体内部是三个可以分别量化的结论。独立计算的96行 × 4项（总L2、内部、背景、IoU）与历史 `confirmation_metrics.csv` 的最大绝对差均为0，胜平负和配对区间也与 `paired_summary.json` 一致。新预算或新选择规则不能套用这些历史区间。[逐物体配对](results/half24_pairs.csv) · [数值验收](results/half24_replay.json) · [公开数组附件](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/srap-cie-array-study-20261121)

### 逐物体看半环的正反例

直接看图册第二行的同色标切片。有相同 GT 和预算的两种情况值得一起保留。

- **双环 `s2d_010`。** SRAP raw 的背景已经干净，但两个环的内部几乎不可分，IoU 为 0。轻校正后两处暗孔均可见，L2 / IoU 为 0.472 / 0.619；Dual 轻校正仍有大块背景伪影，为 2.710 / 0.105。这里能看到校正将一个较好的空间位置转为可分辨结构。[完整 7 阵列原图](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/srap_cie_array_study_20261121/figures/confirmation_atlas_30s/s2d_010.png)
- **偏心壳核 `s2d_007`。** 同一半环上，SRAP 形成宽而过亮的环，中心偏心材料没有按真值分开。轻校正 IoU 比 Dual 高（0.342 对 0.185），材料 L2 却是 11.050 对 2.013。同物体的全环轻校正反而给 SRAP 0.412 / 0.694。色标上限为 6；定量误差用原数组。孔径改变了这条路径上的结果。[完整 7 阵列原图](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/srap_cie_array_study_20261121/figures/confirmation_atlas_30s/s2d_007.png)
