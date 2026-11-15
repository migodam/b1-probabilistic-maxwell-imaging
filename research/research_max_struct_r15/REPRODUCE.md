# 复现：源码、实际入口与恢复范围

本页给出实际使用的程序，不把尚未执行的命令写成历史运行。算法源冻结提交ba6d68a87；所有R15新增文件留在独立目录，R4–R14不覆盖。完整原始文件从发布后DOWNLOADS.json列出的实际release取得，恢复说明见[DOWNLOAD_GUIDE](DOWNLOAD_GUIDE.md)。

## 1. 环境与路径

主264反演和fork56实际在Windows i9执行：Python3.10.9、NumPy2.2.6、SciPy1.15.3、psutil7.2.2；单BLAS线程。12lane使用逻辑CPU0,2,…,22，对应已核实8P＋4E，同对象所有臂在同一核心顺序运行。没有CUDA计算或NN训练。

实际研究根为`D:/AI/MAX_STRUCT_R15_TWOFOLD_DR_20261116`，实际继承库为`D:/AI/MAX_STRUCT_R11_DEV_20261107`，Python为`D:/python/python.exe`。源文件中仍保留原路径元数据，不假称每台机器具有这些目录。公开的19个实际继承源快照在`sources/legacy_runtime_snapshot/legacy_sources/`；也可将这一真实目录作为`--legacy-repo`，而不是重新发明模型。

离线分析/最后绘图在Apple Silicon Python3.13.13、NumPy2.5.3、SciPy1.18.1、Matplotlib3.11.1完成。不同BLAS/硬件会改变同CPU下完成步数，120秒结果不保证逐位相同；同配置固定步重放与物理/材料指标应分别核查。队列需要CPU affinity，不能把不支持该接口的Mac直接伪装成同一固定核环境。

## 2. 真实反演入口

下面是本轮原远端的命令形式；冻结目录已有输出，程序会拒绝覆盖。重新运行应先把源码、configs和public输入恢复到新的空实验工作副本，替换`--root`为该实际新目录，保留manifest内容和任务ID。

```powershell
D:/python/python.exe D:/AI/MAX_STRUCT_R15_TWOFOLD_DR_20261116/src/queue_independent.py --root D:/AI/MAX_STRUCT_R15_TWOFOLD_DR_20261116 --legacy-repo D:/AI/MAX_STRUCT_R11_DEV_20261107 --manifest D:/AI/MAX_STRUCT_R15_TWOFOLD_DR_20261116/configs/twofold_exposed24_v1.json --cores 0,2,4,6,8,10,12,14,16,18,20,22
```

全队列264任务。单任务真实dispatcher为：

```powershell
D:/python/python.exe D:/AI/MAX_STRUCT_R15_TWOFOLD_DR_20261116/src/run_independent.py --root D:/AI/MAX_STRUCT_R15_TWOFOLD_DR_20261116 --legacy-repo D:/AI/MAX_STRUCT_R11_DEV_20261107 --manifest D:/AI/MAX_STRUCT_R15_TWOFOLD_DR_20261116/configs/twofold_exposed24_v1.json --task-id r14_new_0021__twofold_hop --cpu 0
```

`run_independent.py`同进程分派八个stateful臂到`run_stateful.py`，CIE/Dual/DR到原`equal_budget_C2_runner.py`，不重置CPU时钟。原字段优先级见[CONFIG_PRECEDENCE](CONFIG_PRECEDENCE.md)。不能把mf、β、LBFGS8、period16或候选规则改掉后仍称本轮冻结复现。

## 3. 同状态干预与材料反事实

主队列全部结束后，`queue_interventions.py --root <实际R15根> --legacy-repo <实际继承根>`建立固定56任务：四对象checkpoint512的五个64步base，以及0021/0023、三方法、三个实band、正负16步扰动。实际V2选择均为512，没有fallback。该队列默认只在四个已声明核心运行，不与主实验抢资源。

底层`intervene_fullstate.py`接受真实`--source-checkpoint`、`--public`、`--geometry`、`--case-id`、`--method`、`--band`、`--sign`和可选`--cpu`；必须使用同一材料/current/Q/PR/stage/记忆。只给相同m不等于本轮同状态控制。`compare_readouts.py`在相同四个checkpoint运行G19/S19/U11/邻域各一次，输出实际梯度、方向、原forward变化；该16任务在新GT访问前完成。

完整checkpoint含typed状态与分离的固定算子/背景投影缓存。Windows生成的40条fork缓存相对路径含反斜杠。Windows直接重放；POSIX读取器只在解析元数据引用时把反斜杠转换为斜杠，不重写原NPZ。全部68个非null引用在关键选集已齐，见原`LOCAL_PORTABLE_DEPENDENCY_RECEIPT.json`。窄`STATE_016.npz`仅供材料/current/field比较，不含完整PR历史，不能用于声称完整状态JVP。

## 4. 封存、评分和分析

在线入口不接受GT。主manifest、terminal queue、全部候选及结果先封存；原批准记录为`predictions/PREDICTION_SEAL_MAIN_APPROVED_V1.json`。56fork及16读出另有`PREDICTION_SEAL_FORKS_AND_DIAGNOSTICS_V1.json`。首次评分只因缺实际helper导入失败，尚未访问新GT；修正导入路径后评分，全部费用留存。

实际主分析入口为`src/analyze_campaign.py --root <实际R15根> --seal <该根中的批准seal>`；`--make-seal <新路径>`只生成未授权seal，不自动读GT。`analyze_selection_geometry.py`复用同一批准seal与原端点参考。`score_fixed_forks.py`只离线评分固定4对象反事实。`analyze_fullstate.py`不读GT，需全部原960 cycle before/after状态，不是只需快捷Pro选集。

**恢复后的旧seal不能直接冒充新位置的seal。** 原seal记录的是当时精简审阅树的绝对路径、字节数与JSON快照；随后为Pro恢复了额外完整checkpoint/cache。原seal与当时分析回执保留不改。若需重算全部评分，先在新的空副本恢复全raw，生成针对该实际树的新派生seal，记录“旧预测的离线重评分”，批准后再运行分析；原授权JSON不可简单改路径后称原封存。已有analysis目录也不覆盖。脚本的路径/字节规则不是密码学完整性证明，未使用SHA。

GT源为既有R14的24个`sealed_truth/independent24`材料，另在离线Pro包明确标为OFFLINE_SCORING_ONLY；它们不进入`data_public`。恢复full raw不包含旧端点GT衍生诊断；这些已明确放在离线分析包。不要将旧GT参考当独立采集的连续真值证书。

## 5. 图表和网页

`src/render_final_figures.py`从已封存结果的原指标画最终图，没有重评分或修改成功标准。`FIGURES/final/RENDER_ATTEMPTS.json`保存两次实际版式修正。原始完整状态图在`FIGURES/dynamics/original_remote_*`保留；修正图只修图例/零分母显示，不改数组。`plot_spectral_selection.py`复用现有frame，不重算SVD。

`src/build_site.py --root <实际R15根>`将实际Markdown、CSV和图呈现为静态网页，不运行forward、不读GT、不发布。release的匿名下载回读使用R14已存在的`src/verify_public_release.py --root <实际R15根> --tag research-max-struct-r15-twofold-dr-20261116`。网页访问和附件大小/ZIP签名核验是发布检查，不是重做科学验证。

## 6. 准确的复现承诺

公开源码和无GT观测可重新运行算法；完整raw保留全部current、field、PR、候选、阶段、初始化和费用。离线分析及关键Pro选集提供可审查快捷入口，未声称它包含每个原始文件。参数、成功门和在线规则全部冻结，无新增NN、三维或实测实验等待执行。
