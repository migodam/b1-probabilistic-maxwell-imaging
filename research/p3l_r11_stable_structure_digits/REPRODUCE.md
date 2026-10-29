# 复现与完整证据

默认复现只读取已经封存的数组，不支付新的Maxwell RHS。整轮全阶LU预算已失守并STOP，不应在同一账本里重跑物理来补救结果。

## 视觉资格补充的还原

先解压初版完整tar.zst到新的隔离目录，再把R11_report_qualification_update.zip应用到这份隔离副本的报告与分析脚本。它不含任何替换GT或预测。不要覆盖原冻结历史工作树。P3数值fit门仍PASS，但视觉表示FAIL_THIS_FIT；初版数据执行身份和补充报告身份分别固定。

历史R10路径核还需原[p3lr10-bilateral-20261030 Release](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/p3lr10-bilateral-20261030)的shape-update附录输入；本包已包含本轮所有计算输出和数字分析的全部必要输入，不重复整个历史训练数据。

## 环境

实际执行Apple Silicon CPU、float64/complex128、单BLAS线程。现有 `Gaussian/.venv_nn/bin/python`，保留该venv入口，不将符号链接解析为缺包的basePython。依赖为numpy/scipy/Pillow/matplotlib；评分采用明确的7×7 uniform-window SSIM实现，不依赖未安装的skimage。没有CUDA或外部付费计算。

## 证据恢复

Git保存源码、协议、紧凑CSV/JSON、报告和图库。Release的无损tar.zst保留所有原始NPZ、在线y、GT、72份各自的完整预测、迭代标量历史、每次物理调用账本、初始错误及支持输入。按manifest里的仓库相对路径恢复到隔离目录。使用BLAKE2b-256与字节数验证；无SHA-256检查。

inputs/保存原R9精确host setup副本、原P1 prior源码、原配置与物理points，均为无损已有资产复用。历史宿主setup的264 RHS576与1 LU576不是本轮新增费用，也不能因此叫作从零冷启动免费。

## 零新增物理的分析入口

从仓库根目录执行：

```sh
../../.venv_nn/bin/python research/p3l_r11_stable_structure_digits/shared/portable_analysis.py score
../../.venv_nn/bin/python research/p3l_r11_stable_structure_digits/shared/portable_analysis.py cost
../../.venv_nn/bin/python research/p3l_r11_stable_structure_digits/shared/build_site.py
```

前两项只对封存输入重算指标、图和费用归一化；不会改变模型、GT、物理预测或原始逐调用账本。费用脚本只把明确RHS/LU付费事件计一次，不把END累计快照再次求和。portable入口预载归档里的prior与固定cache，不要求恢复大型历史TRAIN。

## 实际执行过的物理入口（仅用于审阅，不自动运行）

- `shared/prepare_digits.py`：构造N/P truth并用原576阶完整模型生成测量。
- `shared/invert_digits.py --dataset N --method FULL256`，P同理；method为FULL256/AGENT_A/AGENT_B/AGENT_C。
- `shared/run_matrix.py`：原始批处理顺序，保留以复核执行故障，**它没有跨job全局LU预留，不能作为下一轮安全runner复用**。
- `agent_a/benchmark.py`、`benchmark_digits.py`、`check_support_birth.py`。
- `agent_b/kernel_qualification.py`及supplement（实际文件见目录）。
- `agent_c/path_campaign.py`，原192节点全局/局部配对。

若今后重新执行物理，应另建独立目录、重新冻结资源，并先接入跨进程原子RHS/LU预留。当前结果必须保留 `FAIL_CAP_EXCEEDED`；不能把现存代码加上保护后声称旧实验当时已有保护。

## 审阅重点

1. `audit/ALL_PREDICTIONS_SEALED.json`在GT评分前已封存72预测。
2. `audit/BLIND_READ_N.json`、`BLIND_READ_P.json`早于解封key及评分；主线程实际调用图像工具查看。
3. GT与在线文件物理隔离，`invert_digits.py`不读取truth或digit_targets。
4. P可表达数字的高latent范数、N的80J截断、同中点属性初值、窗口计时竞争与资源越限不可省略。
5. 公式集新增F-SS-001—012，陈述数学条件、代码与实测资格，未用机器精度数值当新定理。
