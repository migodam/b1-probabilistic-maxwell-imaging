# MAX-STRUCT 结构研究交付

日期：2026-11-08。历史证据固定仓库提交：24ab7e85d750a477ab337e03cdb1a77d5c9ff2c3。

## 文件

- MAX_STRUCT_RESEARCH_REPORT_ZH.md：完整中文研究报告、推导、证据分层、18组文献与创新审查。
- CODEX_EXPERIMENT_PROTOCOL.md：下一轮实施协议；尚未执行其中的论文级实验。
- feasibility/：本次已经运行的4对象 scalar-TM 小模型、后续机制控制、必要代码审查与全部结果。
- figures/：全部材料、真实残差轨迹及端点路径图。
- MANIFEST_SHA256.json：交付文件校验值。

## 本次真正做了什么

已审查已有 R4–R11 及相关研究，没有修改远程仓库，没有重跑历史大规模反演。新做了100材料未知量、6源、24接收的4个固定有限VIE对象。

主探索中，反射在高对比环/壳上明显优于该 direct/AP 实现；桥未到噪声球、物理corrector降低数据残差却使材料更差。事后控制显示从实际错误空气终点也可以得到好环/壳候选。所有公开输出只按真实数据项选择，不按GT选择。

这些结果不是盲测，不是与成熟CIE/G19的公平性能比较，也未建立连续物理、跨网格、3D或等CPU优势。原物理LU/RHS之外的固定背景求解费用必须计算。

## 复现环境

记录的运行环境：

- Python 3.12.14
- NumPy 2.3.5
- SciPy 1.17.0
- Matplotlib 3.10.8（仅绘图需要）

软件版本以已运行环境为准；没有宣称所有旧版本完全复现相同数值轨迹。保存的原结果可直接审查，无需先运行。

解压后，在本目录依次执行：

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python feasibility/consensus_reflection_toy.py
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python feasibility/followup_controls.py
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python feasibility/toy_qualification_audit.py
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python feasibility/make_figures.py

脚本会覆盖同目录同名结果；若希望保留原始运行，请先复制整个目录，再在副本中运行。所有数据和代码均自包含；运行不需要GitHub访问或网络服务。

consensus_reflection_toy.py 先生成固定对象/噪声并运行主试验。followup_controls.py 依赖主试验的NPZ，并标记为POSTHOC_MECHANISM_CONTROL。toy_qualification_audit.py 只检查新增实现的具体风险，不重跑3000步优化。make_figures.py 读取保存结果生成图。

## 代码范围

此代码是固定100维全观测二维模型的可行性原型，不是下一轮生产求解器。特别是实际receiver masks、相关噪声、端口不可达/奇异极限、大规模矩阵自由投影、共享矢量材料及复杂材料均须按完整报告和协议接入，不能仅更改n后宣称已支持。

NPZ包含GT及clean fields，因为本次是公开探索研究。下一轮锁定测试必须把生成器、在线求解器和评分器隔离；不能直接复制本原型的数据访问方式作为盲测设计。

没有用GT生成跨盆地proposal、选择最优候选或在线停止，但脚本记录GT诊断，因此不称完全盲化。相对材料L2不是历史amplitude NRMSE。

## 推荐下一步

先做协议M1：同R11真实错误终点，AP/DR/原物理继续/同预算扰动控制；同时保留M2球内固定点no-go。只在这些机制明确后解释独立24对象和外部锚点，不先扩网络或扫Green rank。

