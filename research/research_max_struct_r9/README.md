# MAX-STRUCT R9：材料残差传播、噪声约束与跨维度反演

**本轮实验完成。** A/A+保留高对比开发正例，但尚无稳定超越成熟CIE的恢复增量。A+相对原A降低部分完整反演成本；公开噪声相容的DPCS修复让普通CIE与A+共同受益。

[简明交接](CHAT_HANDOFF.md) · [中文报告](RESEARCH_REPORT_ZH.md) · [公式](THEORY_DERIVATIONS.md) · [先例](PRIOR_ART_AND_NOVELTY.md) · [英文初稿](PAPER_DRAFT.md)

完成的证据池分别为：

- 113次二维机制反演：47公开开发、66次六个新材料比较；六个独立对象才是确认分母。
- 原NIE Austria两种材料条件、16项Windows对照；另存9项Mac迁移前结果及成本。
- 29次小型全矢量三维反演。
- 原CIE 3λ四立方体九臂，27000个实材料未知、60电偶极源。
- 三个三维起点各八步物理校正，以及Mie、积分Green和界面误差诊断。

Austria ε3.5上普通CIE明显优于direct，原尺寸三维CIE明显优于LSIE；A+没有进一步的稳定增量。三维两个辅助端点经八步校正都由约.423误差降至.249，但Austria更长校正损害材料质量。两种相反停止问题不能靠统一增加迭代解决。

[经典二维](CLASSIC_BENCHMARKS.md) · [小三维](VECTOR3D_EXTENSION.md) · [原尺寸九臂](delegated/original3d_feasibility/ORIGINAL_CM_NINE_REPORT_ZH.md) · [短校正](delegated/3d_extension/physical_refinement/REPORT_ZH.md) · [物理误差](VECTOR_DISCRETIZATION_DIAGNOSTIC.md)

二维远端四P核、原尺寸三维resident CUDA，未为提速改变源/材料/观测合同。R1–R8不覆盖。GT未进入在线选择；本轮无新实测或NN训练。

[数据合同](BENCHMARKS.md) · [远端与CUDA](REMOTE_COMPUTE.md) · [最终成本](COST_LEDGER_EXTENSION.md) · [失败](FAILURES.md) · [复现](REPRODUCE.md) · [完整数组Release](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-max-struct-r9-20261106)

Git保留代码、紧凑结果、报告和图；大型状态、长轨迹和原始数组在同仓库Release无损分包。清单列路径与字节数，不要求SHA核验。复现先读REPRODUCE，提取数组后使用新的输出目录，避免覆盖冻结结果。
