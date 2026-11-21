# 整个研究项目从哪里进入

项目外置工作区 Inv_SLAM 包含多个独立仓库、工作树和早期交付。路径以下均以工作区为根；公开库链接单列。**数据集登记表按可查到冻结对象和观测的位置建立**。这里再给只有论文、协议、诊断包或专门训练材料的系列一个可追溯入口。名字相同的文件夹、打包时附带的实验副本和工作树数量均不能作为样本数。

|研究问题或系列|实际入口|沿用价值及证据边界|
|---|---|---|
|MAX-STRUCT、FBG、A17-X、Z-ALG、Formula Atlas|B1 的 `research/research_max_struct_r*/`、`research/research_fbg_r3_experiments/`、`research/maxwell_formula_atlas/`；工作区 `research/delegated/max_struct_r1/` 保留分轮冻结工作树|按轮用相应代码、算子及观测；R16 后端中装载的早期资料分别登记原队列。Atlas 管公式的假设与来源，手册管具体实验|
|初期 A 系列和论文资料|工作区 `communication/A2_RESEARCH_REPORT_ZH.md`、`Gaussian/README.md`、`Gaussian/A6/deliverables/`、`Gaussian/A8/deliverables/A8_RESEARCH_REPORT_ZH.md`、`Gaussian/A9/A9_EXPERIMENT_REPORT_ZH.md`；`research/delegated/a2_literature/data/` 保存原始论文测量|早期报告可用于找机制、原数据和失败证据。重新使用其中测量时登记具体文件、来源和采集合同|
|A10–A18 的投影、TSOM、current 与交换| `Gaussian/A10/protocols/`、`A11/protocols/`、`A12/A12_FINAL_RESEARCH_REPORT_ZH.md`、`A13/A13_EXEC_FINAL_REPORT_ZH.md`、`A14/paired_tsom_research/`、`A15/DECISION_REPORT_ZH.md`、`A16/CURRENT_SELECTION_R1/`；在同一 Gaussian 下进入 `A17/DIGITS_2D_20261030/`、`A17/SHAPES_3D_BASELINES_20261030/`、`A18_C/a18-physics-native-current-proposal/`|A17 的端点电流秩、盒约束与全梯度等时结果要同时读取。A18 的 current 提议、闭环检查和训练材料应带上冻结版本；仅凭同名方案不能接到新的 GT 分母|
|A19–A32 的反馈与风险| `Gaussian/A19/a19-task-adaptive-feedback/`、`A20/a20-r1-mechanism/`、`A21/a21-two-sided-anatomy/`、`A22/a22-r1-subspace-separation/`、`A23/A23_research/`、`A24/a24-source-anchored-om-feedback/`、`A27/a27-feedback-failure-attribution/`、`A32/a32-material-risk-directed-feedback/`，其余沿实际 A 目录查|适合查条件电流、两侧子空间、材料风险、消融与已知负例。A24、A25 等目录还带后续打包子目录，应查实际父轮和冻结协议|
|B1 的 P3L / GlobalDual / Z-ALG| `Gaussian/B1/b1-p3l-r11-stable-digits/`、`b1-p3l-r12-bilateral-fixedpoint/`；B1 的 `research/p3l_r11_stable_structure_digits/`、`p3l_r12_bilateral_fixedpoint/`、`research_z_alg_asymmetric_validation/`、`z_alg_host_fold/`|自由像素重建、条件材料图册、已知宿主折叠属于不同逆问题。材料非唯一性结论引用具体观测和允许材料类。数字的生成标签或文件名不能作为形状先验|
|B1 的 GN/NN 与 wavelet|B1 的 `research/research_gn_nn_completion_r3/`、`research_gn_nn_preservation_r4/`、`research_gn_nn_feedback_wavelet_r5/`，及该轮原始训练/验证资产|R5 有 128 个训练 parent/256 视图、78 个验证 parent/156 视图；已选的 36 parent/72 视图 corrector 结果保留 VAL 选择身份。辅助场能量、原 Maxwell 残差、材料恢复分别核算|
|B1 的 A/B/T/W 与 Q 系列| `Gaussian/Research_A/`、`Research_B/`、`Research_T/`、`Research_W/`；`research/q5_remaining_v1/`、`q6_selfcal_closure/`、`trispace_self_calibration/`、`agentic_input_trispace/`|背景、校准和 nuisance 未知时，观测含义及信息预算会变。公开通信、成像或标定结果以对应协议为准。Q6 的 Fresnel 资格单列于 `FRESNEL_PROTOCOL_STATUS.md`|
|生成式与 Flow Matching 项目| `Gaussian/C1/c1-physics-state-generative-imaging/`、`Gaussian/C1/c1-joint-p2/`、`Flow Matching/A13/`、`A14/`、`A15/START_HERE.md`|复用训练前需查对象划分、采样后验意义、权重和物理修正。Flow Matching A15 与 Gaussian A15 使用不同工作目录和任务合同|

公开资料从[B1 研究索引](https://github.com/migodam/b1-probabilistic-maxwell-imaging/tree/main/research)、[A17 原仓库](https://github.com/migodam/a17-maxwell-certified-exchange)、[Formula Atlas](https://github.com/migodam/b1-probabilistic-maxwell-imaging/tree/main/research/maxwell_formula_atlas)进入。需要找未逐对象登记的旧资产时，先沿表中原始报告、协议和生成器追到数据；实际纳入新比较时，为所用子集补一行来源和观测映射。这避免将只含论文建议的包误作已经运行的证据。
