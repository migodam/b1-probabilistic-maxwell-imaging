# MAX-STRUCT-R4 邻近工作与证据边界

日期：2026-11-03。范围：材料方向有限步有理/ROM 搜索；完整残差 Hessian 与负曲率；端口状态残差、变量投影/WRI、T-matrix 反演。仅文献审计；没有数值运行、GT 读取或历史结果修改。最终数学、新颖性和路线裁决由主线程负责。

## 可直接用于主报告的结论

1. **材料方向 ROM 线搜索已有直接先例。** Besset 等的 MOR-Tₗ 用方向 Fréchet 导数构造扩展基，重复求解同一背景算子的不同 RHS，并在新材料参数上求解 reduced system。其 §5.2 显式搜索 `θ₀+αδθ₀`。因此“用高阶响应/RHS 构造便宜有限步线搜索”不能独立承担原创性。它未在所读部分使用 Maxwell 接收端可见极点进行方向选择；这项差异仍需明确实现和验证，不能由检索空缺认证。([P1](https://arxiv.org/html/2505.00709v1))
2. **完整 Hessian 的残差二阶项、二阶伴随和利用负曲率均已有。** Pratt 等与 Métivier 等分别给出完整 Hessian/高效计算；Adriaens 等在 FWI 中沿负曲率走到信赖域边界；Royer 等给出通用 Newton-CG 的显式负曲率步。若 R4 只是把这些机制用于现有 Maxwell 目标，应以适配、诊断或对照描述。([P4](https://doi.org/10.1046/j.1365-246X.1998.00498.x), [P5](https://doi.org/10.1137/120877854), [P6](https://doi.org/10.1190/segam2021-3581935.1), [P7](https://arxiv.org/abs/1803.02924))
3. **端口相关 residual geometry 已有很近的文献。** CC-CSI 明确把状态误差映到测量空间；WRI 精确消去状态；van Leeuwen–Yang 明确得到材料相关 Gram 加权的数据残差。R4 的 `S A(m)⁻¹r` 是否形成有用的新目标，要与这些已知结构逐式区分，不能只据“measurement-aware”命名。([P8](https://doi.org/10.1088/0266-5611/32/1/015007), [P9](https://arxiv.org/html/2403.15292v1), [P10](https://doi.org/10.1109/TAP.2017.2673758))
4. **T-matrix 吸收多次散射及用其进行反演已有。** DCTMC 直接使用 `T=(I−VΓ)⁻¹V` 和外部数据 `ATB=Φ`；Jakobsen–Ursin 使用 T-matrix 的 BI/DBI 反演。普通分块、重写或局部消元不能作为新的通用数学。由观测选择可负担的耦合块并取得材料收益，是可以单独检验的具体差异，当前文献审计不替代该检验。([P11](https://doi.org/10.1103/PhysRevE.94.043317), [P12](https://doi.org/10.1088/1742-2132/12/3/400))

## 最接近的证据账本

A = 已读可定位全文；B = 原作者/出版社摘要；C = 仅身份。完整机器可读字段见 `EVIDENCE_LEDGER.json`。

|ID / 原始文献|证据位置与已知机制|R4 需保留的边界|
|---|---|---|
|P1 Besset, Barucq, Djellouli & Frambati (2025 preprint), **MOR-T L: A Novel Model Order Reduction Method for Parametrized Problems with Application to Seismic Wave Propagation**. [arXiv](https://arxiv.org/abs/2505.00709)|**A**，§3 Eq3.1–3.14、§5.2 Eq5.2；本地 `mor_tl2025.txt`。方向高阶导数、QR/SVD 扩展基、同算子多 RHS、材料路径上 reduced forward 与步长搜索。|该展示的观测由同一条已给定方向上的参数生成；搜索 200 个非负步长。它是线搜索能力证据，不是自由材料未知恢复。没有在所读部分发现接收端极点驱动方向选择。|
|P2 van der Sijs, El Gawhary & Urbach (2024), **Padé approximants of the Born series of electromagnetic scattering by a diffraction grating**. [DOI](https://doi.org/10.1103/PhysRevA.109.033522)|**B**，APS 摘要；Crossref 身份。对强电磁散射的发散 vector Born 级数作 Padé，保留逐阶多次散射解释。|正演光栅研究；不等于材料逆问题的有限步接受规则、跨盆地证书或 GT 恢复。它直接排除“用 Padé 处理强散射”这一宽泛原创表述。|
|P3 Druskin, Moskow & Zaslavsky (2021), **Lippmann–Schwinger–Lanczos algorithm for inverse scattering problems**. [DOI](https://doi.org/10.1088/1361-6420/abfca4), [预印本](https://arxiv.org/abs/2101.12317)|**A**，§2 Eq2.3–2.7，PDF p4；`lsl2021.txt:160–225`。测量驱动 resolvent/Galerkin ROM 与 Lanczos 内场近似后代入 LS 反演。|其 ROM 自变量是频率/谱参数，数据合同含多谱点，实点还含导数。不能视为本项目单频、材料步长 t 的同一算法，也不能从单频数据凭空构造该 ROM。|
|P4 Pratt, Shin & Hicks (1998), **Gauss–Newton and full Newton methods in frequency-space seismic waveform inversion**. [DOI](https://doi.org/10.1046/j.1365-246X.1998.00498.x)|**A**，可识别出版社索引全文 §The Gauss–Newton method / §The full Newton method、Eq32/35/42；Crossref 身份。明确分离 GN 与被舍弃的精确 Hessian 项，并给 backpropagation 计算和二阶散射解释。|不是新 Hessian 恒等式。既有声学案例的作用不能无条件移植到当前 Maxwell、多源归一化、box/prior。未保存该文完整 PDF，仅保存网页检索收据。|
|P5 Métivier, Brossier, Virieux & Operto (2013), **Full Waveform Inversion and the Truncated Newton Method**. [DOI](https://doi.org/10.1137/120877854), [HAL全文](https://hal.science/hal-00763702/document)|**A**，HAL PDF p6 Eq1.7–1.8、§2.2、p16 §3.4；`metivier2013.txt:245–278,835–852`。完整残差 Hessian、矩阵无显式构造的二阶伴随 HVP、同 GN 比较。|其 §3.4 遇负曲率停止内层 CG；首步遇到则退回梯度。不要将这条实现规则说成主动负曲率逃逸。2017 SIAM Review 是重刊，不另算独立实验证据。|
|P6 Adriaens, Métivier & Geuzaine (2021), **A trust-region Newton method for frequency-domain full-waveform inversion**. [DOI](https://doi.org/10.1190/segam2021-3581935.1), [机构全文](https://orbi.uliege.be/bitstream/2268/288603/1/ADRIAENS_TRNM_FWI.pdf)|**A**，PDF p2 Trust region、p3 Algorithm2；`trust_fwi2021.txt`。明确沿负曲率方向到信赖域边界，并用实际/预测下降比调整半径。|直接 FWI 先例，强于仅发现负曲率的相邻文献。论文声学 Marmousi 结果不保证任意错误严格局部极小点可逃逸，不是 Maxwell GT 证书。|
|P7 Royer, O’Neill & Wright (2018 preprint; 2019 online), **A Newton-CG algorithm with complexity guarantees for smooth unconstrained optimization**. [arXiv](https://arxiv.org/abs/1803.02924), [DOI](https://doi.org/10.1007/s10107-019-01362-7)|**A**，§3.3、Algorithm3，PDF pp11–12；`royer2018.txt:546–615`。显式负曲率步或 damped Newton 步、回溯与近似二阶最优性。|通用无约束、光滑与 Hessian 正则条件；不能直接套到 box 边界或非光滑 TV，也不保证找到全局或真实材料。|
|P8 van Leeuwen & Herrmann (2016 volume; 2015 online), **A penalty method for PDE-constrained optimization in inverse problems**. [DOI](https://doi.org/10.1088/0266-5611/32/1/015007)|**A**，§4 Eq13/19–22；`penalty2015.txt:365–388`。固定材料精确消去线性状态；约化 Hessian 为联合 Hessian Schur 补。|状态消元、envelope gradient、Schur 补不是新通用机制。有限罚的原始物理验收与材料恢复仍须独立报告。|
|P9 van Leeuwen & Yang (2024 arXiv identity), **A data-driven approach to PDE-constrained optimization in inverse problems**, 已读 v1. [v1全文](https://arxiv.org/html/2403.15292v1)|**A**，§2.3 Thm1、§2.4 Cor1/Eq10，HTML introduction 亦列出 `(I+ρ⁻¹G(c))⁻¹` 数据残差权重；`vanleeuwen_yang2024.html`。|Gram 是材料相关的伴随场/状态结构。其 data-only Gram 估计和凸性极限具有额外假设；不可直接移植为任意非同置、单频有损 Maxwell 结论。按 arXiv ID/version 识别，不能混合后续改题版本。|
|P10 Sun, Kooij, Jin & Yarovoy (2017), **Cross-Correlated Contrast Source Inversion**. [DOI](https://doi.org/10.1109/TAP.2017.2673758), [作者全文](https://arxiv.org/abs/1906.10864)|**A**，§3.2 Eq11–13，PDF p5；`cccsi2017.txt:239–269`。在普通 data/state 项外加入测量空间的 constitutive mismatch。|其传播算子是背景测量算子，输入为一步 constitutive 预测；不等于 R4 当前材料的完整 `S A(m)⁻¹r`。这一区别支持逐式比较，不能据此宣称全新类别。|
|P11 Levinson & Markel (2016), **Solution of the nonlinear inverse scattering problem by T-matrix completion. I. Theory**. [DOI](https://doi.org/10.1103/PhysRevE.94.043317)|**A**，§III Eq8–11，PDF p4；`dctmc2016.txt:225–265`。完整 transition matrix、外部测量线性约束、在迭代中放宽材料严格对角性。|它估计完整 T 并向局部材料拉回；R4 按可观测耦合选择小块不是同一实现。论文理论部分不能替代成像效果证据；不应把放宽局部性的中间解直接视为合法材料。|
|P12 Jakobsen & Ursin (2015), **Full waveform inversion in the frequency domain using direct iterative T-matrix methods**. [DOI](https://doi.org/10.1088/1742-2132/12/3/400)|**A**，可识别出版社索引全文 §3.1–3.2、Eq50–54 与77–80；Crossref 精确题名记录。T-matrix BI/DBI，多次散射更新，指出 DBIT 与 GN 关系。|“T-matrix 吸收反馈回路并用于逆散射”已知；该文的声学、多频、多源合同及具体更新不等于本轮选块法。出版社直接 open 失败，索引全文存在，完整 PDF 未保存。|

## 供主线程裁决的差异，不是原创性结论

- **有限步路线：** 真正需要解释的是方向如何从当前失败端点和可观测高阶响应产生、ROM 对候选的误差如何控制、与普通精确线搜索/已有材料 ROM 的差别。P1 本身已覆盖“低成本多步长评估”；P2/P3 已覆盖强散射有理逼近与 resolvent ROM。没有检索到逐式相同实现 ≠ 不存在先例。
- **曲率路线：** 需要区分 GN 看不到的负曲率、box 的合法切锥、以及 Hessian 已半正定的真实局部极小。文献已有二阶搜索；材料是否改善只能由原始物理目标和独立材料指标确定。不能把 residual 减少写成 GT 证书。
- **端口路线：** 先回答完整 current 消元后与 WRI/Gram-weighted profile 的关系，再谈新的优化路径。某个权重忽略 receiver-null 方向或随材料变化，其后果应由父线程公式裁决；本文不代替证明。
- **块散射路线：** 局部 T/Schur 的代数已知。可检验差异是如何用可观测反馈选块、保留什么隐藏反馈、控制多少 RHS/因子分解，以及同预算原生 CIE 对照。完整逆放入一个大块只是转移计算。
- **评价合同：** 本轮不得借 P1 的 planted-ray 实验方式宣称自由材料恢复，也不得借 P3 的谱域数据合同偷偷加入测量。R3 已公开的 Möbius/rank-two 与有限-rank 有理路径是项目既有知识。

## 检索、版本与缺口

采用 scholarqa-research 的 evidence-first 流程，已读 novelty/retrieval/reporting 参考。首选 ScholarQA CLI：3 查询、6 个 snippet/paper 操作全部 HTTP429，exit1、零候选；一次 batch verify 亦顶层429。没有暴露 Semantic Scholar MCP。依父任务授权，以已识别作者/机构/出版社原文补缺，并用 Crossref/arXiv 核对身份。429 是检索不可用，绝不是没有先例。

`receipts/SEARCH_LOG.json` 保存24条 web 发现查询、3条 ScholarQA 查询、入选/排除规则；`receipts/scholarqa_*.json` 保存失败；Crossref 原始记录与下载/页码记录保存于 `receipts/`。第三方全文和原始网页结果只在已忽略的 `local_only/`；不应跟随项目公开。

P1 的 DOI `10.1016/j.cam.2026.118249` 已有登记，但 Crossref 当前卷期日期为 **2027-04**；本报告仅以已公开2025预印本作为证据，不称其2026正式卷期。第一次该 DOI 的 Crossref 查询429、后续单请求成功，均不隐去。P5 的2017作者页面404，改用成功的2013 HAL原文。P9 的版本题名与网页自动日期存在混淆风险，按已读v1和arXiv身份标记。

本审计不是穷尽或系统综述；没有执行完整 Scideator/SPECTER/40-example 分类器，没有给出二元或全球原创标签。尚未核清：与 R4 完全相同的单频 Maxwell、port-visible material-path pole selection 是否有直接先例，以及一个最终冻结实现能否产生独立材料收益。这些是待解决问题。

