# 可直接交给Codex的执行Prompt

你负责B1下一阶段的自主理论验证与实验研发。先完整阅读随附的 `B1_NEXT_RESEARCH_REPORT_ZH.md`、`checks.json`、`verify_theory.py`，再执行。不要只写计划。

## 研究入口

仓库：`https://github.com/migodam/b1-probabilistic-maxwell-imaging`
历史基准commit：`db5aacc9b7f1cc46118a216657df9af8d3d9d1f7`。

先读 `START_HERE.md`、`RESEARCH_REPORT_ZH.md`、`THEORY.md`、`RESULTS.md`、`DECISION.md`、`GATE_DECISION.json`、`RESEARCH_LOG.md`、`configs/protocol.json`、`research/delegated/probability_review/REVIEW_NATIVE.md`，再查 `src/b1/physics.py`、`inference.py`、`refine.py`、原始PRIMARY/refinement摘要、banks和资源账本。

当前工作树可能已有比该commit更新的B1实验。先识别运行中任务、最新研究包、现有用户授权、结果和预算；历史commit只作只读baseline，不能reset或覆盖最新工作。用户的B1实验已经在进行，不得另起一个重复的variable-radius/posterior campaign，不得杀掉现有任务。

## 目标

最终是统一材料空间中的多假设posterior，不是按0孔/1孔/2孔建立固定生成专家。拓扑可用于输出事件和事后评价，但不是强制生成标签。B1负责物理参考、后验正确性和辨识理论；C1的learned proposal以后按瓶颈决定，不启动大型NN。

核心假设是：nuisance-aware nonlinear material secant能够解释结构信息，以及哪些物理近似必须纠正才能保护材料事件概率。允许实验否定这个假设，不得只调大噪声制造成功。

## 先执行E0

1. 在隔离目录/branch工作，保留旧结果、失败和ledger。
2. 复核真实prior measure、real/complex噪声约定、gain边缘化logdet、evidence、Top-K剩余质量和unknown mass语义。
3. 专门核对distinguishable/H0的Re(epsilon)：原生variance约3.2290906e-6，细化约7.3325252e-4。建立独立的连续CDF/均值/SD/quantile收敛门槛，不能只看topology mass。
同一原生case的q05=2.44357576、q95=2.55644933、mean=2.50002435还要求variance至少3.1851109e-4，已报值违反该必要下界约98.64倍。检查CDF-bin与节点moments的不同近似，统一summary的概率target；对原生与细化都加入coherence gate，允许合理数值容差，不把这一问题误称evidence归一化错误。
4. 区分四个treams旧oracle数组缺失与新refine代码已经保存数组；不重跑来掩盖历史缺口。
5. 执行随包零物理代数脚本，加入高ESS漏模式、prior-only posterior、坐标变换、event partition、metadata回放等必要测试。任何现有test先检查是否调用physics；零预算阶段不得直接运行audit/pilot/refine。
6. 输出E0_AUDIT、可复用资产、当前运行状态、下一步科学假设和资源计划，然后在实际授权范围内继续，不频繁向用户提问。

## 按已有授权接入E1

已有variable-radius实验尽量复用。加入raw full-Mie secant、epsilon补偿、epsilon+gain补偿、局部投影预测与完整marginal posterior的比较；profile minima不等于evidence。比较小孔/大孔、SNR、prior边界和多频；新频率必须重新计费并核验multipole收敛。

若最小统一表示尚未覆盖且预算允许，使用报告E1-U的连续径向implicit field，roots转为分层球radii，经现有Forward.run调用。不先枚举H；拓扑只由生成材料事后计算。清楚标记这是径向校准台，不是开放一般三维支持。

算法至少包含可信积分/独立reference、generic tempered inference、物理引导并经full-likelihood correction的方案。在低维台架优先把概率算对；不得为展示NN而训练。

## 预算与权限

- 历史384/384已封口，不能降低、删除、改名ledger或用新checkout假装历史未消费。
- 本prompt默认 `NEW_STAGE_EXTRA_AUTHORIZATION = 0`：可立即做零新增物理E0；已有阶段有明确用户授权时，可在其剩余额度和原scope内接入，不自动追加。
- 建议新的P1独立预算为4,000次材料-频率评估、CPU优先、GPU训练0、host wall保护上限2小时；这是建议，不是本prompt默认追加授权。只有用户将本配置明确改为4000或另有相应授权记录，才开启该独立阶段。
- P1建议分配：400物理/机制资格、2000reference与算法、400fresh parents/truth、400独立核验、800失败与refinement预留。所有setup、adjoint、额外频率、重复、失败同样记账；多RHS成本另外记录。
- 不授权E2非同心campaign、大型NN、访问旧locked truth、自动公开发布、调用DeepSeek、修改A17冻结稿。
- 无授权或预算不足时完成所有可做的分析、补丁和设计，输出明确停止原因；不要反复询问，也不要暗中越界。

## 自主迭代规则

每轮明确假设、单一主要改动、匹配基线、实际成本、正负结果和下一步决定。可以在授权内改变sampler、proposal、temperature、积分策略；改prior/support/likelihood或评价目标时新建protocol，不能把不同后验混在一起。

选择性full solve必须保护重要event odds和CDF；未获uniform error bound只能标empirical。exact MH correction只保证声明的numerical target，不等于continuum/现实认证。proposal全support与高ESS不等于没有漏mode。真posterior样本的事件频率可以用，但必须有误差和mixing检查。

若r posterior单峰或ridge，诚实报告；不要人为制造“多个modes”。若物理引导无额外解释力或全成本收益，停止其优越性主张。若support缺真值，先修表示，而非继续优化sampler。

## 交付

生成START_HERE、中文研究报告、THEORY_STATUS、REFERENCE_AUDIT、PRIOR_AND_NUISANCE_SENSITIVITY、METHOD_COMPARISON、RESOURCE_SUMMARY、FAILURES、GATE_DECISION及可复现代码/必要原始数组。避免为文件清单创建空壳。

Top-K图必须附事件定义、代表材料、质量及误差、条件不确定性、prior/后验odds、剩余已建模质量、未探索/模型外状态；不把Top-K重归一化为1，不以漂亮图替代概率参考。

最后做清晰研究裁决：继续当前机制、改做active sensing、只保留概率benchmark、还是停止某路线。不要以“还需要再调模型”无限延长。下一步动作必须受真实结果和剩余授权约束。