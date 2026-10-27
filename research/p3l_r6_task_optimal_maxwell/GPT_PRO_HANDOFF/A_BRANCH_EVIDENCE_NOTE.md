# A分支：条件copula、真实Maxwell响应与未解决边界

本说明只整理已完成A实验，供后续研究审阅，不作最终创新裁决。A实际合计**158144 RHS**：分层/局部参考157656、显式向量复测40、copula控制448；旧分析的阶段小计与旧raw/协议/ledger保持原样。模型范围是固定2D/3D、活动Gaussian截断[-4,4]的原P1子模型，不是Full66边缘后验。[阶段说明](../STAGE_A_TASK_INTEGRATION.md)与[实际补充receipt](../raw/jobs/fixedmaterial_fibre_copula_counterexample/CLOSED.json)可追溯成本。

**条件copula机制。** 固定材料V=QD=v且vR>0时，δR=vR/Q，同材料prior fibre可以保留非退化Q不确定性；把两个精确边缘独立相乘也会改变材料均值。原GL128 mode条件矩给出的独立乘积实材料偏移为REGULAR **.00263820**、PURE_OVERLAP **.01802785**。真实控制各抽32个配对、只置乱一次Q相对于D：保持两个经验边缘，材料却改变。它们是固定材料的quadrature-conditioned coupling controls，**不是总体iid product draws、posterior samples、Gamma/PIT校准或多材料模式概率**。条件推论与适用假设由[根理论文件](../THEORY_AND_COUNTEREXAMPLES.md)负责。

**Full-Maxwell实测。** 新增448 RHS/56个forward states/0 adjoint，10个字节相同材料精确复用。单次置乱的材料mean-shift ROI RMS分别为**.0024849243/.0260784558**，variance defect为**.0006234622/.0661445834**；相对同材料joint reference的complex L2响应差除σ，均值为**15.52605/117.57761**。单次32置乱均值不同于精确边缘乘积均值，应分别阅读。完整q、材料、prediction、reference current/field、raw errors及图源均在[补充报告](../track_a/COPULA_COUNTEREXAMPLE_REPORT_ZH.md)；[copula图](../track_a/copula_analysis/COPULA_JOINT_AND_PRODUCT.png)、[响应图](../track_a/copula_analysis/COPULA_FULL_RESPONSE_DISTANCE.png)只展示控制现象。

**host零切向J与非零任务响应可以并存。** 饱和host平台的原decoder/source切向Jacobian为零，不等于Maxwell对任意材料扰动都为零，也不使整个条件材料任务停止响应。已存最细S3 quotient规则在基准σ的within/between任务导数范数约**.0025931/.0009989**，总范数.0035920；host分支内贡献为零，整体within来自非host分支，between来自源质量迁移。导数坐标为sqrt(2)×source-major pack(y)。同一缓存节点的中央差分与协方差代数相符至约1.56e-9相对误差；这是有限求积规则检查，不是全后验导数认证。[分解CSV](../track_a/analysis/WITHIN_BETWEEN_DERIVATIVES.csv)、[差分CSV](../track_a/analysis/FIXED_NODE_FINITE_DIFFERENCES.csv)。

**仍未建立。** global tail/遗漏质量上界没有取得，4D保持NOT_RUN，材料均值与分支/事件概率的可信同精度成本分离为NOT_ESTABLISHED。prior归一化接近1的全局网格仍会漏掉窄likelihood：REGULAR最细source/quotient相对局部GH21的均值误差约.524/3.585个材料SD；overlap全局分支质量也错配。局部GH15→21、内层GL64/128/256及latent功能的稳定都只标Q1/Q2；不因新的copula控制而提升为全局概率、连续PDE或Full66认证。[逐任务资格](../track_a/analysis/TASK_QUALIFICATION.csv)、[独立fibre QA](../track_a/COPULA_FIBRE_REFINEMENT_QA.json)、[先前latent Q说明](../track_a/analysis/POSTHOC_LATENT_Q_REPORT_ZH.md)。
