# R5失败、限制及未运行记录

## 保留的科学负结果

- current-only未复制旧正例：NRMSE .558514；完整记录、corrector和成本均保存。
- material-only虽复制旧正例，但42cb/63ab仍失败；63ab材料误差改善伴随原Phi/残差变差，不由GT选它。
- early-only与late-only自然阶段切换均失败。共同corrector1000步耗尽，不能称已收敛到严格局部极小。
- 全域cross符号、同相模态占比、connection存在、冻结transport大小均未成为可验证的泛化预测器；未拟合阈值救解释。
- R4六确认共同4/6、32网格比较优势消失、平移平均失败、继续corrector未修复全部继承，不用R5旧正例覆盖它们。
- 三维例仅表明dyadic极化结构不同，未验证任何三维成像优势。
- 新独立盲测、R5九数字、Fresnel/eispy2d反演、自由材料三维成像均NOT_RUN；具体科学理由见BENCHMARKS和RESEARCH_LOOP。

## 实现/环境及输出错误

1. 首次trajectory诊断遇NumPy局部索引错误，在异常前已消耗1LU/24complex RHS；已修复，失败源记录与代价保留，计入439LU/9672RHS诊断总额。失败前辅助操作计数未返回，保持null。
2. review指出普通/self解析零transport不应报告浮点残余方向；仅用保存数组修正文案/规范零，raw roundoff保留，无重跑物理。
3. review指出exact-Phi曲率减data-GN还含prior；现已分离prior、exact-data及residual-feedback，无重跑物理。
4. 父线程一次只读评分检查假定了错误JSON嵌套字段，报KeyError；读取实际结构后修复分析，物理费用0。一次系统python3无NumPy，改用现有venv，未安装包；费用0物理。零物理人工检查的独立耗时未计量。
5. ScholarQA发现检索遭429，结果为空；核验9条metadata成功，额外核验遇429后使用作者全文/官方网页/Crossref。没有把空检索称穷尽查新，完整stderr/receipts保留。
6. 图像标题留白和换行修正，无数据变更；发布可达性检查与科学资格分开。

## 封存与后续文档演化

STAGE_I_SEAL在第一次评分前列出693个文件，独立评分当时验证其bytes/mtime未变化。此后THEORY_DERIVATIONS仅追加第7–9节；原6480字节前缀保存在audit/THEORY_DERIVATIONS_AT_FIRST_SEAL.md，重建保存时刻和方式如实写入POST_SEAL_DOCUMENT_EVOLUTION.json。不是声称该副本在最初封存时生成，也不声称最终理论文件仍符合旧mtime。预测、配置和第一次评分依据没有改写。

其余实验各自freeze/score独立记录；后续post-score发展不变成盲确认。没有SHA验证，没有新安装服务、CUDA基础设施或外部GPU占用。


发布整理中一次工具编排变量名错误在启动shell前报ReferenceError，随后修正；物理调用0，未改变预测。
