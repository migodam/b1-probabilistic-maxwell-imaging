# R12 原要求—实际证据—剩余工作

本表是恢复导航，不代替[用户原始任务](sources/user/R12_ORIGINAL_TASK.txt)。更新于2026-11-12；科学结论以封存与评分回执为准。

|用户要求|实际证据位置|状态/还需做什么|
|---|---|---|
|完整继承Pro理论/toy与两Chat、R4–R11|sources/pro；audit/SOURCE_ACCESS.md；THEORY_DERIVATIONS.md|原始ZIP保留；Chat有界访问限制披露；最终查新需保留缺口|
|共享实材料双场投影、固定背景A投影、DR/AP|src/field_current.py；METHODS_AND_ABLATIONS.md|已实现标量/三分量、共享材料、多源mask、grouped优化；不称新DR定理|
|四个Pro对象与成熟CIE/Dual/global/A+/G19|experiments/pro_reproduction；delegated/baseline|重放与40成熟对照已完成；toy结果不代替新确认|
|R11真实坏终点逃逸及同预算控制|RESULTS/basin_remote；BASIN_ESCAPE_REPORT_ZH.md|12远端比较已封存评分；body CPU前缀/超额开销明确；完整states归档中|
|多样独立高/混合对比二维、成本公平|delegated/confirmation|24对象×8=192已封存；另6开发×8=48；离线评分/统计与root回读均完成；不是全部方法都已证胜出|
|NN→direct/CIE/DR、复用小网络|delegated/neural；NEURAL_AND_SELECTION_REPORT_ZH.md|16开发比较完成；新确认含NN→direct；区分训练prior信息，非所有NN组合均独立确认|
|真正vector3D、不同网格与高对比|delegated/vector；delegated/vector_high|36 moderate与12 high/跨网格运行完成；coarseDDA资格与模型误差限制披露|
|Fresnel、Austria、原论文物理合同|delegated/external及external_*；BENCHMARKS.md|实际实测/150→32与matched32诊断完成；不是原作者算法完整复刻；保持全部负结果|
|失败后诊断修正而非任意gate结束|RESEARCH_LOOP.md；external_continued/selection/physics/matched32；SECONDARY_SELECTION_FREEZE.json|续跑、物理匹配诊断、grouped提速及预注册次选择已实施；综合分析已完成，完整负结果保留|
|同噪声球可行错误材料的选择限制|THEORY_DERIVATIONS；vector_high分解；confirmation主次输出|次规则须独立汇报，不能隐藏unavailable或替换主Phi结果|
|成本/所有失败/LU/RHS/FFT/投影|COST_LEDGER.md；各receipt；confirmation执行成本|确认与评分已回填；运输/打包继续更新；缺失保持null，不重复相加|
|先例和原创性|PRIOR_ART_AND_NOVELTY.md；literature|已读DCTMC；部分全文访问受限；需诚实定位，不能宣称首次反射/非物理状态|
|中文报告、英文论文、图表、原始数据、GitHub/网页|各分支报告已成；根RESEARCH_REPORT_ZH/PAPER_DRAFT/CHAT_HANDOFF已整合并恢复验收|科学文稿完成；归档发布未完成，须上传既有仓库并核验链接|

尚未完整覆盖的实验轴应明确列出：独立24主要是固定1%噪声/12源全孔径/24²匹配离散；其它孔径/源数见已公开NN开发对象，跨网格见Austria和vector。不能将分散开发试验说成所有轴的全面独立统计验证。是否继续新增实验须从确认结果和实际机制出发，不按GT筛样。

## 2026-11-12最终发布

完整原始证据已归档至[GitHub Release](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-max-struct-r12-20261111)，12个逻辑包、89分片，共11,416,183,356字节；逐个公开回读通过。新增综合评价与结构复核已纳入，不改变冻结原指标或在线规则。所有新增大文件位于外置盘。
