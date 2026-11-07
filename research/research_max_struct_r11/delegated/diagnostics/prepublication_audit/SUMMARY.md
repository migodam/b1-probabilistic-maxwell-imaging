# R11发表前只读审核

发现4项需根线程决定的实质表述/披露问题，均为P2，详见[findings.csv](findings.csv)：

1. 网页将1728个材料未知量写成“1728个独立实材料”，应避免与独立对象数混淆。
2. 网页“付费组合没有增加成功数”缺比较对象：两组合均5/16，单CIE为4/16；只能说没有优于CIE+Dual控制。
3. 主48的公开状态段遗漏132成功后4attempt工程中断/恢复及unknown端到端CPU边界。保留主48、分开尝试成本，并披露受影响成本比较的敏感性。
4. 配对网格补充缺主入口/明确状态：8独立材料、2条件、64任务的数据已生成，逆解未派发，封存/评分未执行；不能仅写“未评分”或混成16材料。

没有在已完成小表中发现16对象成功数、CPU均值/比例、Global增量身份、吸收型8DEV负结果或3DDEV12指标的数值不符。Global已覆盖G19相对CIE的两个壳成功，G19相对Global仅多新环，无成功回退。吸收型K19/K0各1/4且均仅DEV；3D12行确属2曝光对象×6臂，全部320 accepted cap，没有将其作为fresh50结果。理论对有限rank/full derivative、噪声球不等于统计不可区分、face Hessian浮点证据与严格证书、current界不等于材料恢复均保留了限制；未发现需另列的实质过度主张。

检查6份指定公开文档的47个Markdown/HTML链接或图片目标，本地缺失0；GitHub/raw URL仅映射到本地同路径检查，未联网、未宣称远端已发布或HTTP可匿名访问。配对网格缺的是入口/状态披露，而不是其现有协议文件不存在。证据在[LINK_EVIDENCE.json](LINK_EVIDENCE.json)。

范围：RESEARCH_REPORT_ZH.md、THEORY_DERIVATIONS.md、PAPER_DRAFT.md、README.md、index.html、CHAT_HANDOFF.md；核对已授权16对象80行小评分表/组合表、吸收型8DEV摘要、3DDEV12行材料表与公开报告，以及中断和新网格组公开回执。未读任何未解盲GT/大数组/trace、未运行solver/scorer、未重复独立48代码审计、未改公开文档或活跃代码。根线程决定是否修改；本审核不延伸新研究。
