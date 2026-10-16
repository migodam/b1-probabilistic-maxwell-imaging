# 复现入口与封账保护

实际执行源码：b87413bbe06100600988ed87f69b48e21dfe3230（F1.5及首轮）；4d45c3b（一次舍入修复后的恢复）；bd60ee6（新增9中点诊断）。完整ID可由仓库git rev-parse解析。后续只增加文档/绘图、legacy路径和背景路径fallback，不再调用物理。首轮源码副本在evidence/source_executed_attempt1。

环境复用Gaussian/.venv_nn/bin/python，NumPy/SciPy/Matplotlib已有，不安装GPU或新服务。后台原A30 source717d3d26a290248bd2b9aeb55f3d03694b407b84与A26生成器2603f3bd16c48e48587fedc4363975779def17f5；full物理接口不另写求解器。已有B1源码包含prior rasterizer快照fallback。

零新增物理：tools/research_f3_prior_audit.py可重做商先验/decoder/缓存矩阵审计；tools/research_f3_results.py只读结果做统计、离线GT与绘图。请指向新分析目录或保留本轮产物以免覆盖记录。前者默认旧缓存路径，后者需要Release原始数组；未读取或应用LU。

付费代码是python -m research_f3.run；恢复及中点代码分别有--recover和research_f3.path_audit。正式目录STAGE_CLOSED及同一本ledger阻止误重跑，不要删除/改名它们绕过额度。新的物理复现必须获得独立授权并使用新明确目录与预算；代码、固定协议和原始结果已完整交付。

Release包含55个新endpoint的完整材料/current/field/prediction/LU/pivots、8个新Jacobian、新方向Hessian/补偿数组、MH保留轨迹及所有proposal、两份旧F1付费state/J、原已支付背景G/S/incident、在线观测和离线GT。

若本地没有原F1 worktree，B1_F3_LEGACY_ROOT可指向本Release的inputs/legacy_cache（保持results/F1与inputs/online布局）；模型优先读取本轮inputs/background。该路径fallback是发布可移植性修正，未改变实际运行数据或物理结果。

数学对象与费用：closed weight支持、canonical roundoff clip、stable logprior、原uniform shape、7q密度2L*p_rho、fixed diagonal reflected MH。原noise/gain、material amplitudes、knownconditional globalgeometry完整在configs/PROTOCOL。

不做SHA256核对；没有第三方全文或私人对话。旧账本、Gate与冻结A17稿未修改。
