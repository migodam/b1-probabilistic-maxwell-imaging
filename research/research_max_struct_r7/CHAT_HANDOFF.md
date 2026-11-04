# MAX-STRUCT-R7 给Chat入口（2026-11-04）

先读[中文完整报告](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/research_max_struct_r7/RESEARCH_REPORT_ZH.md)、[网页与所有图](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_max_struct_r7/)、[推导](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/research_max_struct_r7/THEORY_DERIVATIONS.md)、[论文初稿](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/research_max_struct_r7/PAPER_DRAFT.md)。

## 已实际执行

原30个新材料×6方法=180；后选家族20新材料×4=80，分母分开。另有初始化18、旧九数字9、公开开发20；三个旧失败锚点15个校正端点（含控制）；真正矢量3D初始16+长校正16；真实Fresnel2对象×3=6。没有NN/GT在线选根，所有新GT评分在预测全封存后。R1–R6未改写。

|原30材料成功|CIEsingle|Dual|Global|A|A+|direct物理|
|---|---:|---:|---:|---:|---:|---:|
|成功/30|21|20|20|21|21|10|
|CPU秒|1543.60|1539.15|1512.97|2495.14|1711.34|838.18|

A+是RA=β(I+DzGoff)(J−DzH)的条件Möbius状态二次驻点+原prior-MM接受；不是精确正则化条件最小、不换prior，保留完整loss fallback。普通global-current不改。相比A，原30成功集相同，CPU少31.4%；相比成熟Dual仍贵11.2%。A-global平均NRMSE改善CI负，但成功率差CI跨0，不能说普遍胜出。

原30两个大阳性obj12曲条/20壳体；obj20 single也成功，obj12 single失败。obj09是IoU边界loss（.6966vsGlobal.7126），obj17是相对single边界。等成本按公开Phi incumbent检查点，不按GT选迭代；较低Phi可以丢失材料IoU。保留预算分辨率/float32限制。

## 最后的家族确认

在原30后选家族、再生成20独立材料，所有算法固定；不要和原30合并成事前50对象总体胜率。A+/Dual/global/single为12/9/9/8；曲条9/8/8/5，壳体3/1/1/3。壳体A+和single成功集合相同；曲条target03仅A+成功。A+对Dual平均NRMSE差−.06831[-.14213,-.01242]，成功率差CI[0,.30]。总CPU1089.69/1061.57/1027.61/855.19；对single等CPU mean误差CI跨0。八个A+失败完整保留，均有更好离线GT-Phi见证。

## 其他路线和资格

- 精确rank2/8/32 finite response端点正确，未获得优于同corrector原anchor续跑的独有材料恢复；15端点0成功，只涉及3个旧锚点。低秩动作不够不是公式错误。标准Woodbury/Gram/MM/Möbius/GENOME不是新发现。
- projected8丢失b20d，说明内块下降不保外路径；最终A+保留两旧bar材料恢复。九数字A+8/9，旧Dual9/9历史引用，不重新标盲测。
- 真3D矢量DDA375/648维、125/216材料DOF已经实际反演；共同长corrector后方法接近，无持续A优势。binary Mie refinement失败，一次体积分数修正后误差仍约8%；不能称连续3D资格/成像优势。3D用独立prior和相对L2 NRMSE，不是二维A+测试。
- 官方Fresnel原始文件已取得，六次真实反演+留出预测完成；三臂同终点、触box、天线校准不足，无可信像素GT或A优势。原数据/全文不再分发，公开获取脚本及derived结果。
- NIE/LUCIE未独立成熟复现；MLP/LNO/LSL未运行。不把未做算法伪装baseline，也不证明学习普遍无用。

最适合写A-material低成本条件读出与适用边界论文，不能称全局极点学习框架或泛3D/实测胜利。下一步优先独立物理/实测照明校准，保留有条件曲条收益并解释残余错误结构；别加无依据gate、扫β或用大网络掩盖。

最终配置configs/aplus_final.json；完整复现REPRODUCE.md；费用COST_LEDGER.md、全部否证FAILURES.md、内部审查audit/。大型原始轨迹在[本轮Release](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-max-struct-r7-20261104)，所有分卷解压到仓库根。原30和后20的GT现在公开，此后都只是开发数据。
