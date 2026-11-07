# BothBlocks16独立确认：事实汇总

冻结后的16个新对象、五臂共80次反演全部执行并封存，0失败、0重试。80个单臂预测、64个预声明组合选择、80个共同物理LU选择，均在GT评分前封存；主线程核对后明确授权离线评分。算法、生成seed、原prior、初始化、双beta、预算、候选及在线Phi规则均未按新结果调整。

主比较不支持BothBlocks优于材料-only G19：BothBlocks材料成功5/16，G19为7/16；均值amplitude NRMSE分别.173574与.144847。paired均值差(Both−G19)=+.028728，high/mixed分层、object为单位的10000次bootstrap95%区间[+.000176,+.084961]；NRMSE 5胜11负。材料门0改善、2回退；材料成功率差的区间包含0。该有限样本结果与曝光DEV桥的成功必须分别报告。

|方法|材料PASS|CN99 PASS|mean amplitude NRMSE|median amplitude NRMSE|实际完整child CPU秒|
|---|---:|---:|---:|---:|---:|
|cie_single|4/16|12/16|0.231343|0.112186|767.393173|
|cie_dual|5/16|13/16|0.210683|0.108894|777.117487|
|global_dual|6/16|13/16|0.213442|0.108553|758.285589|
|fullG19|7/16|14/16|0.144847|0.103577|973.533201|
|fullGreenBothBlocks|5/16|13/16|0.173574|0.108908|1005.468849|

BothBlocks在high组0/8材料通过、mixed组5/8；G19相应为2/8和5/8。BothBlocks对CIE有11/16 NRMSE胜、均值差−.057769、材料门多1例；这不能替代预声明的BothBlocks对G19机制比较。两个G19→Both门回退须区分：

- r11_bothblocks16_0000 annulus: G19 NRMSE/IoU=0.100323/0.723684；Both=0.549488/0.296774。
- r11_bothblocks16_0002 slotted_shell: G19 NRMSE/IoU=0.106830/0.702247；Both=0.108897/0.664894。

前者是大材料误差差异；后者主要是固定IoU在.7门两侧的变化，不能把两者都描述成同等大的材料恢复崩溃。CN99是数据一致性门，不能替代材料恢复。

预声明付费部署组合以originalPhi选材料、完全相同Phi时保留CIE，所有选择先于GT。CIE+Dual与CIE+Both均5/16材料通过；其mean NRMSE为.210005与.211445，二者paired差+.001441的区间跨0。实际部署CPU逐对象求和后总计1544.510660与1772.862022秒。两组合共享campaign中的CIE执行，不能把重叠组合成本再次叠加为研究总成本。

所有单臂实际自然结束均早于240秒，因此此次240+240组合视图与完整端点视图相同；允许上限同为480，不等于实际花费相同。单CIE480只有一次实际child花费；两起点组合始终支付两次。共享两lane campaign wall不当作一次部署wall。

共同physics表是新增的已预封存次级视图：每对象在五臂已达到、CPU合格的物理LU中取共同cap，按cap内最小Phi材料评分；不是更改主CPU端点。五臂材料PASS仍4/5/6/7/5。Both mean NRMSE .173374，G19 .145319。只匹配物理LU不等于匹配总计算；全部辅助小LU与setup仍单列计费。

PAID_INCUMBENT_FRONTIER.csv包含158626个真实付费Phi-incumbent节点的CPU/LU/RHS与材料分数；ACTUAL_TERMINAL_SPEND_POINTS另列80个完整child实际花费点。曲线连接仅为阅读引导，不进行状态插值。

确切wait4计账：80反演child总CPU 4281.798299秒；生成child 1.372750秒；controller 0.448456秒。共享campaign wall 2182.972411秒，不与CPU相加。生成物理cost与反演cost分开；postseal scorer/frontier/画图的外部CPU见POSTSEAL_EXTERNAL_PROCESS_COSTS.json。新16不与Windows48合并rawCPU优势。

本汇总只报告冻结合同和已测事实；理论有效性、唯一因果解释、创新性及后续方向由根线程判断。无新增优化或forward用于评分和绘图。完整16材料与CPU/LU图位于FIGURES/bothblocks_confirmation。
