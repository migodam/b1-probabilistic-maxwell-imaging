# G2：本轮没有启动

状态：NOT_RUN_PREREQUISITE / NOT_RUN_RESOURCE。没有新30/60-parent locked确认，没有新A0优化，没有reopen速度对比，没有新acquisition。

旧A0三法0/24、5/24、9/24配对成功率及补优化均未达到一般可靠在线锚点；这些结果保留，不作为充分收敛的普通FWI上限。本轮没有将oracle anchor换入在线算法。

后续若G0/G1成立，强基线需用普通每源加权伴随、合理连续收敛与可逆支撑搜索；同样给予背景LU、Schur、goal-ROM、snapshot和cache权限。达到call ceiling就是NONCONVERGED，不能画成MAP。2×2有限离散prior不表达孔洞，3×3孔洞合同必须独立冻结。

当前静态全部y的representation update不增加因子。KEEP_ACTIVE/KEEP_MULTIPLE/ASSIMILATE/REOPEN必须保留原参数、prior条件核、joint权重和历史版本；task超合同、误差界失效或stability变差则恢复原分支并替换旧likelihood。新增独立静态通道必须另作面板，不归因成BGGA调度收益。

任何未来速度比较都应给真实cold `setup+query+certificate+reopen+fallback`，并和普通cache/ROM比较。cached decoder/小核心时间、较少reopen或少参数不构成加速证据。目前没有被实测安全避免的FE或RHS，COLD_COST_NOT_ESTABLISHED。
