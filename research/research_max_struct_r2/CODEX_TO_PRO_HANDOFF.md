# R2 交接：不要把负面证据改写为新算法胜利

先读RESEARCH_REPORT_ZH、THEORY_DERIVATIONS、独立审查和成本。R1未变；R2用户明确授权GitHub发布。

已完成3新独立对象×4方法盲测，全部12预测封存后评分。成熟native/global都2/3，profile1/3，direct0/3；新环和双岛恢复，弯条都失败。global弯条NRMSE下降而residual变坏，不能称同数据新增成功。profile损失了双岛。三维另有一份实际新增模拟观测，在已知一维张量族内得到相容估计；不是自由张量成像、不是CIE优势。

最后一次事后开发control将失败native的最终mf10改mf12，1500步、β、prior、corrector相同。data residual只从.100735到.097693，仍失败。保留它单独的预算与POST_HOC标签，不改3对象盲测统计。没有继续调β、先验、初值或噪声。

full-profile的可变残差度量有vanLeeuwen–Yang2024直接先例；精确状态消元、envelope、affine差有更早先例。globalPR也是既有机制。可用成果是实现、失效分解、信息补偿实例和反例；PAPER_HOLD。不能以理论成立/更低auxiliary声称成像优势。

未解决：噪声DPCS硬锁冲突；corrector驻点/预算边界；一般被动张量连续同一分支；M的兼容Hcurl；T的统一k和delta；候选族覆盖/连续分离认证。下一步只选一个可判定问题：在同一原似然下修正DPCS强模态，并和已知CC-CSI/普通软释放做同成本控制。Atlas已有这一方向；需要新的量化机制，不能只换算法名。

新鲜确认必须重新生成对象，本轮truth现已公开。禁止测试GT选根。不要重跑已有恒等式到机器精度，不为更多图把有限预算掩成收敛证明。所有历史和新增费用分开；代码和原始结果可从本目录及同仓库依赖复现。
