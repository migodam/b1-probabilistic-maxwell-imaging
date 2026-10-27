# GPT Pro handoff：当前真实完成边界

## 可以使用的新事实

固定公共target ROI的P3 mesh已实现，源项按散射场语义纠正，cut-polygon目标mass/解析incident积分已实现。实际mesh 28,903 DOF，公共2,257 ports；示例目标仅79活跃DOF。面积/质量/矩阵/有限概率单测实际通过，完整数量见CODE_CHECK.json。没有新PDE。

118个全部wall groups的完整公共Green eager准备需266,326 RHS，已超过G1建议上限约13.3倍。普通Schur复用同一paid cache，不再二次支付。不能用这个架构直接宣称经济BGGA。

## 严格不能说

- G0的独立field精度通过：未运行。
- 跨目标全波cache合格：未运行。
- 有2950原子新posterior或interval省查询结果：未运行。
- 新一般unknown-support anchor：未运行，旧A0失败身份保持。
- BGGA新调度比ordinary goal-ROM、Schur或joint smoother快：未证实。
- 资源不足导致科学假说为假：不成立。

## 下一步最小可执行工作

资源授权后先跑tools/g0_qualify.py，固定12条件、同源/噪声/PML、独立interface-fitted参考，第一条白化差>.1即停。数值误差必须小于相应task margin，.1不是posterior证书。完整四类各3个Maxwell stress仍需实现与运行，不能以几何脚本替代。

如G0足够精确，优先寻找真正低成本按需物理界或普通goal-ROM，之后再决定是否值得比较G1各方法。默认[0,1]界可能因高维likelihood太保守而求全表；应接受CERTIFICATE_COST_NO_GO，不用reference weights“帮算法”选择。

原始两个Git来源与缓存地图已保存；旧证据不重写。不自动扩G2，不训练网络，不新开锁定确认，不把SLAM动态模型带回这个静态合同。
