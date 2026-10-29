# 强普通Full-J基线与费用审计

## 物理与统计

24×24 scalar TM/VIE、576 full-grid、256 ROI、8源×32复接收、原被动材料与C² quintic decoder、Gaussianized66 prior、sigma=0.00119211798514519。复数据source-major。白化实残差为$\sqrt2\operatorname{pack}(F-Y)/\sigma$，$\Phi=\|r\|^2/2$。

R7材料半径与历史竞争门直接读取冻结源码。E2 Arm A在全部16个开发parent-radius对上复现已保存的R7 Full-J端点Phi。该结果是复用基线一致性，不是新科学发现。

## 全部J的成本

完整forward为8 RHS，接收端adjoint构造完整latent Jacobian额外32 RHS。若当前corrector已付费生成$J$，动态weak SVD不再收一次40 RHS；其CPU分解时间单独记录。固定weak frame因此没有自动省掉一个物理Jacobian。

## exact ROI强后端

新E5全部方法共用精确256格ROI后端，不能只给固定方法更快solver。forward与adjoint费用按256维计，setup的264个576维RHS只支付一次并显式摊销。原48次匹配forward+J kernel计时包括八对象两个合法状态、三次旋转重复：direct共2.145907秒；ROI在线0.870254秒；setup0.016972秒。setup-inclusive kernel节时58.65%。

这不是完整成像系统节时58.65%，也不是证明每种新采集都同样快。E5另报anchor、corrector、分解、decoder、保存、失败和全作业费用。不能把不同系统规模RHS直接相加成等效物理工作量。

## E5普通控制

RANDOM66_FULL与B/C共享同一anchor、材料半径、最多forward/J数、signed progress和原decoder。它保留更多corrector自由度，属于较强的普通物理控制。随机提案、固定提案和动态提案都必须付拒绝与失败费用。首个GT可用时点在全部32对象预测封存以后。

E5采用R7 direct Probit64 anchor，不是R5-B Native96+Probit64；不能把两者端到端成本或材料风险当同配置直接比较。
