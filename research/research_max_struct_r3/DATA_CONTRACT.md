# R3物理、数据与评价合同

## 共同二维开发

读取同仓库R2公开 `data_public/obs_b20d5a55677d4eca.npz`；不是新生成。完整24×24 disk-pulse scalar-TM，576个共享实材料，12外部源、48外部复接收，k=2π，exp(−iωt)，背景1、μ固定1。完整resolvent求解，无Born截断。box [.02,20]，host初始化.1；本轮续跑/有限路径使用已保存CIE端点，逐run指明。

Y含CN噪声：σ是复幅度标准差，实虚各σ²/2，576复数据对应99%半径25.165422。最终白化范数||F−Y||/σ；原在线目标仍归一化数据+.0001 smooth/.0001 smoothedTV（差分19.5，epsilon.01）。等价白化目标需同倍缩放prior，不换先验。

NRMSE=||m−mGT||/(12√576)；支持评分为预测自身256-bin Otsu，GT mask仅评分时读取。继承材料.15、IoU.7及噪声门，但报告实际值，不以gate代替分析。Otsu near-tie浮点bin可有小差别，冻结实现不改；所有开发联合结论不受该差别影响。

所有online候选由原物理Phi排序。GT只用于已经公开的固定材料诊断和预测封存后评分。clean_native读取已公开Y_clean，是明确的信息反事实，不能与含噪算法收益混账；没有clean Gaussian门。

## 各类证据分离

- 11个开发控制是同一独立材料对象，不计为11个独立样本或新盲测。
- 单圆盘控制是公开1维真实Hankel有限模型；clean/no prior，mGT1.2，边界初值20，box [.02,20]。它不是576参数成像。
- Z-ALG fiber例使用原6-cell/2-source/1-receiver geometry与旧公开材料对；新增uniform i损耗只用于资格诊断，不给反演新增观测。原实部[.5,4]，被动虚部0/1。
- 三维例是真实矢量CM+RR DDA：8 cells，三源横向，6外部xyz near receiver；实SPD锚点，一cell的vvᵀ/I更新t1/2/4。无noise、无拟合、无GT恢复、无连续精度主张。
- noise-floor分布要求固定材料/模式/状态权重且源噪声独立；未知test材料不按GT计算online门。数据选择后的统一界未证明。

`data_public/`只保存依赖说明；既有R2公开NPZ引用不复制算子制造新dataset身份。`sealed_truth/`本轮没有新增隐藏truth；旧truth已公开不能重新标blind。准备的生成器和DRAFT协议从未执行。
