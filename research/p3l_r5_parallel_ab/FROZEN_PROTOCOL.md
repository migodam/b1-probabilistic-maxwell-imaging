# R5 冻结研究合同 v1

冻结日期：2026-10-26。基线提交：85b7bcbf39f0eeab7062ecce17531538151082bf。独立分支 research/p3l-r5-parallel-ab-20261026。当前用户任务是执行和预算授权的依据，Pro ZIP 是数学与设计输入。

## 共同目标

原 P1 Gaussianized X~N(0,I66)，原 quintic C² 饱和 decoder，24×24 TM/VIE，256 ROI，8源×32复接收器，sigma=0.00119211798514519。每个 source 按32实部、32虚部打包，f=sqrt(2)pack(F)，Phi=||F-y||²/sigma²。无额外材料 Gram prior，不将 noise-ball 当后验硬支持。

历史只读。本轮真实求解通过 shared/physics_adapter.py 和追加式账本；旧成本另列。每个 endpoint 立即保存 material/latent/prediction，指定 reference 才另存 current/field/LU。采用无损压缩和有界缓存；单物理 worker、BLAS单线程。

## Stage A

三个事先固定二维模型，活动坐标为(0,1)、(0,64)、(12,27)，其余坐标公开固定为0。它们是新的二维子模型，不是 full66 marginal。采用归一化截断标准 Gaussian [-4,4]²；截断是模型定义，域内遗漏和数值覆盖仍须评估。每模型两个独立 parent，三档噪声使用共享标准 innovation，并分别生成 d_sigma=f(x*)+sigma epsilon，比例1、1/2、1/4。只在实际发现时声明饱和折点。

解析测试通过后才运行 Maxwell。原坐标积分与响应坐标积分分别核对 evidence、矩、弱边缘、分支和覆盖。响应目标保留 prior、determinant、normal residual、sigma^rank、branch及重叠权重。未建立全域严格包络时，数值参考和经验曲率不称认证。

门：TV<=.05；可计算双向KL<=.1；branch误差<=.02；遗漏<=.01；主要材料均值误差<=posterior SD的1%，退化任务绝对floor=1e-8；参考误差小于比较效应1/5。扩大另需同等可信精度下相对最强基线成本降低20%，或基线达不到而新法可验证的精度。未通过不扩大Full66。费用和失败按 AUTHORIZATION.json 硬上限计入。

## Stage B

复用 joint TRAIN32768/VAL2048，昂贵context TRAIN1024/VAL256。最多8表示臂×3训练seeds：2026102601/02/03。旧R4测试GT已暴露，只作历史证据，不参与本轮选择或冒充新盲测。

核心比较：K32+冻结head；确定性矩+冻结head；确定性矩+匹配重训head；无复杂density直接矩/RRR。恢复原x坐标并检查协方差后，才能使用独立特例；一般相关Gaussian保留属性与场的相关。四通道严格匹配R4实/虚均值和方差。原K32 ddof=0与population variance差别显式记录，不假定去随机性必然改善非线性head。

RRR rank=2/4/8/16，ridge=.001/.01/.1/1，仅用TRAIN交叉拟合与VAL选择。GT残差PCA不称可预测风险谱。真实Maxwell特征与raw-y、decoder-only、same-spectrum随机控制比较，不能同时配平时报告Pareto。

开发门1：保留Gaussian fusion相对strongCNN新增收益至少80%，且额外provider推断成本不超过原1/3。仅删K32、保留encoder不能称删除provider。开发门2：真实Maxwell比最强控制相对MSE改善至少2%，独立配对区间支持正收益。

开发门通过且原强anchor合同能在32768确认RHS内执行，才运行两批独立256-parent确认；否则标RESOURCE_LIMITED或NOT_RUN，不降级anchor却声称同一基线。先冻结最多3个finalists、checkpoint与推断参数，再完成全部在线预测，最后由隔离evaluator读取GT。单位为parent，seed平均loss与图像ensemble分开。4000次paired bootstrap，三个主要比较采用98.333%区间；其余探索性比较单列。

## 修订与证据

实质修订保留原协议和失败，并写ITERATION_LOG。不改变noise/prior或评价支持救结果。A/B互不作为停止门；独立结果成立且剩余授权充分才考虑融合。概率样本、风险动作、合法代表分开。物理费用权威为 shared/cost_ledger/LEDGER.jsonl；阶段训练费用分开，根账本只做明确去重汇总。
