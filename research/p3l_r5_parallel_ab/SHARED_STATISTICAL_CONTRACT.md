# 原 P1 统计与物理合同

使用原源码的66维 Gaussianized latent、固定噪声 sigma=0.00119211798514519、原饱和decoder、24×24 TM/VIE、8×32复观测。每源依次实部和虚部，实白化系数sqrt(2)/sigma，Phi=||F-y||²/sigma²。材料统计为原latent law的推送，不增加材料Jacobian体积作为prior。

A研究明确归一化的二维截断子模型，其他64维固定，不是66维边缘。B输出监督材料风险动作，不因Gaussian context积分正确就宣称真实后验正确。noise-ball仅用于单独物理代表决策。

完整预注册、独立性、门槛、费用与失败处理见[FROZEN_PROTOCOL.md](FROZEN_PROTOCOL.md)及shared/cost_ledger/AUTHORIZATION.json。
