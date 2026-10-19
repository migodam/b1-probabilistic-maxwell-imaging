# 概率状态：没有P1可信完整reference

完整目标为p0(z64,u2) exp(-Phi)，Gaussian×uniform原prior，gain固定，sigma原值。没有变换成材料Lebesgue均匀prior。

C温度链的pCN+反射property proposal为prior可逆，使用完整likelihood MH；交换目标也正确。代数检查最大prior抵消误差2.49e-14。只有16 sweeps，mixing不足；不报告posterior mass、KL、CDF credible coverage或calibration。

Flow/Gaussian raw delta-z log probability通过显式Gaussian/affine-coupling密度定义；这只是提案密度。material retraction、固定距离和deterministic corrector后的频率不是posterior。概率全部null，不用residual softmax。

Package剩余已建模mass、未探索mass与模型外mass均无法量化，不把展示K项重新归一化到100%。已有Joint P2 Mie reference属于另一target，未用于这里的P1概率背书。
