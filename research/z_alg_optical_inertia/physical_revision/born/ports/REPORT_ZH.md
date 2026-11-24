# 固定端口消光判据的弱反差和暗室复核

固定的全环端口权重在 96 份强散射读数上给出的 Born 距离下界为实测范数的 **40.23%–61.56%**，全部高于各自 1% 或 2% 的噪声半径。这 96 份读数来自六个物体，在不同的发收数和噪声水平配对重复。只看原生十发二十收，12 份读数的下界为 **49.09%–61.56%**。同一十发二十收权重对六份弱散射对照均给零下界，其中包括较灵敏的逐观测 Born 锥搜索已经排除的低反差双圆盘。零下界使这六份读数在固定判据下仍待判定。 [全阵列消光结果](../extinction.csv) · [配对阵列源谱](../../array/source_spectra.pdf) · [逐观测弱散射对照](../controls/RESULTS.md)

## 弱散射对照

| 生成模型 | 固定端口带符号读数 / 实测场范数，1% / 2% | 固定端口距离下界 / 噪声半径，1% / 2% | 独立结果 |
|---|---:|---:|---|
| 超弱椭圆，全波 | −7.23% / −7.26% | 0 / 0 | 两档球内各有显式 Born 拟合 |
| 双圆盘，峰值 0.04，全波 | −6.84% / −6.85% | 0 / 0 | 1% 与 2% 的逐观测 Born 锥距离下界分别为 \(5.348r\) 与 \(2.748r\) |
| 同几何双圆盘，Born 生成 | −9.39% / −9.36% | 0 / 0 | 两档球内各有显式 Born 拟合 |

弱散射结果沿用同一个原生 \(S,B\) 和同一个端口权重，连权重相位也未重选。带符号值为 \(\operatorname{Re}\langle Z,Y\rangle_F/\|Z\|_F\)，再取正部给距离下界。原 96² 算子与阵列扩展中的十发二十收算子逐元素一致。 [核对结果](results.json)

## 暗室复核

按每个源的实际接收位置建立权重。对每组 \(g\)，以 1% 相对奇异值截断从 \(S_g\) 和空场校准的 \(B_g\) 得到权 \(C_g\)。在全部源组上先逐像素加总 \(Q_i=\sum_g\sum_t(S_g^*C_g)_{it}\overline{(B_g)_{it}}\)，再按原规则选择共同系数 \(\eta=\max(0,\max_i-\operatorname{Im}Q_i/\operatorname{Re}Q_i)+2\times10^{-8}\)，令 \(Z_g=(i-\eta)C_g\)。这一步只用几何和空场参考。区间核验每个像素的合并实内积非正，因而对所有共享的非负实材料成立。在线读数是各组 \(\operatorname{Re}\langle Z_g,Y_g\rangle_F\) 求和并除以共同的权重范数。缺测端口不参与，也不补零来推算共同照明。

| 数据 | 保留端口和频率 | 距离下界 / 所用源接收场范数 | 高于两种经验标尺 |
|---|---|---:|---:|
| Manitoba | 奇发偶收、偶发奇收，各 12×12；3–6 GHz 七频 | 28.16%–88.30% | 25/28 |
| Fresnel | 每个源自己的 241 个实收点；4 GHz | 38.36%–44.85% | 3/3 |
| Fresnel | 同样的原实收掩码；6 GHz | 43.58%–48.64% | 0/3 |
| Fresnel | 同样的原实收掩码；8 GHz | 39.34%–49.93% | 0/3 |

全部受测组的共同相位都通过了逐像素符号核验。Manitoba 使用旧报告的两套金属参考差、空场漂移加留出金属残差。两种标尺需同时对照时，唯一未过的三行是 3 GHz 小尼龙柱、4.5 GHz 大尼龙柱、5.5 GHz 大尼龙柱，下界分别为 34.06%、28.16%、59.63%。Fresnel 复用逐源空场拟合与旧报告的留出空场均方根迁移值，也复算加入诊断复增益的放大尺度 \((\delta+\varepsilon)/(1-\varepsilon)\)。其 6 GHz 三行均超过简单空场标尺，但低于含复增益的标尺。8 GHz 只有外嵌介质体超过简单标尺。 [Manitoba 原校准曲线](../../../validation/manitoba/calibration_margins.csv) · [Fresnel 原标尺和曲线](../../../results/fresnel/SENSITIVITY.csv)

这些读数检验的是选定校准和离散模型内的无损非负 Born 类。经验标尺由原校准对照取得，尚无覆盖真实测量的统一校准和三维模型误差上界。Manitoba 尼龙另有小损耗，实测偏差可能包含线性有损响应。因而上述计数是对误差预算的条件比较，不能单凭它们判定暗室对象必须发生可分辨的多次散射。

## Reproduction

Run the following from the project root with NumPy, SciPy, python-flint and threadpoolctl. The replay reads the original Manitoba private calibrated port archive and original Fresnel raw and derived records. The measured replay reads local caches built from the collectors’ data. Its weights depend on the empty-field calibration and stay with those local files. The public archive includes the derived tables and reconstruction entry points. [Data access and local preparation](../../REPRODUCE.md#暗室数据的本地复算)

~~~sh
python -B research/z_alg_optical_inertia/physical_revision/born/ports/replay.py
~~~

The single-threaded replay took 22.39 s wall and 22.00 s CPU. The original synthetic geometry and six weak readouts took less than 0.2 s. All Manitoba geometries and readouts took 0.39 s wall. Fresnel took 21.83 s wall including reading the original files, repeating the empty-field calibration, forming nine masked weights and checking the pixel signs over 3,600 cells. [Replay](replay.py) · [Interval results](results.json) · [Original extinction weighting](../extinction.py) · [Nonlinear source curve](../../../submission/curve/noise_curve.csv)
