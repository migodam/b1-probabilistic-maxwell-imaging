# 物理、电流和材料空间的区别

采用原P1离散合同，$D=\operatorname{diag}(m)$、$E=E^{inc}+Kj$、$j=DE$、$Y=Sj$。固定背景和采集固定的是$K,S,E^{inc}$。Born只取$Y_B=SD E^{inc}$；完整解为$Y=SD(I-KD)^{-1}E^{inc}$。后者不可因为$K$固定而把全部材料反馈替换成线性投影。

## SOM与投影式twofold

$S=U\Sigma V^*$给出固定receiver/current分析空间。TSOM式二次限制可写成$P_{weak}V_{internal}a$，即投影后的生成空间，不必构造集合交集。列不再正交时需保留Gram或QR；实际有效rank可能下降。详见公式F-OP及F-OPT系列。

receiver-dark只保证$S\delta j=0$。多源共享材料还要求$\delta j$与内部总场共同满足逐cell构成方程；把不同source电流分别扰动一般得不到同一个$m$。本轮选择的八个开发current扰动，材料lift相对残差约0.52—0.75，均记MATERIAL_LIFT_FAIL。这个反例只否定所实施lift，不否定所有SOM/TSOM。

## Fourier与Wavelet

完整换基若可逆，改变坐标不改变材料集合。截断、投影到原Gaussian64+properties支持则另有表达误差。本轮固定Fourier/Wavelet方向先通过真实材料Jacobian映回原latent，投影误差约0.22—0.70，不能以原始模式图的外观宣布材料可达。

## 固定背景的精确消元

令$P$选取256格ROI，$D_h$为固定host，$A_h=I-KD_h$，$B=A_h^{-1}KP$，$C=P^TB$，$E_R=P^TE_h$，$O=S(P+D_hB)$。对任意原支持内非均匀ROI增量$D_\delta$：

$$h=(I-D_\delta C)^{-1}D_\delta E_R,\qquad Y=Y_h+Oh.$$

这是完整非线性VIE的块消元，未改变材料prior、噪声或网格。一次host setup后，变化系统从576格缩至256格。所有ROI格同时变化时不是rank-one更新。均匀复幅度$D_\delta=vI$可以进一步利用scalar参数结构；任意非均匀对角矩阵不与$C$交换，不能套用同一种固定特征对角化。

该公式来自标准resolvent/Woodbury/Schur结构，本身不是新定理。真正要判断的是它是否组成有效成像链、覆盖是否改善、实际总费用如何。
