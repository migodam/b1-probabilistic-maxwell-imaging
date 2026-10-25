# 概率审计：没有用风险改善代替posterior正确性

完整目标是原P1 latent posterior，不是局部切片或profile。模型明确区分完整可归一化Spline/高斯/三角密度与CFM数值生成law。对风险图、合法代表、几何candidate，不输出伪posterior权重。

## 固定生成样本的必要物理检验

若X'是真实posterior draw，Y来自原prior predictive，则Phi(X';Y)应服从Gamma(256,1)。每parent固定index0，不用best residual筛选。实际两批分别为：

|模型|512平均Phi|256平均Phi|结论|
|---|---:|---:|---|
|raw Spline seed1|8,371.95|8,662.42|严重失败|
|anchor普通CFM|36,282.18|40,924.89|严重失败|
|anchor几何CFM|22,139.70|23,557.96|比随机改善，绝对仍失败|
|anchor随机控制CFM|47,633.59|49,493.64|严重失败|

原期望为256。物理PPC共7,104个实际状态、56,832RHS、无失败，且未根据GT或response挑样本。主Spline的PIT超过.975为510/512、255/256。生成law没有经过FULL资格。

边际coverage偶尔接近标称不抵消此失败；plateau原子与ties使用随机化rank单列，不能将全图含大量宿主单元的覆盖率当ROI或完整joint校准。全部latent/material/projection/phase分层指标见evidence_summary/DENSITY_*.json。

K64重要性采样仅固定前16parent，ESS中位1、最大2.273。它是权重集中度，不是独立链ESS，不构成posterior teacher。短校正链保留所有拒绝重复状态，形式MH比正确，但material函数ESS/MCSE仍UNQUALIFIED。

## 仍然成立的相对概率证据

在同一joint law、原latent参考测度和独立parent上，显式模型proper NLL差等于未知posterior的期望forward-KL差；未知conditional entropy相消。由此可判Spline比同容量Gaussian更接近，不能说已知绝对KL或完全覆盖。CFM只用样本energy score，不给不稳定数值密度加资格。

数据依赖rank在固定前64有实际truth prediction的对象上进行；没有缓存truth prediction的其余对象标NOT_RUN，不新增物理补齐后伪称全768 rank审计。Gamma固定draw本身覆盖全部768。

## 输出身份

风险估计器：监督material-risk action；生成样本：qθ的材料推送；约束代表：candidate-only合法世界；几何端点：可达性见证。四者分开保存。posterior probability/mode mass默认null，禁止residual softmax。


## 封账后的数学解释补充（零新增物理）

固定样本的PIT尾事件还可通过标准数据处理不等式，给出平均条件材料TV/KL的保守下界，不需要MCMC teacher。完整26模型×2批校正后，rawSpline seed1的两批TV下界为0.953/0.939；这是原模型与理想独立随机流前提下的post-hoc解释，不是每观测保证、预注册主结果或新Maxwell定理。详见[数学证明与全部限定](THEORY_CALIBRATION_ERROR_BOUND.md)。
