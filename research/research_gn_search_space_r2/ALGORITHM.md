# 算法与边界

输入是冻结R4的固定G/E/S、公开复观测Y、合法box和原Native/Reciprocal材料端点。所有材料均为576维共享实变量；不输入数字模板、GT方向或额外照明。

## 空间分析（E1，离线机制入口）

1. 从上轮缓存读取固定Green谱，强门σ/σmax≥0.1；复用固定接收核，不重算其SVD。
2. 构造暗电流空间正交补U=orth([Vin,strong,Vob,nonzero])。
3. 当前材料做一次完整VIE因子分解，求总场Q、完整材料Jacobian和L的逆动作W。全部照明共享LU；显式W付费576RHS。
4. 实化U*W diag(h)Q，直接SVD判定严格数值核；不对平方Gram小特征值作零核认证。
5. 使用两个实材料Gram的广义谱求泄漏≤1%、≤10%空间。Euclidean正交化后计算材料投影。
6. 对1/3/6/12源分别报告当前J、Born、公共背景J和离线GT J的谱。GT仅用于后诊断，不构成可部署算子。
7. 所有基按物理谱排序，再计算GT误差覆盖。完整结果同时保留256维之后的陡增，避免稀疏采样图误导。

主要入口：`src/space.py:prepare_geometry`、`experiments/run_geometry.py`。返回J谱、泄漏广义谱、Euclidean近似基与费用。

## 公共数据搜索（E2/E3）

```text
输入 observation, anchor, common prior, budget
开始计时；加载原物理Model；anchor加入在线候选
按方法准备当前J或当前Green材料拉回，并付费
对两级Haar细节、两符号、两RMS尺度：
    用同一原始提案p
    Haar: h=p
    Green: h=Papprox p
    Soft: h=V diag(lambda/(s²+lambda)) Vᵀp
    沿整条射线缩步至box内，不逐像素裁剪，不重新放大
    原完整Maxwell评估seed；按共同目标更新incumbent
    全576维TRF校正，最多3秒/8次函数评估
    保存before/after/拒绝与所有成本
剩余预算继续校正在线incumbent
封存预测与原始物理目标排序
独立评分器读取已披露GT，计算top1和best@K
```

E2直接从同一anchor继续成熟SciPy有界TRF；smooth/TV同时乘1、0.1、0。E3全部使用原共同权重，两种过滤不获得额外正则或更有利初值。全部臂45秒、500LU、30000复RHS；正常提前收敛即结束，不反复空转耗满预算。重复E2/E3原权重clean对照只计一次物理计算。

Haar指当前材料图中的两级细节方向：图像减2×2均值，以及2×2均值减4×4均值。不是遍历所有Haar基原子，也不是从GT提取方向。因此本轮失败不能外推为所有wavelet proposal都无效。

## 排序与物理资格

在线分数为0.5||F(m)−Y||²/||Y||²加原smooth/TV。原anchor永远保留。所有seed、corrector迭代均实际物理评估；保存的oracle子集为anchor、top1、before-seeds、per-seed校正最优，至多18槽，不是18个独立解。

历史成功门：clean数据残差≤1e-3，noisy≤.003；联合还要求NRMSE≤.15、Otsu IoU≥.7。E1的clean1e-4只是相对anchor的数据变化诊断容差，与历史成功门不同。

## 费用与可复用范围

固定G/S谱属于公开物理几何的offline准备；本轮不重新收费一个未运行的分解，也不把其冷启动时间说成零。每个搜索臂重新计算当前J/拉回，不能免费调用E1缓存。完整J利用源与接收RHS、共享LU；没有576次逐像素差分。

 dense材料J谱、576维广义本征与W准备都是当前实现费用，未声称最优实现。候选corrector用成熟TRF；本轮不是新的GN定理。

## 失败模式

- 电流暗核的材料拉回可为零，不能临时放宽门制造根。
- 相对电流泄漏小不保证绝对数据漂移小。
- box活动点使可行集成为锥；正负射线可能阻断。
- 弱空间只提供不敏感性，不能提供GT方向或正确符号。
- 公共正则可能使在线分数更好的材料反而远离GT。
- 原gtol对微小归一化目标可能提前停止；严格停止诊断单列，不替换主实验。
- 原生数值worker失败有未恢复的费用；失败保存，不宣称完整加速。

本轮没有训练NN，没有新概率解集，也没有验证全局材料解纤维或TSOM截断因果关系。
