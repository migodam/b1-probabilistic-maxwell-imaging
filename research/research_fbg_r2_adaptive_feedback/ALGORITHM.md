# 可运行算法与公平比较

输入仅为公开 G、E、S、匿名复数 Y、材料 box [.02,20]、公共初始材料和预算。没有数字类别、支撑、字体、GT 或隐藏生成状态。主物理始终为 576 单元完整非线性 scalar-TM disk-pulse VIE；8→12→24 改变材料参数化，不改变物理矩阵维数。

## 1. 共同物理与选择

复用 `research_z_alg_cuda_pipeline/src/device.py` 的常驻 complex128 Green、共享 LU、多照明 RHS、解析伴随及 Full-J。四个独立 CUDA stream 在同一进程中运行；大型 TRF SVD 使用已有非阻塞 CPU/GPU 混合后端。

每个阶段支付其原始物理计算。在线 top1 按

`0.5*||F(m)-Y||²/||Y||² + R_common(m)`

选择。共同正则与 R4 相同：差分除 19.5，二次和平滑 TV 权重均 1e-4，TV epsilon=.01。辅助目标再小也不能替代这个分数。保存各阶段材料候选和完整目标轨迹；这些点通常不是独立 Maxwell 根。

## 2. 比较方法

- Native TRF：成熟 SciPy 有界 TRF，解析完整材料 Jacobian，`x_scale=jac`。
- Reciprocal TRF：粗材料系数的倒数参数，精确链式导数 −c²，相同材料空间。
- Data L-BFGS：原始数据伴随目标，材料宽度缩放，相同阶段及重启分配。
- GD-late / GD-all：原始 GlobalDual 分别仅在细层或各层使用；随后原数据校正。
- Fixed-Gram：固定公共 .1 参考场的 Gram，区分动态场权重与普通源加权。
- Scalar feedback：B=βI，固定参考尺度；Stage 1 包含请求 β=.05,.2,1。box 限定后实际缩放写入 receipt，不能只按方法名称推断实际 B。
- Current-fit diagonal / block：从公共空间方向的完整电流切向拟合，Cholesky 源归一化，ridge=平均行能量的 .05 倍；块为 2×2 像素（4×4 矩阵），批量小系统求解。
- Shuffled：循环平移对角系数 37 个格点，保持系数多重集和谱/范数；不是保持场相关性的反馈。
- Pilot-selected：有限候选各做一次有界辅助梯度试探，按完整原始物理目标选择。试探、失败和拟合全部收费。它是短视的经验筛选，不是 basin 证书。
- Q-beta：用 F09 的公共均匀参考梯度变号阈值自动选择标量 β。规则在 Stage 2 前冻结，不能使用该阶段 GT 选择 β。

## 3. 阶段式反馈

外层到达一个材料分辨率时，选定 B 并作整个 box 的可逆性限制；内层优化时冻结 B，因此不遗漏 DB 项。当前没有基于停滞的连续重拟合触发器，不能把阶段切换称作已实现的动态失效检测。

```text
prepare cached Maxwell operators and public reference field
for each anonymous observation, method, public start:
    for stage in 8,12,24:
        prolongate material coefficients
        choose/fill B at the anchor (or calibrate Q-beta once at .1)
        certify/restrict B over the original material box
        freeze B and reference normalization
        optimize feedback Gram using one original-Maxwell adjoint
        correct using original data objective
        retain best original-physics incumbent, including earlier stages
    return online top1 and paid candidate/trajectory archive
```

所有材料仍在原始 box 中。对角反馈做精确区间极值；一般块反馈一次缩放使 `20*max||Bb||₂<=.8`。原物理 Gram 失败会保留失败记录；Cholesky 成功不等同 Gram 条件数良好。

## 4. 多初值预算

Stage 2 的 8 个公共初值包括常数 .1/1/4、平滑随机场、相关随机场、区域块。每对象各方法完全相同。

固定单起点协议每起点 45 秒、最多 2000 原物理 LU / 60000 复 RHS；K=1/4/8 从同一个 K8 运行的前缀分析，不重新跑相同的 K1/K4。

固定总预算协议将 45 秒和上述物理上限分配给 K4/K8 个起点，每个短起点独立执行。不能把单起点协议累计的 8×45 秒成功当作 45 秒部署收益。闭账时可能为返回 incumbent 多付一次物理作用，实际超支保留。

原数据 top1 与离线 oracle best@K 分开。材料与数据成功分别报告；独立统计单位是形状对象，初值、噪声和照明相关。

## 5. 费用与局限

Stage 1 是 64 项、同 45 秒 wall 上限的开发筛选，未固定相同 LU/RHS。新版本再加入 LU/RHS 上限；各 source freeze 明确区分。固定标量反馈不需要完整方向切向，后续版本移除了 Stage 1 为机制记录而支付的冗余切向诊断，不能直接把跨版本 wall 变化归为新的优化机制。

本轮完整切向拟合使用原 LU，不是 A17 Schur 加速。Schur 有单独的误差/投影诊断；未解析尾项保留。联合电流 CIE 原型为另一个明确标记的标准方法探索，不等同当前 Gram 目标。

运行入口：`python -m src.adaptive_solver stage2_remote.json`。精确本机、远端和数据路径见 REPRODUCE.md。
