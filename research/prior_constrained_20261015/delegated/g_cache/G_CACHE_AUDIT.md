# G 既有缓存读取审计（E0，零新增物理）

读取来源：同级 `b1-research-g`，实际 HEAD 为 `54d7f463f9d1bad1cc5625f7ddfd8f449e0b7978`。本审计没有写入该 worktree，没有实例化 solver，没有读取 locked truth/geometry/clean-response，没有训练或重跑 G。

## 可复用精确数组

`inventory.json` 列出全部 291 个已付费 full-state 缓存：3 个已验证基础解，以及 72 条轨迹各自的 raw、correction1、correction2、fresh validation。每条含绝对数组路径、chi/prediction/current/field 的精确 NPZ member、observation、background、原始 sigma、case、method、round、parent 与 raw/endpoint 关联。可直接供本轮 E0 读取，不需要从报告反推数值。

仅访问三个历史 exposed development 条件：

| case | sigma | G base 白化残差 | base ell |
|---|---:|---:|---:|
| p250069_weak_a1_snr20 | 0.001270913712869431 | 20.929164702218824 | 219.01496756630118 |
| p250073_weak_a1_snr20 | 0.00119211798514519 | 20.891919121358494 | 218.23614228669234 |
| p250077_weak_a1_snr20 | 0.0011846108030513216 | 22.994967915752223 | 264.38427472323707 |

从精确 prediction/observation 及固定 sqrt(2)/sigma 重算，291 条 residual 与历史 metadata 的差均未超过 1e-8。没有按候选响应重估 sigma。历史 291 states/6936 RHS 只列 provenance，**本次新增 states=0、Maxwell RHS=0**。

`results/CONTINUOUS_REPRESENTATION.npz` 的 19 张 toy phase fields 是软件表示示意，没有 Maxwell prediction；明确排除出真实候选/竞争见证。

## 原 G 的先验与 corrector 合同

原 G 使用历史材料 BOX：Re chi∈[-.5,4]，Im chi∈[0,.5]，gain=1。全局重启为 base-centered、固定 Haar chart 的相关随机扰动，不是原 P0 的采样。局部加/减以已知 background 为参考，gaussian mask 连续位置/尺度，新增 target 为 background+(1.2+.12i)。没有原固定增量复相位约束，也没有原生成器 retraction。

corrector 是 Research A 冻结实现的 **32 个实参数 Haar GN/SLSQP**，最多两步，trust=.35，单步拒绝后 trust×.3；约束是逐像素 BOX 与 trust。它不是任意 1152 参数 full-space GN，但它覆盖全域并允许 DC/低频补偿；不保护 proposal 的结构，不检查原 P0 membership，也不约束增量相位。最后另付一次 full-state validation。普通 exact material Jacobian 已包含 resolvent feedback，不能将此实现称为人为删除 A17 feedback 的 baseline。

因此本审计的原先验标注为 `PRIOR_UNKNOWN` / `BOX_ONLY`；“原算法没有强制 P0”不等于逐个端点已经证明不属于 P0。具体 PHASE 与 P0 判定应由本轮统一先验审计进一步计算。不能将 G 的 box-valid/noise-ball-valid archive 升级为 P0_EXACT 或 posterior 高质量区域。

## Noise ball 与点似然的区别

- 72 个 raw 都未进入原 95% noise ball。
- correction1 有 2/72 进入；correction2 及 fresh validation 为 18/72。
- G 候选索引有 21 条 success entries=3 base+18 validation，包含 first/quality 重复。
- `EXACT_DUPLICATE_GROUPS.json` 以同一 observation 内完整 chi 的 `np.array_equal` 比较：291 缓存对应 111 个精确材料；72 validation 对应 36 个精确材料，其中 9 个 noise-ball-valid。不能把两轮重复验证算作新世界。
- 18 个通过 noise ball 的 validation 相对**同观测 G base**，Delta ell 最小 13.37066803742303、中位数 38.62147606053949、最大 64.45368925945664。全部超过 log(10)，也超过 log(100)。

最后一项是原 G 内部的历史诊断。G base 本身不具备 P0_EXACT 证据，不能直接冒充新轮正式 prior-compatible reference；但它明确说明“在同一 noise ball 中”并不表示两者点似然接近。上述没有计算 posterior odds。

## 补偿与结构保留

以下是 validation chi 减 raw chi 的物理数组量，不是先验验收标准：

| 指标 | 全72轨迹中位数 | 18个有效轨迹中位数 |
|---|---:|---:|
| 复材料改变量 RMS | 0.3000347689083427 | 0.279589907899162 |
| DC energy fraction | 0.08183326611746147 | 0.07744220065079954 |
| DCT 总 degree≤2 energy fraction | 0.3276505322043164 | 0.3291116735078606 |
| 改变像素数（abs delta>1e-10） | 576 | 576 |
| 沿原 raw-parent edit 的原 G 保留投影 | 0.015716130859402547 | 0.006687061556447644 |

G 原保留投影是完整复材料方向投影，不是本轮的多阈值/持久结构保留测试，不能直接代替本轮 eta=.5；但原正/负极小投影与全域补偿提示，raw novelty 与 corrected novelty 不能混算。对于局部动作，另在 inventory 轨迹记录 raw 改动超过峰值1%的描述性支撑与其外部 correction RMS；它不是 core hard-mask，也不是新增资格门槛。

## 现有结论边界

原 `GATE_DECISION.md` 已给出 `NO_GO_SEARCH`：只有一个近预定材料门槛、集中于一个 parent 的新增 world，未建立 QD/bandit 相对普通方法的持续优势；不是新孔洞拓扑或已证 disconnected branch。Best-of-K/coverage 没改善，T-guided、A17 新消融和一般非同心 3D 未运行。本轮保留原结论，不改写其 ledger 或 gate。

需要父 E0 继续完成：统一固定相位/幅度约束审计、原生成器支持证据、与 T/N3/pools 的共同 likelihood 对照，以及 P0_EXACT cached pair screen。G 的数组本身不足以建立本轮竞争解释的存在性。
