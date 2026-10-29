# E3 墙返回表示资格：实现交接

状态：**第一轮已在主线程释放后完整运行：708 RHS、19 FE，64 行全部保留。** 另做一次固定小矩阵代数核验（PDE=0，最大误差 3.52e-16），不作为科学证据。未写旧 GQ、vendor、合同、预算或发布端。

第一轮 raw 输出在 `results/e3_return_operator/`；wrapper 回执为 `results/jobs/return-operator-frozen/RECEIPT.json`。主线程随后明确授权仅修复复对称算子的右侧投影；独立脚本 `tools/return_transpose_repair.py` 已复用第一轮缓存完整运行，实际追加 256 RHS/3 FE，并保存于 `results/e3_return_transpose_repair/`，第一轮未覆盖。两次表示运行合计 964 RHS/22 FE；实测 review 见本目录 `RESULT_REVIEW_ZH.md`。

实现文件：`src/r2/return_model.py` 与 `tools/return_operator.py`。唯一物理入口为新 `r2.physics.Eval`；其 `factor` 先预约 LU，返回因子的每次 `solve` 自动先预约 RHS 并逐列核验 backward error。本实现不重复预约、不使用原生 SuperLU 绕过计费。

## 冻结范围与顺序

- `public_wall(0)` 与 `public_wall(2)`，相同宏观公共网格；`f=1`、S0 两源、16 接收器。
- 合同中的 2 train + 前 2 heldout 目标。每个目标保留原 cut-polygon 解析入射积分 `bX`，不替换为 `V * nodal_incident`。
- 全公共目标 ROI 的自由 FE DOF，维数由后端给出；不使用当前目标的 touched mask 作为端口。
- 每墙独立 `H=P(Gw-G0)P^T`，固定 seed `28102931`，32 个复标准 Gaussian 列，oversampling=0；同 seed 重置使两墙使用同一 Omega。
- `Q=QR(H Omega)`，规范复相位，保留随机列顺序；秩 4/8/16/32 为嵌套前缀。未做事后最优秩选择或重新排序。
- 另付 `Gw Q`、`G0 Q` 构建 `Hcore=Q* H Q`。保存的奇异值仅属于采样 action 矩阵与小 core，不是完整 H 的谱。

已经阅读原 Pro TeX §8.1–8.4；实现按第 8.1 节的算子顺序保留完整目标自作用：

```
q   = bX + V uw
j0  = q + V (L0-V)^-1 q
KQ  = VQ + V (L0-V)^-1 VQ
j_r = j0 + KQ_r core_r (I-Q_r* KQ_r core_r)^-1 Q_r* j0
F_r = dWall + Bw j_r
```

端口注入与限制在代码中显式完成。没有完整 C0 缓存，也没有把自由目标自散射投影到 Q。

## 对照、真实反馈与指标

`Bw=R Gw P^T` 由每墙 16 个 receiver-adjoint RHS 构造。每墙/目标另做 fresh direct full `(Lw-V)u=bw+bX` 参考。每个方法/秩再付 2 RHS 计算真实 `Gw j_r`；这些 action 只用于校验，没有回流修改基底或 rank。

报告完整 FE residual、`j-(q+V Gw j)` 电流闭合缺陷、参考电流误差、自由 FE 与公共 ROI 场误差、目标诱导场误差、白化响应 norm/RMS、目标增量响应相对误差，以及接收伴随映射和真实重构的一致性。小闭环矩阵的全部奇异值、2-condition、inverse norm 为数值诊断，不是物理 resolvent 上界。

整体 C 对照使用**相同 H 拟合的 Q** 与已付费 `Q* Cw Q`，求解 `j=q+V Q core (I-Q* V Q core)^-1 Q* q`。它不保留精确目标自散射；不是最优 whole-C ROM、独立 C rangefinder 或一般普通算法下界。

QR 与 action/core 奇异值始终采用 **Euclidean FE coefficient coordinates**。另保存固定公共目标 ROI 的稀疏 P3 L2 Gram 矩阵，仅用于场误差，没有质量白化、L2 算子谱或周期角谱结论。

## 资源

|部分|复杂 RHS|FE/LU|
|---|---:|---:|
|共同 free factor、两墙 factor/源/伴随/范围与 core|292|3|
|每墙每目标精确 free-target 自作用|272|8|
|8 个直接 full 参考|16|8|
|64 方法/秩行的真实 feedback 校验|128|0|
|合计|**708**|**19**|

每墙有自己的 Q，因此未采用最初低估的 4×34 模板 RHS；计为 2×4×34。两个实验阶段独立，本子任务不启动 cold Bayes。

最多同时持有 free、wall、当前 template 或 reference 三个大因子；每因子实际存储限制 1e9 bytes，活动 campaign RSS 限制仍由 4 GiB wrapper 控制。没有 m×m Green cache。执行前检查 E3 与总账本余额可覆盖 708 RHS/19 FE；失败预约与已完成行保留。

## 输出与执行

默认调用只输出 PLAN_ONLY，无 PDE。主线程释放后应通过新的 R2 `tools/run_job.py` 包裹执行：

```
<python> tools/run_job.py E3 e3-return-operator <python> tools/return_operator.py --execute-frozen
```

默认 `results/e3_return_operator/` 必须尚不存在。输出包含 sparse L0/Lw/V、公共 mesh/port/ROI-L2-mass、Omega/Q/原始 paid actions/cores/singular values、uw/dWall/Bw、精确模板、直接参考、每方法每秩的 current/field/response/closure，以及 PER_RANK_ERRORS JSON/CSV、SUMMARY 和失败记录。

进程实际 CPU 在 child 摘要中保留本地计时；唯一收费 CPU 为 wrapper 在 child 结束后写入的 RECEIPT。不会把运行中账本 CPU 增量 0 当成实际成本 0。最终科学解释、是否支持 §8 与任何后续概率阶段资格由主 Codex 复核。
