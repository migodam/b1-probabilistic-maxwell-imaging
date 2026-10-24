# 前瞻 R4L1 四臂：只准备

入口 `tools/run_r4_interval_chain_locked.py`。未来独立 Book **IV**，实验 `interval_chain_locked_v1`，固定 `R4L1_00000`–`00007`，seed `20261025171`。每父对象每臂 1,024 步，full66 probit pCN beta=.002；目标、likelihood-only MH 与两端势区间、四次 refinement、固定精度 .05/.05 条件均沿用现协议。

四臂 Full / Richardson / global GMRES / fixed-accuracy global GMRES。独立进程、同初值/创新/U、共同 stored/canonical endpoint-once 证据合同。近似臂各自付 8×576 逆列 RHS，之后相同 R0/G/S BLAS 权利；每次真正的小型复 LS 预付 1 RHS，所有 Full 回退按原 8 源收费。组合上限 300,000 RHS，固定子预算 65,536 / 78,154 / 78,155 / 78,155，各臂 1,200 秒、3GiB。没有隐藏借用、重跑或增加 refinement。

准备文件 `FROZEN_PROTOCOL.json` 与 `configs/INTERVAL_CHAIN_LOCKED_RELEASE_TEMPLATE.json` 已存在，后者 **released:false**。没有创建真实释放，没有物理启动，也没有等待输入后自动启动的逻辑。

父线程未来需要先完成 `configs/FROZEN_LOCKED.json`、`results/contexts/locked/R4L1_*.npz`、对应公共 `results/contexts/cache/` 与原 `inputs/online/`。上下文取得必须有已关闭 COMPLETE 的 StageIV receipt。然后由父线程另建 released:true 文档，精确填写 prepared_protocol_identity 与 context_receipt。

先运行 `preflight --release <真实释放文件>` 检查，再显式 `pair --release <真实释放文件>` 启动一次。缺标记/数组、目标变化、probit/缓存/观测不一致、上下文作业未关闭、释放身份不匹配等均在 Book 前停止。协调器保存输入身份；各子臂也核对该初次身份，避免臂间上下文/释放变动。已有结果或首次尝试 marker 均禁止重放。

四臂都关闭后才调用共享比较，检查 Full 参照的 flags、x 路径、所有势 proof 和固定精度停止条件。float64 包络不是定向区间证书，数值相同不证明后验参考、混合或 MCSE，也不预设 Richardson 必须胜普通 Krylov。

fast 引擎只增加可覆盖的 ledger stage / context directory / cache directory，默认仍为原 II/DEV；旧 global 入口及封账 source_snapshot 不变。未来 Book 会保存当次 hooks、包装器、global 与 fast 源码，不以事后源码替代旧 as-run 证据。

零物理测试与当前缺失门禁记录在 `AUDIT/physics_mechanism_tables/`。这些准备检查不等于 locked 实验完成。
