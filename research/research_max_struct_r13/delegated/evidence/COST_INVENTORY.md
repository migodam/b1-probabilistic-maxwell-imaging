# R13 已知费用清单：来源、重计规则与未知项

入口：`src/build_cost_ledger.py`。最新机械快照为 `delegated/evidence/cost_inventory_snapshot_v3/`，包含 `COST_INVENTORY.json`、`COST_INVENTORY.csv`、`COST_SUMMARY.json`、`COST_SUMMARY.md`。当前669条来源/区间/计数记录，28条应纳入消耗范围但CPU计时未知；这是运行中的元数据快照，不是最终完整费用。

脚本只读取实际队列、运行receipt中的成本字段、analysis/COST_TOTALS以及设置/验证/压缩等回执。不读取评分数值、材料/预测/真值NPZ，不调用正演、伴随或优化器，也不修改活动runner/config/queue。支持在最终归档后向**新目录**重新生成：

```text
python src/build_cost_ledger.py --root R13_ROOT --out R13_ROOT/delegated/evidence/cost_inventory_final
```

## 计时与去重规则

- 完成任务的外部CPU以终态queue jobs为权威；内部 `whole_process_CPU_seconds`只是覆盖到导出附近的引用，不再加入外部CPU。
- C1 `development_v1`已合并36条最终任务。`development_v1_recovery`的18条任务已包含在该合并记录内，全部作为引用保留，不再计一次CPU、物理调用或controller时间。
- 有终态回执就使用终态；否则只使用最新不可变STATUS快照，标记pending/running/无本地回执。未退出子进程的CPU不是零，也不以预算上限替代。
- `analysis/*/COST_TOTALS.json`只交叉核对相同online任务，不能与queue再相加。队列wall是并行历时，不能把child wall求和冒充它；C1干预间隔与active-controller窗口的不同口径保留原说明。
- 设置、online、controller、物理/机制验证、scoring/plotting和release分列；Apple与i9不汇合为一个速度数。
- 物理LU/solve RHS按矩阵阶数分列。`full_factorization_dimension`常是累计dimension计数，不能当成单矩阵阶数。576/1296生成调用已按实际计数拆分为48/12个LU；矢量sphere使用3N阶。
- projector H/port、profile QR/SVD、FFT、Green/receiver动作和current块RHS保留为单独的非加总引用，不冒充额外物理LU RHS。已有aggregate RHS优先，source/adjoint角色计数不再次相加。
- 生成/NN/诊断的outer区间已经包含case、setup、preparation或route子计时，子计时不重复加总。没有outer计时但有明确inner区间时，只纳入这个已知区间作为下界。

## 明确保留的缺失与失败

1. **C1四条中断尝试**：实际partial目录与任务ID保留，原process handle丢失，CPU/wall和完整物理总计未知。不能从重试或后来成功结果反推。
2. **dual-pressure序列化失败**：回执明确已完成3 LU、84 complex RHS（36 source+48 adjoint），CPU未捕获。这些调用是独立失败尝试，保留并计入已知物理消耗。
3. **最初本地psutil导入失败**：父线程披露发生在模型建立前；0物理调用，CPU/wall未知。原始计时未持久化，来源明确写作父线程披露，不伪造回执路径。
4. **生成脚本早期两次失败启动**：语法/缺NumPy发生在生成前，0物理调用，原CPU/wall未知；日志和LAUNCH_ATTEMPTS仍在。
5. **既有NN训练**：新的推理无再训练，但原监督训练CPU、GPU小时、训练样本正演费用没有在本轮重新测量，保持未知，不算零。
6. **scorer、旧图与表格聚合**：只有真实存在的计时字段才计入，例如mechanism_public/PLOT_RECEIPT；其余scorer/plot/bootstrap计时为空，不以“离线”视为免费。
7. **下载与部分runtime包装**：没有匹配计时回执时，仅承认文件存在，时间未知；字节数不推出下载时间。已知runtime打包与archive压缩分别列出，不能与online解算混合。
8. **Mie refinement首摘要与最终摘要重合**：first volume-only文件声明相同row receipts后来保留在带staircase控制的摘要中。清单暂按最终组合摘要一次计入，旧首摘要作为重叠上下文保留；没有擅自宣称两次独立的volume refinement。如父线程有独立重跑证据，可据确切attempt来源补账。

## 已核对的不变量

当前快照中，C1合并jobs的已知外部CPU与其权威终态记录一致（2800.71875秒），recovery引用没有再计；四个中断尝试仍为null；失败dual-pressure仍为3 LU/84 RHS/null CPU；生成维数拆分为576阶48 LU和1296阶12 LU。这些是成本记账检查，不是新数值实验或成像质量验收。

快照自己的metadata扫描CPU/wall也在JSON顶层保留。所有总计应解释为**按host与阶段分列的已知下界，外加明确未知项**。最终科学/费用与发布验收由父线程完成。
