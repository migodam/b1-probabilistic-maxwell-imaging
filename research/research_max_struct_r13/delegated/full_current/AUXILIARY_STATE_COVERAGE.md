# 辅助状态留存覆盖与已知缺口

结论：当前主 runner 留存了规定物理快照、全部已审计材料候选及其 metadata/cost，但**没有持久化完整 CIE/G19 辅助轨迹、每个 current/material 更新状态或 PR gradient/direction**。不能声称所有辅助候选都已保存，也不能承诺从这些产物精确续接原优化器。

本说明只核对源码与此前已暴露的产物约定；没有读取 fresh 材料质量/GT，没有修改活动源码，没有 replay 或重新运行在线批次。该缺口应作为交付范围保留，而不是通过事后重跑隐藏。

## 实际持久化覆盖

| 状态/证据 | 当前主批次留存情况 | 不能由此推出的内容 |
|---|---|---|
| initial/raw/corrected/primary/secondary 物理快照 | 有相应候选时保存 m、原 Maxwell 物理 j/Q、Y_prediction 和选择 metadata | 这些 j/Q 是对该材料重新求解的物理状态，不是 CIE/G19 当时的辅助 J/Q 或 DR iterator |
| 全部已付费、过预算物理候选 | `CANDIDATES.npz` 保存材料 m 池；`HISTORY.json` 保存各次物理审计 metadata、Phi/D/R、时间、合法/eligible 标记和累计成本 | 未选中候选的完整物理 j/Q/Y 数组并未逐一保存；也没有其完整辅助状态 |
| CIE/Dual/global-PR Dual/A+/G19/full APCS/full current | source/config、DPCS、denominator setup/basis、warm projection信息、current_spaces、累计 operation counts、termination及物理审计记录 | 主 runner 没有保存返回的 `stage_states`、native auxiliary `history`、每步系数 A/J/Q/m、old_grad/old_direction |
| NFFT-equations | 正常返回时，最终 joint J/m 作为 `FINAL_AUXILIARY_STATE.npz` 保存 | 并非所有 joint trials/accepted iterates；没有完整 L-BFGS curvature/line-search/gradient restart state；预算异常跳出时可能没有最终 joint state 文件 |
| 后续 physical Full-J TRF | residual trial的已付费材料池、物理快照、时间/成本；正常返回时message/nfev/njev/optimality摘要 | 没有逐步 residual/Jacobian矩阵、trust-region半径/active-set/内部线性代数restart状态；预算中断可能没有正常返回摘要，不能据此宣称精确TRF续接 |
| DR/AP/事件辅助 iterator | 正常完成固定循环后，最终 X/Q 保存至 `FINAL_AUXILIARY_STATE.npz`；另有物理审计、gap/cost/event/forecast metadata | 中断时可能尚未执行 final_states 赋值，因此没有终点 X/Q；也没有逐步完整 iterator/shadow/normal 历史 |
| 独立 exposed diagnostics | 其各自 NPZ 明确保存所声明的 c/d/e/p/correction 或短轨迹字段 | 它们是独立暴露例子的指定诊断，不能当作在线主批次缺失 PR/DR 状态的补档 |

“raw”在主 runner 中指 proposal 阶段按原物理 Phi 选出的**物理材料 incumbent**，不是原始 auxiliary J/Q。`HISTORY.json` 名称也不意味着它等于 CIE 返回的逐 outer-update auxiliary history。

## 源码证据

1. `src/run_experiment.py:196-198` 调用 scalar `run_full_current` 后只取返回 status；没有将 `result['stage_states']` 或 `result['history']` 写盘。periodic physical audit 经 `g19_audit` 重新记录物理状态，而非保存 auxiliary J/coefficients。

2. `research/research_max_struct_r8/src/cie_fft.py:306-319` 在内存中构造逐更新 scalar history；`:320-322` 的 iteration callback 可提供 m、auxiliary J、coefficients；`:332-350` 在正常结束路径构造并返回 stage-end m/J/Q/coefficients。`src/full_current.py:183-191` 主要用该 callback 调度物理审计，未持久化传入的 auxiliary数组；`:297-306` 返回结果，也不等于主 runner 已存储这些字段。

3. native PR 的 old_grad/old_direction 是内层局部变量，标准返回值中也没有完整 PR restart state。即使补存 stage_states，单凭 stage-end m/J/Q/A 仍不自动构成每个更新时刻的精确 PR continuation checkpoint。

4. `src/run_experiment.py:295-303` 只把规定物理 snapshots 的四个数组写入独立 NPZ，把材料池与物理/projection/event/forecast metadata写盘。这里没有遍历 scalar CIE 的 returned stage_states/history。

5. NFFT 正常返回后的 final joint 保存赋值在 `src/run_experiment.py:160-161`；DR/AP 的 final iterator 赋值位于固定循环之后 `:280`。预算异常通过 `:282` 跳出，不能把未执行的赋值当作已保存状态。已有文件是否存在须逐 run 核对，不能从支持该代码路径泛化为每个 run 均完整。

上述行号对应本次查阅版本，后续文件编辑可能移位；语义界面与字段名是本说明的主要定位依据。

## 当前仍能支持的审计范围

- 根据冻结 source/config、DPCS/Fourier/current-space设置、warm/cold政策、付费成本、termination、材料池与物理快照，审查 full-current 对照的实现、预算与材料选择是可行的。缺失 PR 状态并不将已保存的材料/物理结果自动作废。
- `audit_predictions.py` 可验证已有物理 snapshots 与 public G/E/S 的关系、初值一致性、材料索引、主/次/prefix选择和 metadata覆盖；它不能审查没有保存的 auxiliary J/A/PR gradient，也没有重算全部历史物理目标。
- 固定单 β quadratic profile 诊断是一个独立的指定 current-block 最优化诊断，不是对真实 legacy PR 轨迹、双 β G19 共同目标或缺失主轨迹的事后证明。
- 不能据此声称每个 PR current step 都有完整存档、已逐步验算全部辅助损失下降、可重建每次 readout 输入、或能从任意 CPU 截断点无损续算。

## 未来交接要求

若下一批研究需要逐步 current/readout 归因或 continuation，必须**在启动前**冻结所需检查点方案及存储/CPU成本，明确区分物理审计候选和 auxiliary iterator。可声明保存的 clean-boundary字段应包括 stage/update、native frequencies/coefficients、JD/DPCS与固定 source metrics、m/J/Q、beta、PR gradient/direction、相关计数和事件 RNG state；还需说明未完成 inner solve/line search 的处理。

仅保存 m/J/Q 并重启 L-BFGS/PR 应称为 warm restart，而非精确优化器 continuation。SciPy L-BFGS 的内部 curvature/work buffers等未保存时，不应承诺精确续接。完整逐步数组可能有显著I/O与存储成本，必须在相应预声明协议里计费，不能默认为免费观测。

本批次保持原始产物与已知缺口。没有为了使覆盖声明更漂亮而重跑或覆盖任何主批次；未来若确有科学理由重做，应登记为新的、事前instrumented批次，并保留当前批次与费用历史。

## 后续 TRF source-only补充

Parent后续 source `2cd5516cb`（ready `32bfdfa4c`）加入原 Full-J TRF分支。`run_experiment.py:127-143` 只把物理 residual trial送入普通record，并在正常返回时保存停止摘要；通用NPZ导出仍只保存 `m/j/Q/Y_prediction` 与材料池。`model.jacobian` 在当前物理 state内缓存的完整解析 Jacobian并未被导出，不能因为采用“Full-J”方法就写成“所有Full-J矩阵已存档”。这里的内部 `state['J']` 是观测对材料的解析Jacobian；物理NPZ的 `j` 是current，NFFT终点NPZ的 `J` 又是joint current，必须按文件schema区分。

该审查没有读取新TRF活动队列或任何fresh材料质量。TRF exact residual objective/prior的源码核对见 `BASELINE_SOURCE_AUDIT.md` §7；它不消除本文件的状态覆盖缺口。只保存最终材料后重新启动TRF应称为新的warm restart，而非恢复原trust-region内部状态。未修改活动代码，也未为补档而重跑任何批次。
