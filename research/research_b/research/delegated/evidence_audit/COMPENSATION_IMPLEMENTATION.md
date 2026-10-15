# E1 补偿机制模块实现合同

这是尚未执行的实现说明。模块由主线程在已授权数值作业内调用；本子任务不执行插值优化或付费物理。

- 入口：`compensation_experiment.run(response, paid_predict_callback, output_dir)`。
- Callback：`paid_predict_callback(latent, frequency, kind=...)`，返回单材料32实通道数组，或既有forward的tuple首项。Callback负责统一物理预约、历史cache、wall guard及全部失败计费；本模块不导入solver。
- 结构固定为 r=.08、r=.25中心空腔、h=.08环形空腔，reference为solid ε=2.5+.1i/gain=1。只优化ε与gain，几何不可改变。
- 按已冻结四个新频率1.2/2.0/2.8/3.2，与旧1.6联合。每个世界可以有不同合法gain；这里固定reference世界gain1仅为局部机制见证。
- 每个优化使用三个预定义starts、最多80次function evaluations；保存全部attempt。old-only/new-only/joint-shared/故意independent-frequency gain/calibration+old/calibration+old+new均明确输出。independent gain负控制仍共享ε，禁止每频独立材料。
- full验证只选全部12个joint-shared endpoints，固定顺序、不据收益筛选：每个2频=24callback；5个canonical reference频率共享，合计计划29，硬cap30包括失败。非joint优化只为proxy结果。
- raw nominal材料若没有full缓存，不冒称full raw验证。full callback所得joint endpoint可额外计算同一ε下old/new各自所需gain及校准约束后的gain，不新增forward；这类诊断不是重新优化ε的全局界。
- noise固定为每频canonical solid full响应real RMS的.12倍；不得按candidate缩放。Calibration使用冻结已知solid旧频field导出的gain statistic标准差σ/||Fref||。
- objective为无噪声reference数据的profile距离，不是evidence或posterior，且所有最小值仅为bounded best-found。calibration+frequency为2个新blocks，对比频率单block必须披露费用不同。
- 使用append-only callback attempt日志；成功响应恢复缓存，失败不自动重试，30cap跨同output_dir恢复保留。已有最终结果直接读取。主线程仍是整个campaign预算authority。
- 报告：`COMPENSATION_MECHANISMS.json/.csv`、`COMPENSATION_OPTIMIZER_TRIALS.jsonl`、`COMPENSATION_CALLBACK_ATTEMPTS.jsonl`、`COMPENSATION_RAW_PREDICTIONS.npz`。

执行前限制：所有四新频率必须已存在response bank；callback失败只产生null/失败记录，不允许fallback替代full标签。只读旧P2、A和source package。
