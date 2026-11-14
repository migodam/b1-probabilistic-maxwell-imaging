# 限频事件日志：最后一次窄范围只读检查

**没有找到指定的固件限频/热关机事件；事件证据未能进一步区分约1.6GHz的原因。** 这不证明没有热、功率或固件限制。

查询时刻：**2026-11-14 04:56:22.701 UTC**（新加坡12:56:22.701）；起点为**2026-11-07 04:56:22.701 UTC**，向后7天。只查Windows **System**日志，限定provider、ID和起始时间；每个窄查询最多20条，实际总返回0条，没有读取整份日志或个人文件。

| Provider | ID | 结果 |
|---|---|---|
| Microsoft-Windows-Kernel-Processor-Power | 37 | 0条，NoMatchingEventsFound |
| Microsoft-Windows-Kernel-Processor-Power | 35/54/55（补充电源诊断ID） | 0条，NoMatchingEventsFound |
| Microsoft-Windows-Kernel-Power | 86（临界热事件关机） | 0条，NoMatchingEventsFound |

PowerShell正常结束；返回的是“没有符合条件的事件”，没有权限拒绝或查询语法错误。由于没有事件，没有可列出的事件发生时间，也没有可与当前科学任务对应的事件。查询窗口包括已确认的真实P8任务观测 **2026-11-14 04:38:22.707–04:39:27.413 UTC**以及之后仍运行的112例比较阶段。它不能排除事件未产生、未记录、被清理或不属于所选ID/通道；没有向其他日志通道扩展。

Intel说明，Event37涉及平台限制OS可使用的处理器频率范围；本次未出现这个可支持固件限制的记录。[Intel Event37解释](https://www.intel.com/content/www/us/en/support/articles/000090904/processors.html) 微软列出System/Kernel-Power ID86为临界热事件关机；其缺失也不能排除未达到关机程度的热节流。[微软热事件诊断表](https://learn.microsoft.com/en-us/windows-hardware/design/device-experiences/examples--requirements-and-diagnostics)

OEM元数据直接复用同轮已完成的窄范围只读结果，没有再次扫描：HP只暴露 `Fan Always On=Disable`，当前OMEN性能/风扇档位仍未知；该设置不能被解释成满载风扇关闭。有效AC模式已确认为MaxPerformance；其余已知配置见 `CAUSE_REPORT.md`。CPU温度、封装瓦数、硬件限制原因仍保留未知，不能依据日志缺失填0或宣布正常。

本次只读SSH往返约**1.16秒**，具体CPU与成本见 `EVENT_LIMIT_EVIDENCE.json`；设置变更0、压力测试0、新科学求值0。当前冻结P8比较继续原资源/设置。没有新增可据此立即实施的OS修复或硬件操作。
