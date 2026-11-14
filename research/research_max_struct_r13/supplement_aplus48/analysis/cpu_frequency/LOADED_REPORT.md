# 真实 P8 科学任务中的低频旁路观测

后续原因审计见 `CAUSE_REPORT.md`：已确认有效AC模式为MaxPerformance，8个worker进程为系统管理QoS（未显式开启EcoQoS）。本文件的候选原因与建议须结合该更新读取。

**这次忙P核的计数器也偏低，不能再仅用“空闲E核拉低全机平均”解释。原因仍未定位。** 2026-11-14 04:38:22.707–04:39:27.413 UTC，旁路观察已授权的112例冻结版本/接口比较队列；没有另起科学实验、压力测试或改设置。13个Windows性能快照，间隔约5.39秒、跨度64.71秒；同时记录研究根目录匹配的Python进程CPU累计时间和亲和性。

| 读数 | 结果 | 能说明什么 |
|---|---:|---|
| 8个选中P逻辑核停车 | 104/104读数为0 | 这8核均未停车 |
| P8 PercentProcessorTime | 中位100%，范围14–100% | 多数取样点忙，少数任务切换/等待 |
| 根目录匹配worker CPU时间增长 | 各稳定观测段约0.83–0.98核 | 科学worker在实际消耗CPU，非空队列 |
| P8 ProcessorFrequency | 中位1466MHz；均值1544；范围1466–2200 | 忙P核也报告低频，不能只怪全机平均 |
| P8 PercentProcessorPerformance | 中位69%；均值76.06%；范围64–152% | 大多数忙P核低于名义性能参考 |
| 2.2GHz × performance% 的估计 | 中位1.518GHz；均值1.673；范围1.408–3.344 | 明确以P基频2.2GHz为参考的**计数器推算**；未取得直接有效时钟遥测 |
| 全机 ProcessorFrequency | 中位1646MHz；范围1533–1734 | 与用户约1.6GHz现象相符；未直接截取任务管理器显示 |
| P8 PerformanceLimitFlags | 全部0 | Windows此计数器没给出限制原因，不能证明无PL/热限制 |
| P8 PercentPerformanceLimit | 全部100% | 不能替代硬件节流原因传感器 |
| AC供电 | 已连接 | 非当前DC状态 |

进程亲和性与冻结配置相符：worker各自绑定0、2、4、6、8、10、12、14之一。观察期内有任务退出、新worker接替，因此每进程CPU平均值按其自身可观测区间计算，不把不同PID连成同一任务。HT兄弟全部停车，E核停车比例约76%；E核仍有系统/其他工作，未枚举无关进程。少数Python管理进程允许全部32线程但CPU增长接近零。

本机之前已读出：最大处理器状态100%、PROCFREQMAX/MAX1为0、boost mode2、AC EPP33。**未发现传统方案的显式MHz封顶、99%上限或禁用boost**；因此这三种常见解释缺少支持。仍可能涉及OEM性能档位/Windows电源模式覆盖、温度/封装功率/电流限制、XTU已有策略或固件约束；本次数据不能区分它们，也不能宣称其中一种已确诊。现有XTU服务不等于已开放可读温度接口。温度、封装瓦数、PL1/PL2值与热/功率限制原因均保持null。

后续先用现有XTU/OEM界面的只读监测核对：忙核有效频率、封装功率、核心温度及当前限制原因；再按实证决定是否调整OEM性能/风扇档位或AC EPP。若确认无热/功率约束，低EPP可作为可回退的单项候选；若已触限，先处理供电/散热/并发。不要关闭保护、强制固定核心频率或无依据提高功率限额。

用户已授权后续合理利用E核和HT；**本轮仍保持原P8配置，是为版本/接口比较保持相同CPU资源**。换32线程会同时改变吞吐/功率/温度，妨碍当前公平比较。下一次工程吞吐选择可以使用已验证的E16–31及HT兄弟，仍需实际整批耗时与成本判断，不能由核心数直接推出加速。

采样成本：成功SSH往返69.877秒；远端Python采样CPU1.547秒、远端采样墙钟69.159秒，本机控制器CPU0.001341秒、SSH子进程CPU0.009702秒。另有一次54.448秒的失败采样及3.087秒的自有采样器回收；总往返时长相加127.412秒，其中回收与失败采样部分重叠，不能当作顺序新增等待。失败是采样输出管道处理问题，只停止了唯一命中采样专用命令的自有PowerShell进程，没有停止科学worker；失败回执与成本已保留。PowerShell自身CPU未单独计量，不将它记为0。新增物理求值0，设置变更0。

结构化证据：`LOADED_EVIDENCE.json`，原始逐逻辑核快照/进程与成本留在 `private_operational/cpu_frequency/LOADED_RAW.json`、`LOADED_SAMPLER_COST.json`、`FAILED_LOADED_*`、`RECOVERY_COST.json`。之前的轻载报告是历史快照，应结合本报告读取。

参考定义：[Intel本机型号的2.2GHz P基频、1.6GHz E基频及睿频/功率规格](https://www.intel.com/content/www/us/en/products/sku/235995/intel-core-i9-processor-14900hx-36m-cache-up-to-5-80-ghz/specifications.html)、[微软EPP的性能/节能方向](https://learn.microsoft.com/en-us/windows-hardware/customize/power-settings/options-for-perf-state-engine-perfenergypreference)。最高单核5.8GHz不是持续P8频率承诺。CPUperformance百分比推算与Windows ProcessorFrequency是不同观察量，均未被提升为硬件直接有效时钟证据。
