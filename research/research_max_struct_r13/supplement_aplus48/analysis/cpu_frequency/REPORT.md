# i9 低频与核心利用：只读审计

后续真实P8忙载观测已完成，见同目录 `LOADED_REPORT.md`。忙P核也呈低性能计数器读数，更新了下文轻载阶段的证据范围；本文件保留初次轻载审计历史。

2026-11-14 04:28:58 UTC（新加坡 12:28:58）开始取样。没有增加科学实验、压力测试或新负载；没有改电源、BIOS、驱动、亲和性或其他会话。

## 结论

这台机器实际是 **i9-14900HX：8P + 16E、24 物理核、32 逻辑线程**。Windows 的 CPU-set 拓扑直接确认：P 为逻辑 0–15，HT 兄弟为 (0,1)、(2,3)…(14,15)；E 为 16–31。只选偶数 0–14 是 8 个 P 物理核，只选 0–15 才包含 P 的超线程。此映射来自本机当前 Windows 拓扑，不能套用于其他机器。[微软 CPU-set 说明](https://learn.microsoft.com/en-us/windows/win32/api/winnt/ns-winnt-system_cpu_set_information)

**本次不能确认“P 核工作时被压在 1.6GHz”**。它捕获到的是很轻的自然负载：总 ProcessorUtility 3%，CPU-set 中 27/32 逻辑处理器 parked；两秒 psutil 利用率合计约 1.37 个逻辑核。全机 ProcessorFrequency=1669MHz，与 P/E、空闲及停车状态混合有关；单个忙 P 逻辑核 6 同时报告 PercentProcessorPerformance=167%，而其 ProcessorFrequency=2200MHz。WMI 的 CurrentClockSpeed/MaxClockSpeed 都是2200MHz，也不代表硬件睿频上限。不同窗口、不同计数器不能互换为实际忙核有效频率，更不能据此直接宣布热降频。

## 当前电源设置

已确认接交流电，电量100%；HP OMEN 16-wf1xxx；只有一个传统电源方案，当前为“平衡”。传统方案读数不能替代 Windows 电源模式覆盖层或 OEM 控制中心状态，本次尚未读出这两者。

| 设置 | AC | DC | 判断 |
|---|---:|---:|---|
| 最大处理器状态，类0/类1 | 100% | 100% | 没有99%限制 |
| 最小处理器状态，类0/类1 | 5% | 5% | 允许空闲降频；不应直接改100% |
| PROCFREQMAX / PROCFREQMAX1 | 0 / 0 | 0 / 0 | 此方案没有设置显式MHz上限 |
| PERFBOOSTMODE | 2 | 2 | 已允许 boost，非禁用 |
| PERFEPP / PERFEPP1 | 33 / 33 | 50 / 50 | 有节能偏好，可能影响响应；不是1.6GHz硬限制 |
| CPMINCORES / CPMINCORES1 | 4 / 0 | 4 / 0 | 允许停车；轻载停车正常 |
| CPMAXCORES / CPMAXCORES1 | 100 / 100 | 100 / 100 | 未设置最大可用核比例限制 |

微软说明：[最大频率按MHz设置](https://learn.microsoft.com/en-us/windows-hardware/customize/power-settings/options-for-perf-state-engine-maxfrequency)、[boost模式](https://learn.microsoft.com/en-us/windows-hardware/customize/power-settings/options-for-perf-state-engine-perfboostmode)、[EPP越低越偏性能](https://learn.microsoft.com/en-us/windows-hardware/customize/power-settings/options-for-perf-state-engine-perfenergypreference)。本机枚举提供0值，因此不能把文档中的非零范围误读成0MHz。

本次性能计数器的 PerformanceLimitFlags=0、PercentPerformanceLimit=100；这是一次轻载读数，不能排除此前或重载时的PL1/PL2、温度、电流或OEM限制。Intel 的 5.8GHz 是最高单核条件值，不是32线程持续运行承诺；官方给出 P 基频2.2GHz、E基频1.6GHz、P最高5.8GHz、E最高4.1GHz、基础功率55W、最高睿频功率157W。[Intel 规格](https://www.intel.com/content/www/us/en/products/sku/235995/intel-core-i9-processor-14900hx-36m-cache-up-to-5-80-ghz/specifications.html)

## 传感器与证据缺口

- OpenHardwareMonitor、LibreHardwareMonitor WMI 命名空间均不存在；这表示未提供数据，不是温度/功率0。
- ACPI thermal zone 存在，原始3010/2832十分之一Kelvin，约27.85/10.05°C；没有证据它们是CPU核心或封装温度，不能拿来证明散热正常。第二个 CriticalTripPoint=0按缺失/未知处理。
- 已有 XtuService / XTU3SERVICE 正在运行；未调用调参接口，也未安装任何工具。存在服务不代表已读取CPU温度、封装功率或节流原因。
- CPU封装功率、P核有效时钟、核心温度、thermal/PL1/PL2限制、风扇/OEM档位目前均未知。
- 当时仅有两个低CPU累计时间的Python进程，普通优先级，亲和性都允许全部0–31；没有观测到正在持续运算的科学worker。因此这不是之前P-only满载时的重现。

## 可执行顺序

1. 在下一批已授权的真实计算旁路记录忙核计数器、CPU利用率、worker完成数/秒及线程亲和性，区分等待、单线程瓶颈、BLAS嵌套并行和持续低有效频率；不另开压力测试。若已有XTU/OEM工具能只读导出封装温度/功率/限频原因，优先复用；缺失则保留未知。
2. 独立单线程任务可将8个P物理核与16个E核纳入有限就绪队列；HT兄弟可作为额外候选。可使用0–31，但32线程不等于32物理核，不保证32worker最快。避免每worker再开启32线程BLAS；以整批吞吐、内存与尾延迟判断是否加HT/E，而非以总利用率或GHz判定科学/工程成功。
3. 若忙P核确实持续低频，且无热/功率限制证据，再考虑低风险、可回退的AC EPP性能偏好或OEM性能档位；一次一项并保留原值。现有最大状态100%、MHz无上限、boost已开，盲目“解锁”这些项不会解决已知限制。把最小状态/不停车强设100%也不是本证据支持的默认措施。
4. 若实际功率/温度已触限，应先检查供电、OEM散热/风扇与整批并发；不关闭保护、不提高PL1/PL2、不保证能维持最高单核频率。CPU+GPU共享笔记本功率/散热，加入CUDA后CPU频率可能下降，仍以整任务耗时与物理残差判断收益。

## 审计成本与文件

三次只读SSH往返约18.51 + 11.12 + 12.28 = **41.91秒**。第二次服务查询过滤表达式错误，AC/机型/方案仍成功；第三次只读修正成功。具体本机控制器CPU和传输CPU由已有 `audit/REMOTE_METADATA_COSTS.jsonl` 记账，未新增物理求值。网页资料查询及本地整理成本不包含在41.91秒中。

公开审计证据见同目录 `EVIDENCE.json`；私有原始回执与采集脚本留在 `private_operational/cpu_frequency/`，不含密钥内容。未承诺已解决历史降频；下一步需要真实忙状态的并行观测。
