# 低频原因：EcoQoS、有效电源模式与OEM可见性

**没有查到可以立即修正的已确诊OS配置错误。机器实际已处于最高性能电源模式；8个科学worker没有显式开启EcoQoS。** 忙P核约1.5–1.7GHz的计数器证据仍然存在，热、功率、OEM策略或自动QoS的具体原因没有被当前接口区分。

此次于2026-11-14 04:42–04:44 UTC只读审计112例真实队列。没有修改进程、线程、方案、OEM档位、BIOS、驱动或科学预算。

| 项目 | 现场证据 | 解释与边界 |
|---|---|---|
| 实际8个worker的进程节能状态 | 04:44:02 UTC确认8个单核亲和性worker，覆盖0/2/4/6/8/10/12/14；全部GetProcessInformation成功，Version1、ControlMask0、StateMask0 | 系统管理，未显式选择EXECUTION_SPEED节流/EcoQoS；不是显式强制HighQoS，也不排除系统自动QoS或线程级策略 |
| worker优先级 | 全部Normal | 独立信息；不能由普通优先级推断EcoQoS关闭 |
| 用户配置AC电源模式 | `ded574b5-45a0-4f42-8737-46345c09c238` | 已为Best Performance |
| 有效电源模式通知 | 版本2注册成功，当前回调4；正常注销 | SDK中4=MaxPerformance；已核对别名取值，非将文档列表序号误作枚举数值 |
| 当前供电与省电 | AC=1；SystemStatusFlag=0；电量100% | 接交流电、省电状态关闭 |
| DC配置 | Best Efficiency GUID | 当前未用DC，不能解释这次AC工作状态 |
| 传统方案 | 平衡 | 与有效覆盖层不同；先前单凭方案名称推断未开最高性能不成立 |
| 最高处理器状态 | 类0/1/2 AC/DC全100 | 无99%限制 |
| 显式最大频率 | PROCFREQMAX/MAX1/MAX2 AC/DC全0 | 本方案无显式MHz上限；本机拓扑仅见效率类0/1 |
| boost | mode2 | 已允许boost |
| 传统方案EPP | AC33，DC50，类0/1/2相同 | 方案静态值，不能替代实际MaxPerformance覆盖层后的有效硬件请求 |
| 旧ProcessorPowerInformation | 所有逻辑核MaxMhz=2200、MhzLimit=2200；CurrentMhz约1466/1512/2200 | 旧接口名义参考，未报告低于该参考的限制；不是“硬件最高只能2.2GHz”的证据，更不排除别的功率/热限制 |

微软的进程节能接口中，EXECUTION_SPEED掩码为1。ControlMask=1、StateMask=1才是明确选择EcoQoS；ControlMask=1、StateMask=0是明确关闭该节流；本次0/0是系统管理。微软明确说明，未显式选择时系统仍可根据启发规则判定QoS。因此本次排除的是“worker明确请求EcoQoS”，不宣称排除了全部后台QoS行为。[进程信息类别](https://learn.microsoft.com/en-us/windows/win32/api/processthreadsapi/ne-processthreadsapi-process_information_class)、[节流结构](https://learn.microsoft.com/en-us/windows/win32/api/processthreadsapi/ns-processthreadsapi-process_power_throttling_state)、[QoS行为与HighQoS示例](https://learn.microsoft.com/en-us/windows/win32/api/processthreadsapi/nf-processthreadsapi-setprocessinformation)

电源模式取自有文档的 `PowerGetUserConfiguredACPowerMode` 与 `PowerRegisterForEffectivePowerModeNotifications` 当前值回调；不是读取传统方案标签。回调只订阅当前状态，随后注销，不改系统功耗策略。[AC模式说明](https://learn.microsoft.com/en-us/windows/win32/api/powrprof/nf-powrprof-powergetuserconfiguredacpowermode)、[有效模式通知](https://learn.microsoft.com/en-us/windows/win32/api/powersetting/nf-powersetting-powerregisterforeffectivepowermodenotifications)、[SDK枚举及相同数值别名](https://github.com/microsoft/win32metadata/blob/main/generation/WinSDK/RecompiledIdlHeaders/um/powersetting.h)、[省电标志](https://learn.microsoft.com/en-us/windows/win32/api/winbase/ns-winbase-system_power_status)

## OEM与XTU

HP筛选后的WMI设置仅返回 `Fan Always On=Disable`。这不是当前满载风扇停转的证据，也不能代表OMEN的Performance/Balanced或风扇Auto/Max状态。HP Omen/App Helper服务和XTU服务均运行。`hpqBIntM`及相关缓冲类存在，但没有直接返回当前OMEN档位的明确属性；没有调用其不透明BIOS命令或任何设置方法，也没有读取密码字段。

HP官方说明，OMEN Performance档位会调整CPU/GPU功率，手动风扇偏好也可改变可用性能；因此它是值得读取的缺失状态，但尚未观测，不能预先确诊为静音模式。[HP温度与性能控制说明](https://support.hp.com/lt-en/document/ish_4052969-4031879-16)

XTU已知公共日志目录存在，最近文件的更新时间为2026-10-28 22:03 UTC，早于此次约16天，不能充当当前温度/功率遥测。仅列最多12个文件的元数据，没有全盘或个人文件扫描，没有解析未文档化的PerformanceSettings.bin，也没有安装或启动驱动。没有找到可直接复用的现成实时传感器API；CPU温度、封装功率、PL1/PL2、热/电流限制原因、当前OEM档位与风扇RPM均保持null。[旧处理器功率结构说明](https://learn.microsoft.com/en-us/windows/win32/power/processor-power-information-str)

## 后续判断

没有证据支持现在“切最高性能”“把最大状态改100%”“打开boost”能解决问题，因为这些状态已经存在。也没有依据在当前公平比较中强设HighQoS或调EPP。若之后需要验证自动QoS，应在另一个明确记账的工程对比中，保持工作量、核心和算法相同，只改变一个可回退的QoS请求；这仍是待验证假设，不是本次查出的修复。

最有价值的下一项证据是**现有OEM/XTU界面或已支持只读导出中的实时封装温度、功率及限制原因，并同时读取实际OMEN性能/风扇档位**。没有这些，不能区分散热限制、低功率档位、电源适配器识别、XTU/OEM策略和后台QoS。当前冻结版本比较继续原P8及原设置，避免中途改变CPU工作预算；E核/HT的吞吐选择属于后续工程配置，不由本次核心数直接承诺加速。

## 成本、失败与证据

三次只读往返合计10.373秒，本机控制器CPU合计约0.006秒、传输子进程CPU约0.029秒；远端CPU没有完整单独计量，保留未知。第一轮尝试的两个旧overlay导出入口在ctypes绑定时返回访问异常，被作为无效观察保留，随后改为上述有文档API成功读取；没有用失败读数判断模式。一次进程扫描遇到任务交接只看到7个worker，第二次确认完整8个，均有原始回执。设置变更0，新增科学求值0。

结构化证据：`CAUSE_EVIDENCE.json`。私有原始采集：`CAUSE_RAW.json`、`CAUSE_FOLLOWUP.json`、`PROCESS_CONFIRM.json`及对应脚本，位于 `private_operational/cpu_frequency/`。本报告更新了轻载/忙载报告对有效电源模式未知的范围。
