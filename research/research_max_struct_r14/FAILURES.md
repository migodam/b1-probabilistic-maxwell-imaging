# 失败、偏差与修正记录

## 1. 远端首次验证漏载既有vendor依赖

`wiring_v1_CONTROLLER.log`在导入threadpoolctl时失败，尚未执行数值模型。新入口加入已有legacy/vendor，未安装依赖、改物理或优化算法。v2通过原A+-Arc接线及重入核对。失败wall14.2958125秒；该失败进程CPU未捕获，写null而不是0。

## 2. 8²局部机制首次配置被原接口拒绝

原arm_config只接受mf12 continuation，不能直接用mf4。保留首次失败及已付1次64阶LU/6RHS；改为先继承原默认配置，再显式设定声明过的单stage mf4机制实验。没有改变材料/方向/epsilon来寻找结果。

## 3. 主四组合corrector预算偏差（评分前发现并公开处理）

最初书面合同写proposal最多72 CPU秒、corrector最多48 CPU秒。新2D runner实际在proposal完成后设置`deadline=120`，会重分配提前结束proposal的空余时间。旧C2 runner则设置`min(120,start_corrector+48)`。这是实现与计划的真实差异，不能只用相同120标签掩盖。

处理：保留192条主四组合原始预测与全部费用，按实际“72 proposal / 总120剩余可校正”的统一四臂实验解释，不称与旧R13结果停止规则完全相同。另从已付费池报告旧48秒corrector口径的保守/区间前缀，标明非重新停止部署。不得用预算外/未来候选比较。

新独立24尚未运行预测、GT未读时，统一全部10方法使用同样总120策略。成熟C2数值代码只作一处预算延长：`deadline=min(cap,time.process_time()+cap*(1-proposal_fraction))`改为`deadline=cap`；原文件完整保留，逐字差异公开，其他参数、prior、选根、初始化、表示及停止容差不改。这给成熟基线同样的额外corrector机会，避免削弱它。独立比较同时报告完整CPU与原48口径前缀。

3D六臂另已统一冻结严格最多48 CPU秒corrector，不受此2D修正影响；跨物理不直接比较绝对费用。以上处理在主48新GT评分和独立24预测之前确定，不是根据恢复结果修改规则。

## 4. 扩展任务metadata的预测前审查

准备器的一行缩进把普通四pair的pair字段更新误置于stage候选if内；旧methods/task_count metadata也未包含父线程追加的global-Dual。父线程在任何扩展预测前阅读实际代码发现并修复，补齐所有任务显式pair/metric字段、实际10/3/2方法清单及计数；没有失败逆解或GT后调整。复现准备器也同步修正。

## 5. metric资格结果输出的重复关键字

解析whole-box界计算完成首项后，构造输出dict时重复scope关键字导致TypeError。改输出字段kernel_scope后完成2D/3D所有既定参数资格；没有改数学、阈值或参数。首次失败toolwall0.409276秒，进程CPU未捕获，0forward/LU/RHS。

## 6. shared-readout控制首次缺少psutil

本地既有numerical Python没有psutil，16个首次任务在import时失败，0LU/RHS。保留attempt01日志及OS child CPU合计0.514267秒；attempt02只复用已安装的ABI3 psutil包，不安装新依赖、不改变solver。16次成功结果、各自完整CPU和原失败分开保存，不能把同task_id的两次执行混成一次。

## 7. 仅评分和运输的问题

- 首次解包遇到本地/远端配置JSON的换行差别，拒绝覆盖；随后确认解析对象完全相同，保留本地原文本并在解包回执逐项记录。材料数值没有变更。
- `analysis/exposed48_v1`完成评分，但导出native行时重复`pair`关键字，产生192个导出issue。修正字典合并后新目录`exposed48_v2`重做离线评分，0issue；v1原结果及评分CPU保留，不重跑求解。原saved views、选根和预测seal未变。
- 审查extension scorer发现physical counters被规范化两次，会让RHS丢失；在任何extension评分前删掉第二次规范化。在线算法和原始费用回执未变。
- 原主gallery逐行colorbar与末列重叠。另生成`FIGURES/exposed48`专用colorbar列，保留原图；只是可读性修正，数组、指标和选根未变。

## 8. GitHub draft的tag查找

创建未发布draft后，按tag读取REST接口返回404；改从已认证release列表取得真实draft ID。原始证据上传仍面向同一draft，不将draft链接称作已经公开。没有科学重跑。

## 9. 非科学流程的失败和效率修正

原档规划首次对已知生成的多GB数字JSON做冗余文本扫描，尚无ZIP时由主线程停止；已测CPU下界473.875秒、最终尾部未知，原件未删、0新physics。后续跳过已知数值目录/JSON/CSV扫描，完整原始bytes照常归档。另一轮在公开psutil URL上的扫描误报已修正边界，失败规划保留、CPU未知。

独立评分第一次runpy启动漏了本轮src导入路径，在读取truth前import失败；0physics、CPU未捕获。补齐入口后完整264评分完成，0issue，原预测和选择不变。

一次只读Git尺寸盘点触发partial-clone历史blob下载，因与研究无关而仅停止本任务该读取进程及其子进程；没有重设分支、删除原件或改历史，CPU/网络byte未知。随后只枚举本轮文件名/实际大小，不再为尺寸盘点拉取历史。

事实复核纠正文献[4]题名，精确补充0019阈值案例，并把“必须早期11”限于已测schedule；这些不是新的科学运行。最终公开发布时态以真实Release/Pages匿名回读为准。
