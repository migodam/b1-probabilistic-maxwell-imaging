# 固定 observation 000 离线成像图库

四组图为实际保存参数/完整响应与冻结 decoder 的离线展示。完成门槛先于 GT 读取；只读 GT parameters 第0行，不读真值F或生成计划。

候选世界列质量属于两个公开墙分支与目标位置区间事件；不是单个材料点概率。材料均值/方差使用有限联合posterior，不能视为合法第三世界。

目标图的95%区间来自有限原子权重，不宣称连续覆盖率；CDF保留原子step，连续参考只显示独立积分检查点。预测均值使用相关性保留的E[gamma F]。

650材料节点仅流式累计目标覆盖的场矩，公开墙各解码一次，没有场历史堆叠。模型外质量为null，一般unknown-wall在线成功不在图的主张范围，原minimum_feature语义保持未解决。

具体来源与主要数字见SOURCES.json、MAJOR_NUMBERS.csv及其余表；GT访问范围见GT_ACCESS.json，CPU/墙钟/峰值RSS见RUNTIME.json。
