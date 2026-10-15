# 下一阶段：先共同物理接口与结构发现复现

1. 只读导入T最终Git提交的完整physical handoff；版本化核对forward、单位/packing/噪声、合法base、corrector/validator和proposal callback。不得拿正在写的工作树输出自动开新实验。
2. 冻结同一corrector后补齐global multistart、连续尺度static、T-guided、QD与轻bandit五臂；对目前p250073边界witness作独立起点/跨对象复查，不能只复跑同seed。
3. 区分raw未产生结构、保留但不相容、修正折返、有限chart不可达与archive错误剔除；一次只变proposal或archive，不能把改corrector带来的收益归给search policy。
4. 保留原Best-of-K和几何coverage目标，短K如实报告。p250077真洞尚未覆盖，应成为困难暴露回归案例，不按truth指导在线提案。
5. 当前不提交PPO/GRPO训练方案。先有简单策略在多个对象上持续发现有意义新结构，并显示prefix历史可用于选择动作，再讨论学习成本摊销。

本轮E2大factorial、E3 RL、E4一般3D均NOT_RUN。Pro建议512/1944/4096等更大预算没有授权，本轮余额也不自动变为下一阶段预算。这里是未来方案，未执行。
