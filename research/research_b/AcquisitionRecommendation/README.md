# 可复用的真实材料交接

从[E3/INDEX.json](E3/INDEX.json)选择四个固定真实A parent与三种B输出。完整数组与raw日志在本阶段[证据Release](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-b-20261015)中，解压后目录与索引保持一致。

每份目录保留原A材料、原观测、gain-applied prediction、去gain的旧频F、单位/噪声/solver合同、新action与完整更新history。候选权重为none；outside-support质量未知。`HANDOFF_COMPLETION_RECEIPT.json`说明后续corrector已完成，并索引已安装的双频fullvalidated材料与全部失败端点。原history的PENDING字段是历史通用状态，不能替代外层完成回执。

不要按残差softmax生成posterior，不要再次计入旧观测，不要为新频重置同一世界的gain。Research A源目录没有被B修改，B输出也不代表A已接收新数据并自行重跑。

接入实现见`code/ab_adapter.py`；原B schema保留，A原生execution profile与B独立六频菜单用明确版本adapter转换。软件fixture测试不等于真实材料或新成像成功。
