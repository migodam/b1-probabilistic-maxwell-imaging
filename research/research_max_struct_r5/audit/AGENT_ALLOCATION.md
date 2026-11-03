# 本轮原生Agent分配

主线程负责推导、研究问题、科学转向、最终原创性和结论。工作子任务均显式选择gpt-6.1-sol，不继承默认高成本主线程配置，不启动DeepSeek或更换项目provider。

| 任务 | 模型/思考 | 职责 |
|---|---|---|
| state_inventory | Sol / medium | 真实状态库存、封存后独立GT评分 |
| backend_inventory | Sol / medium | 算子/九数字/历史证据提取、成本和复现 |
| literature | Sol / medium | ScholarQA和原作者文献核对、公开数据资格 |
| trajectory_replay | Sol / high | 精确历史接线、有限位置消融/交接、既有矢量后端小例 |
| diagnostics | Sol / high | 共享实材料导数、Gram/HVP、同状态反事实 |
| review | Sol / high | 独立公式/实现及最终主张审查 |
| figures | Sol / medium | 只读封存数组绘图与视觉检查 |

这些任务分批复用有限槽位，并非同时运行七个worker。理论判断、算法方向选择与最终审查仍在主线程。Agent费用/token没有可读计量，不估造金额。实测新求解全部本地CPU，原远端i9/4060队列没有占用。
