# 交付快照与当前G接口

主证据zip保留实际执行时的代码/原始数组快照（源码5c0e8a9），不修改任何已付费材料、Gate或账本。公开发布后检查G接口时修正了一个纯类型分派问题：NumPy材料数组也有.data属性，不能据此当作已求解state。修复同时覆盖validate/correct/compensation_score；真实搜索轨迹一律传owned state，既有实验和48份材料不受影响。必要单项回归已通过，新增state/RHS为0。

G应使用仓库当前main的src/research_t，或在提取主证据后应用同Release的research-t-20261015-interface-overlay.zip。该overlay只包含当前代码、说明和纯接口检查，不包含或改写实验数组。

推荐从当前仓库checkout读取：解压证据时使用unzip -n（不覆盖已有源码/报告），只补齐未纳入Git的大数组。不得在T封口账本上重跑求解。
