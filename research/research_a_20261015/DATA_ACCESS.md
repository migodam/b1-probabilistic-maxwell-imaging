# 公开代码、原始数组与复用资产

本阶段代码、中文报告、全32parent原生/common exchange、图库和主要原始结果保存在现有B1仓库本目录。完整证据资产包含全部新的Mie endpoint缓存、原始失败/旧metadata、所有完整材料/预测、复用的P2 bank与固定MDN、A30冻结后端和本轮8条件必需的旧背景/N3 LU/pools。

[Research A Evidence Release](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-a-20261015)

[下载完整证据包](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/research-a-20261015/Research_A_Evidence.zip)

ZIP根目录保留`src/`和`research/research_a_20261015/`；复用资产在`reuse_assets/`，每项标明历史来源，不重算为新增物理。C1的p2源码从原C1仓库8305688复用，或使用包内对应源码与同一.pt。A30的src/vendor和精选inputs按原相对路径保留；运行时显式设置`A30_ROOT`，不要重建一个相似Green冒充原target。P2 acceptedreference仅供离线参考，绝不作为online feature。

科学运行已封口。读取JSON/NPZ和图库不需要新物理；重跑full campaign会消费新授权，不能删除本次账本、改名run或利用旧P2剩余额度启动。解压到同一版本目录时保留已有原始数据，不覆盖历史Gate。没有摘要值核验或额外smoke流程。
