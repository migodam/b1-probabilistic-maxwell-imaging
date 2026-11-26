# 本轮数据登记

主数据族 bmtm26_reciprocal_blind 含 8 个 32×32 数值材料与另外 2 个 96×96 数值材料。每个目标有 3 或 4 个独立改变位置、角度、形变与空间介电常数的散射体，含凹口、空心结构和非均匀区域。随机种子在正演前确定，未进行 Born 误差筛选。几何生成及事后真值评分在 ports/private/；盲求解只能读取 ports/public/。本轮协调者设计比较，执行者分别负责数据生成、盲重建和矩阵界核查。外部原生 CIE、Fourier 及脉冲算子继承同仓 MAX-STRUCT R8 工作树，源码随包给出。

| 读数类别 | 识别符与单位 | 使用范围 |
|---|---|---|
| 初始拟合与留出 | bmtm26_n{32,96}_cNN_f400_a08_c48，完整 ID 逐条见 cohort_manifest.json；负载波幅单位 √W | 8×8 拟合，另 8 个接收位置留出。噪声范数为清洁读数的 2%，发布的硬球为 3%。盲阶段只读给定 Y 和 r |
| 布局、互耦、标定和双频配对 | 各材料的 _mimo.npz，键含 f400/f460、a06/a08/a12/a16、c016/c048/c144 及误差代码 | 同一材料改变观测条件。仪器单通道绝对噪声地板固定。与初始 8×8 重复测量的噪声方向不同 |
| 成像区域 | imaging/results/材料/stationary/baseline.npz 的 rois 及 background_tiles 或等效存储项 | 在收敛图上自动生成；全路由按最终掩码逐格核对。真值仅在评分时打开 |
| 标定机制材料 | ports/results/mechanism_witness.npz 与 mechanism_class_collision.npz | 专为 6×6 旧测量寻得的竞争物理材料及其两个类。与十个独立盲对象分开报告 |
| 改换收发方向 | results/aperture_control/ 的几何与复数矩阵 | 复用同一八端两弧阵列及匹配网络，A、B 为等功率等曝光等通道条件 |
| Manitoba 原始微波实测 | 2010 年 7 月 University of Manitoba EIGOR 数据集 1，4 个实物，3–6 GHz，频率间距 0.5 GHz | 两次空室和公开的目标扫描，仅用于分析误差尺度；数据来源与原文见复现说明 |

[逐个完整观测 ID](ports/results/cohort_manifest.json) · [阵列合同与生成脚本](ports/DESIGN.md) · [盲反演冻结设置](imaging/FREEZE.json)

交错双弧复测中，A、B 两份读数各自保留完整的通道矩阵和各自半径。联合结果按一次候选须同时解释两次测量统计。阵列数与标定回放的盲宽度比较统一引用最初稳定 CIE 图的第二区域，逐阵列读数和主 8×8 读数的噪声方向分别冻结。
