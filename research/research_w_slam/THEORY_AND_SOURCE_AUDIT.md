# 理论与来源核查：事实草稿

本文件只整理来源、代码接口和保存证据，不作最终科学门判断或新颖性判断。主线程负责核查科学解释。研究包的数学陈述不能因代码可运行而自动获得证明资格；离散恒等式、连续形状公式、局部信息和有限非线性结果须分别阅读。

## 来源边界

| 来源 | 精确定位 | 本文件使用范围 |
|---|---|---|
| 用户提供的 SLAM 理论包 | `source_package/Research_W_SLAM_Report_ZH.md`、`source_package/References_Main.md`、`source_package/Literature_Audit.md` | 理论主张与文献线索；本轮文献访问返回429，不能标成当前在线核验通过 |
| 历史条件见证 | 提交 `958a478739d9553269316ded8125fda9df225c0c`，`research/research_w_p0r_to_p1/` | 固定支撑、联合补偿机制；不自动属于新网格先验 |
| 历史全波后端 | 提交 `483d260c6b5b111db9fb814b096685fccb362448`，`research/research_w_continuation/src/wq/fem.py` | 新后端的明确复制来源；数值资格不随复制自动迁移 |
| 新算子适配 | `src/wslam/fem_backend.py`；适配说明 `research/delegated/assets/ADAPTER.md` | 新积分域、PML物理域、耗散上界；保留原 PDE/LU/残差与预约逻辑 |
| 小矩阵信息与高斯运算 | `src/wslam/information.py` | 标准线性代数的实现；不能据此宣称原创数学或真实后验认证 |
| 当前数值来源 | `results/e0/`、`results/e1/dynamic/` | 实际已保存的计算，不以缺失数据补足主张 |

上述相对路径均以 `/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/Research_W/research-w-slam` 为根。历史 Git 对象可由本地 `../publish-w-p0/.git` 读取。物理历史缓存实际位于 `../../B1/b1-w-p0r-to-p1/research/research_w_p0r_to_p1`；具体绝对输入清单见 `research/delegated/assets/inputpaths.json`。

历史见证的已保存数值为 `d_joint=0.7987161019689071`、`d_cal=3.5605119112114565`、目标位移 `0.1`。其原记录明确是 `ORACLE_CONSTRUCTED_PHYSICAL_WITNESS`、`conditional_fixed_support_geometry=true`、`strict_original_prior_GO=false`、`posterior_probability=null`。因此不能写成 calibration 近等价，不能分配新先验质量，也不能把这个旧点对当成本次宏网格研究的新发现。若独立 calibration 进入统计实验，不能继续只使用 joint 距离。

## 当前实现与证据的分工

`src/wslam/physics.py` 的点输出转置场使用复双线性转置 `T`，不使用复共轭；`point_coefficient_jacobian` 与 `surface_shape_jacobian` 最后才由 `packj` 转成实参数信息。`FEM.vjp` 则实现实内积意义下的共轭伴随。这两种操作服务于不同接口，不能互换。

墙材料、目标位姿、目标内禀参数分别占 `forward.Evaluator` 的 Jacobian 列。当前 `left_right` 列名是 `eps_left,eps_right,pose_x,pose_y,target_contrast,target_loss,target_radius`。墙耗散固定为`.08`的 E1 实验不能说明墙耗散也被联合识别；目标材料导数也不能被称为目标位姿导数。导航平移另由固定材料下的天线坐标有限差分构造，不虚构转角规范自由度。

| 保存证据 | 支持读取的对象 | 不能直接替代的对象 |
|---|---|---|
| `results/e0/qualification/RESULTS.json`、`results/e0/refinement01_executable/RESULTS.json` | 原精度失败及固定修复配置的实际指标 | 任意连续材料、任意几何的统一误差界 |
| `results/e0/four_blocks/RESULTS.json`、`BLOCKS.npz` | 同一 FEM 离散中的四块恒等式、已付费端口准备与 Born 对照 | 连续 Maxwell 认证、加速收益 |
| `results/e1/dynamic/pose*_f*.npz` | 实际点观测、物理 Jacobian、源/接收器坐标与列名 | 广泛先验上的后验或全局可辨识性 |
| `results/e1/dynamic/full_*.npz`、`minimal_*.npz` | A/B、C/R/DC、秩与保留率的原始小矩阵 | 单端口机制向完整孔径歧义的推广 |
| `results/e1/dynamic/PROFILE_*.npz` | 固定墙步长下实际有限端点与局部预测 | 所有优化迭代、真实后验质量或 A0 测试表现 |

E0 的初始失败、修复前成本与失败文件必须继续保留。`results/e0/gate.json` 是父线程给出的有限离散资格记录，本文件只引用其存在和范围，不重新作门判断。当前图表来源及原始量范围见 `figures/actual_e0_e1/FIGURES.json` 和 `SCHEMA_AND_MISSING.json`。

## 原始量缺口与记账

E1 保存了频率/位姿的完整点数据与 Jacobian、导航差分数据、信息小矩阵、最终补偿端点和声明的高斯精度矩阵。混淆模型 NPZ 没有自带干扰参数列名向量，顺序须对照 `tools/e1_dynamic.py`；没有完整优化试探轨迹或每次试探复场；没有单独的 E1 全场文件。现有代表复场来自 E0 初始导数网格，图注已注明。不得把这些缺口用重建值伪装成原始保存数据。

绘图曾在独立进程完成而未采集精确 CPU receipt，形成真实记账缺口。父线程已说明按保守180 CPU秒上界计入 preparation 成本；这不是测得的绘图 CPU 时间，也不能再次收费。本文件不以零成本掩盖该缺口。源包后续提取通过 `tools/run_job.py` 串行预约，记录于 `results/jobs/audit-source-extraction-001/`。

## 主要文献 DOI 线索：未完成当前在线核验

| 研究包/主线程线索 | DOI | 当前核验边界 |
|---|---|---|
| iSAM2 | `10.1177/0278364911430419` | 重试返回HTTP429 |
| IGG/SPO2026 | `10.3390/robotics15050087` | 未包含在四项网络重试中，不能写成429或核验通过 |
| Watson2024 | `10.1088/1361-6420/ad5b83` | 重试返回HTTP429 |
| Wang–Amin2006 | `10.1109/TSP.2006.879325` | 重试返回HTTP429 |
| Yao | `10.1109/TIM.2024.3502875` | 重试返回HTTP429 |

四项网络重试证据为 `research/delegated/literature/rawbundle/verify_network_retry.json`。这些 DOI 是原始文献检索起点，未在本轮成功读取正文、核对版本或核查定理条件；不构成本轮新颖性结论。
