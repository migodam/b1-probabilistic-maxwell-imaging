# 完整成本与可比较范围

| 新增工作 | states | complex RHS | LU | 受控wall |
|---|---:|---:|---:|---:|
| F0 decoder/cached-J代数 | 0 | 0 | 0 | 单列GEOMETRY_TESTS |
| F1 ring | 12 | 480 | 12 | 2.725971s（条件内部） |
| F1 weak | 12 | 480 | 12 | 2.752128s（条件内部） |
| F1合计含共享准备/持久化 | 24 | 960 | 24 | 5.534289500s |
| F2四臂 | 0 | 0 | 0 | NOT_RUN |

每个新state为8个源RHS，共192；每个完整receiver transfer为32个adjoint RHS，共768。没有用低维latent维数偷换物理RHS。共计960，单次full-state只计一次state，不将其伴随算作免费。所有失败和拒绝本应付费，本轮24次计算均成功。

旧背景G/S/incident、background解和旧Jacobian属于已支付资产；本轮读取和方程矩阵审查的时间计入执行，历史生成成本没有重记新增。24份新LU/Jacobian从零求取并计费；不是免费缓存部署速度。没有做cold/warm推断吞吐基准，5.534秒是本机受控小实验总wall，不代表一般Maxwell反演时延。

总授权64states/2048RHS/1200s；未用40states/1088RHS不因余额存在而开启Gate2失败的F2。不与已关闭45/1880诊断、LXZ96/896、B1旧384合并。没有A17/ROM/full-Krylov的新同成本消融，不能宣称其正负加速。
