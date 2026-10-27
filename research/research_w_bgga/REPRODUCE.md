# 复现与执行边界

本目录位于既有B1仓库 `research/research_w_bgga/`。依赖版本见requirements.txt。默认CPU单线程，Python3.13；记录的Apple M4结果不自动等于其他平台数值资格。

先在独立Python环境安装requirements，再从本目录运行。工具会把`src/`与固定`vendor/`加入子进程路径；本地已有依赖目录仅作library复用，不改旧环境。

```
python tools/run_job.py PREFLIGHT reproduce-tests python -m unittest discover -s tests -v
python tools/run_job.py PREFLIGHT reproduce-mesh python tools/mesh_benchmark.py
python tools/run_job.py PREFLIGHT reproduce-summary python tools/report_preflight.py
```

job label必须唯一，避免覆盖旧收据。测试是矩阵/有限概率/几何积分，不执行PDE。已有原始mesh与quadrature在`results/preflight/FIXED_MESH.npz`，图片从这些真实数组生成。

`tools/preflight.py`用于首次重建提出的合同，不应在已有授权/科学运行后重写freeze。固定vendor源码已随仓库发布，无需再次从main取源码。来源commit列表在CONTRACT.json及source_review/PREFLIGHT.json。

物理入口：

```
python tools/run_job.py G0 g0-initial-qualification python tools/g0_qualify.py
```

**当前这条物理命令会因缺少BGGA独立额度而拒绝，且不会产生FE/RHS。** 只有人类批准新scope后才把实际批准内容登记到AUTHORIZATION.json；不能把建议cap或旧stage余额当成批准。runner拒绝后可保留日志，不通过修改门或绕过wrapper运行。

初始资格驱动不是完整12-parent信息/风险stress，也不是G1/G2。后两阶段方法/速度实验未运行；不提供会偷偷启动它们的“all”命令。旧原始数据路径见`research/delegated/source_audit/CACHE_MAP.json`；在线算法不能直接读其中reference目录。
