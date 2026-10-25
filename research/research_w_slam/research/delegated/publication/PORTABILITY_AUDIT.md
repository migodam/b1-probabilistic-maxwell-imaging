# 公开副本可移植性核查

范围：静态读取发布helper、复现说明、imports和保存路径；未执行导入、安装、打包、网络或物理求解，未修改活跃bank、核心科学模块或其wrapper。以下是路径与包装核查，不是新的科学门或跨平台数值认证。

| 组件 | 本次发现 | 处理与仍需注意 |
|---|---|---|
| `src/wslam/` | 只依赖本campaign相对模块及NumPy/SciPy/skfem；`forward.public_mesh`延迟导入Triangle | 不需要旧`wq/wr/B1` Python模块或旧墙缓存；PDE未运行 |
| `tools/run_job.py` | 开发源含相邻`research-w-continuation/.deps`的优先路径 | 开发源保持原样；发布helper只改公开副本，使用stage `src/vendor`与环境安装包，并记录差异 |
| 普通tools兄弟导入 | `a0_wide→a0_run/a0_score`、`a0_adjoint_polish→a0_intrinsic`、绘图/摘要→`polish_artifacts` | 必须保留整个`tools/`目录；按`python tools/name.py`方式运行，不把单个脚本单独搬走 |
| 旧一维负对照 | `e3_negative_controls --analytic-source`默认本stage `legacy_inputs/p0a/analytic_witness.py` | 文件已vendored；通过AST只提取指定函数，不能直接import其带绘图/目录副作用的整个legacy脚本 |
| 依赖版本 | requirements记录旧continuation版本；原job解释器路径属于历史来源 | 未做全新安装可用性或跨平台复算；wrapper是Unix接口，Windows未验证；PDF工具另需Pandoc/XeLaTeX/CJK字体 |
| 保存路径 | ROOT由源文件位置确定；历史JSON常含旧绝对路径；输出多为exclusive目录 | Release须保持stage相对路径；改变job标签不改变科学输出目录。禁止用缺失缓存触发未声明的新求解 |

## 已实施的包装层修复

`tools/prepare_publication.py`原先未复制`inputs`中的JSON、原共享账本、顶层results摘要、job receipt/log，以及有限bank的checkpoint/最终费用文本。现补齐这些文本来源；数组仍不混入Git小包。`--partial`保留其运行中语义，最终有限阶段文本仅在非partial路径纳入。未改写过去的科学JSON、预算或失败记录。

公开wrapper副本移除旧相邻`.deps`，源工作目录和正在执行的bank不变。`PUBLICATION_FILES.json`记录公开副本差异并使用实际目标字节数；README明确来自START_HERE。额外生成`RELEASE_LAYOUT.json`，只枚举NPZ/NPY和未进Git的原始图片的相对路径及大小，**不创建或发布归档**。公开准备若由`run_job`包裹，内部时间仅作信息事件，避免子进程CPU再重复记账；直接调用时保留自身进程CPU记账。

## 必须保持的公开布局

```text
<existing repository>/research/research_w_slam/
  src/wslam/                  # 本地科学模块
  tools/                      # 含兄弟脚本与polish_artifacts.py
  requirements.txt
  legacy_inputs/p0a/analytic_witness.py
  inputs/online/*.npz          # 来自额外Release；JSON小文件在源码包
  inputs/offline/*             # 合成离线标签/均值，禁止在线拟合加载
  results/e0/, e1/, a0/        # JSON/CSV在小包；NPZ保持路径来自Release
  results/e2e3/finite_bank/online/
    FINITE_GRID.npz, FINITE_BANK.npz, OBSERVATIONS.npz
    JOB_COST_FINAL.json        # complete=true后才是完整wrapper准备成本
  results/e2e3/inference/POSTERIORS.npz
  results/jobs/*/{RECEIPT.json,OUTPUT.log}
  results/budget_ledger.json
  PUBLICATION_FILES.json, RELEASE_LAYOUT.json
```

数组归档成员采用stage相对路径，解包位置是`research/research_w_slam/`，不多套一层stage目录。物理bank的paid_rows、代表复场、阶段端点和后验权重按实际layout保留，不以伪造空文件替代。helper当前只生成布局清单；Release上传、匿名取回和完整新环境缓存读回均仍未由本次审计执行。

最终发布前应由主线程基于终态receipt重建清单，确认缓存与图表所需文件实际包含在Release；运行中快照不能代替该事实检查。本审计未增加hash、备份、科学门或新实验。
