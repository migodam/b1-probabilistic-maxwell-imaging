# 可复现入口

1. 从交付ZIP解压本目录；保留 src、vendor、tools、configs、inputs、results、receipts。
2. 使用 requirements.txt 中的依赖版本；原执行环境精确信息见 ENVIRONMENT.json。
3. 实际新部署使用 RUN_ALGORITHM.md 的单命令入口，明确给出自己的本地时长并输出到新目录。
4. 无需PDE的复查：`python tools/score_all.py --results results`，再运行 `python tools/summarize_cost.py`。
5. 原理对照入口：`tools/absorption_controls.py`、`tools/spectral_pilot.py`、`tools/independent_validation.py`。它们会执行新物理调用，受原RUN_CLOCK和新额度隔离，不应把已结束额度复位。
6. 后处理图均来自交付真实数组；主物理格7×7，所有主图使用相同contrast/loss色轴。封存点以完整测量目标加冻结先验选择。

原始观测只含y；真值在 inputs/private 下，只属于生成和封存后的独立评分。源码与配置中的相同公式和数组顺序，不表示一般连续PDE误差认证。与同离散FEM的精确一致、独立更细网格经验误差和实际目标恢复指标分别报告。

本轮不自动提交或发布。原R2/R3代码、账本和裁决均未修改。

## 已执行与未执行的复现检查

实际批命令已在当前环境完成全部B0/B4数字运行与Full-GN补充，所有原始输入、选中终点和停止状态均保存。`vendor/`保留所需12个历史源码模块，不需要旧2950原子bank。`tools/deploy.py`的新预算封装通过静态路由检查；全新虚拟环境重新安装、从新部署命令再次重复求解均为NOT_RUN，未冒称已验证。闭账后原预算拒绝新物理；读取原始数组、重新评分和作图不调用PDE。

## 不重算Maxwell的复查

```bash
python tools/score_all.py --results results
python tools/summarize_cost.py
python tools/finalize_delivery.py
```

这会刷新派生报告，保留封存endpoint和原始物理账本。请先阅读 `COST_ACCOUNTING_SCOPE.md`：旧Book计时未覆盖Python启动/导入，后续whole-child计时单列且不得重复加入。原时间戳、绝对路径和机器版本是出处证据；实际运行命令中的输入输出可改为本地相对路径。
