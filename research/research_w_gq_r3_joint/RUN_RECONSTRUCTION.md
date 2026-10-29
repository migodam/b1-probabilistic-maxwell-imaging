# 从测量到二维联合成像

主入口读取测量与公共配置，输出连续椭圆目标材料图、六宏格未知墙体、完整多频预测及所有优化分支。不读取真值，不从旧目标标签表挑选。

## 可移动运行包

`delivery/portable_runtime/` 包含冻结的 R3 科学源码、实际用到的12份历史后端源码与依赖版本。在 Python 3.13 环境安装该目录 `requirements.txt`。本机使用已经安装的科学环境；全新环境安装是否实际执行以最终复现记录为准。

从 `delivery/portable_runtime/` 目录运行：

```bash
PYTHONPATH="$PWD/vendor" python tools/deploy_r3.py \
  --deployment-budget-seconds 1800 \
  --measurement ../../inputs/online/blind_normal_r1/observation.npz \
  --config ../../inputs/online/blind_normal_r1/config.json \
  --output ./my_reconstruction \
  --method FULL_GN
```

输出目录必须不存在或为空。独立运行显式声明自己的时间预算和独立账本，不延长已经结束的三小时研究账。每个进程单线程、内存上限4 GiB。其他主方法名称为 `FIXED_CURRENT_CHART` 和 `ADAPTIVE_CURRENT_CHART`，它们不具有预先保证的准确性或加速优势。

测量 NPZ 只需 `y`：复数数组 `(3,4,16)`，顺序为频率、发射、接收。公共配置固定源与接收位置、频率、每实/虚坐标噪声0.01、材料容器与参数约束。不能把任何外部设备的数据不经归一化直接当作本合同测量。

## 输出

- `reconstructed_target.png`、`target_material.npy`：目标ROI图与 `(y,x,2)` 数组，通道为实介电率和参考频率损耗参数；坐标见 `reconstruction.json`。
- `reconstructed_wall.png`：真实输出数组中的墙体区域图；`endpoint.npz` 保存全材料场与连续参数。
- `predicted_measurement.npz`：完整 Maxwell 复数预测。
- `reconstruction.json`、`SUMMARY.json`：按完整数据目标函数与冻结 prior 选择的分支；不会按真值误差重新选。
- `optimization_trace.csv`、`TRACE.jsonl`、`STEPS.jsonl`：实际候选、接受/拒绝和停止原因；非收敛照实保留。
- `runtime_cost.json`、`deployment_accounting/PHYSICS_LEDGER.jsonl`：此次运行的CPU、内存与物理调用。历史两个计费收尾失败不影响新入口，原缺失值仍保留为缺失。
- `ONLINE_SEALED.json`：在线输出封存标记。它不表示材料正确或已取得全局最优。

## 独立评分与原研究复查

在完整 R3 研究目录、使用相同科学 Python：

```bash
python tools/score.py --output results/online/blind_normal_r0/FULL_GN \
  --case blind_normal --repeat 0
python tools/aggregate.py
python tools/plot.py --all
```

评分器单独读取 `inputs/private/`，且要求在线封存与方法冻结。在线进程的文件访问日志由 `tools/guarded_online.py` 生成并拒绝读取 private/scoring 文件；这是可审计的 Python 进程接口隔离，不宣称操作系统级安全隔离。

原三小时实验入口 `tools/run_r3.py` 受 `RUN_CLOCK.json` 原期限限制，结束后不能直接重启原账。未来复现请用上述显式独立部署入口。主要源码、实际执行与冻结修订的对应关系见 `METHODS_FROZEN.json`、`DEPLOYMENT_INTERFACE.json` 与最终报告。
