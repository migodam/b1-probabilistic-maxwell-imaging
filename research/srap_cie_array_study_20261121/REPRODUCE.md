# 复现实验

二维使用 Python、NumPy、SciPy、matplotlib、threadpoolctl。请从主仓库检出本目录及相关 R2、R4、R6、R8、R10、R12、R15、FBG-R3 和 Z-ALG-R4 代码。在 Inv_SLAM 工作区也可直接复用历史工作树。三维使用带 CUDA 的 PyTorch，Windows/i9 与 RTX 4060 的实测成本见三维报告。每次新研究应在空输出位置执行；现有同名预测表示已经封存，脚本会跳过该臂。

所有下列命令均在本研究目录执行。单线程环境应在启动 Python 前生效，脚本也会设置常用 BLAS 环境变量。实测输入取 Fresnel 2005 原始七列 TM 文件 FoamTwinDielTM.exp、FoamDielIntTM.exp 和 FoamDielExtTM.exp，六个旧设置及两个匹配掩码共享原始文件。[原实验来源](https://doi.org/10.1088/0266-5611/21/6/S09)

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export FRESNEL_2005_DIR=/path/to/Fresnel2005
export SRAP_HIST="$PWD/legacy_inputs"
python prepare_legacy.py
python data2d.py
python data_followup.py
python run_study.py --cohort legacy --arms dual,u11 --workers 6
python run_study.py --cohort legacy --arms remaining_legacy --workers 6
python run_study.py --cohort development --arms dual,srap,dual_cap,srap_cap --workers 2
python run_study.py --cohort confirmation --arms dual,srap --workers 4
python run_study.py --cohort measured_extra --arms dual,srap --workers 2
SRAP_DATA="$PWD/data_followup" python run_study.py --cohort followup --arms followup_all --workers 4
```

两个直线阵的接收谱放大很大，主评分还需要两条截断臂。可在本机一次完成；正式运行把 sorted(object_id) 的零基索引模3余数2分给Mac，余数0、1分给XINAN，所有方法在对象内配对。脚本接受 --partition 2/3 与 --partition 0,1/3。跨机器复现能固定物理输入和迭代，墙钟前缀会随硬件、Python/BLAS版本和并行负载变化。

```bash
ids=$(python -c 'from run_study import cases; print(",".join(c["case"] for c in cases() if c["cohort"]=="confirmation" and c["array"].startswith("LS_LINE_")))')
SRAP_OUTPUT="$PWD/results/line_cap" python run_study.py --cohort confirmation --ids "$ids" --arms dual_cap,srap_cap --workers 4
ids=$(python -c 'from run_study import cases; print(",".join(c["case"] for c in cases() if c["cohort"]=="confirmation" and c["array"].startswith("LS_ARC_")))')
SRAP_OUTPUT="$PWD/results/arc_cap" python run_study.py --cohort confirmation --ids "$ids" --arms dual_cap,srap_cap --workers 4
ids=$(python -c 'from run_study import cases; print(",".join(c["case"] for c in cases() if c["suite"]=="LA_FRESNEL"))')
SRAP_OUTPUT="$PWD/results/lazy_fresnel" python run_study.py --cohort legacy,measured_extra --ids "$ids" --arms dual,srap,srap_lazy --workers 2
ids=$(python -c 'from run_study import legacy_cases; print(",".join(c["case"] for c in legacy_cases() if c["suite"]=="LA_SYNTH"))')
SRAP_OUTPUT="$PWD/results/legacy_mac" python run_study.py --cohort legacy --ids "$ids" --arms dual,srap,u11 --workers 2
```

评分按 cohort 检查预测齐全才打开真值或真实留出场；公开开发干预可以独立分析。所有数组为 NumPy 无 pickle 格式。每个二维预测保存 raw / light 的两个预算材料、按时间登记的物理目标、材料轨迹、路由诊断、计数和候选资格。公布的图与表从封存结果产生。

```bash
python score_study.py --cohort development
python score_study.py --cohort legacy
python score_study.py --cohort measured_extra
python score_study.py --cohort confirmation
SRAP_DATA="$PWD/data_followup" python score_study.py --cohort followup
python score_replay.py
python score_lazy.py
python aperture_noise.py
python mechanism_stats.py
python paired_stats.py
python write_tables.py
```

三维在独立目录内运行。实际 GPU 脚本默认使用 CUDA 设备；空间单位为背景波长。归档内保留 V1 共面圆环的开发输入、预测和实际费用，正式比较读取 plane_grid、hemi。几何、极化、物理误差检查及计时口径见 vector/FREEZE.json。图与统计使用 CPU 或 GPU 的说明见 vector 报告。

```bash
cd vector
python generate.py
python qualify_finer.py
python invert.py
python score.py
python report.py
```

二维逐区域材料误差和读图用的同色标图在原始输出封存后计算。若要重画全部全景图，运行下列脚本；二维统计按物体而非按通道或曲线点重采样。

```bash
python spatial_balance.py
python cost_summary.py
python visual_contact.py followup PW180_10x20
python visual_contact.py confirmation PW180_10x20
python figures2d.py --section schematics legacy development measured_extra fresnel_lazy legacy_cross_host confirmation
SRAP_DATA="$PWD/data_followup" python figures2d.py --section followup
```


`python correction_effect.py` 从已经评分的轨迹读取相同付费种子，并输出每个共同校正前后的原物理目标和材料、支撑指标。30 秒 raw 终点与其 24 秒种子有时不同，直接相减不能隔离校正作用。

### 从已算完的预测重新分析

[同仓库 Release 的单个数组附件](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/srap-cie-array-study-20261121) 包含二维合成公开输入、离线真值、各预算材料与完整付费轨迹、开发干预场，以及三维数据、材料和运行记录。解压目标为本研究目录，保留相对路径。读取封存材料或重画合成图无需重新做耗时反演。真实数据所需的原始测量与测量输入包由上列 Fresnel 附件取得和校准；归档含反演后的实测材料及留出预测误差。
