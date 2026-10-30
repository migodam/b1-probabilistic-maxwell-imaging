# 从实际复测量到未知墙、未知目标材料图

依赖为 `requirements.txt` 中实际环境版本。主实现为二维 scalar-TM P3 FEM，Apple Silicon CPU，单线程。`vendor/` 包含本算法需要的12个历史 Maxwell 后端源码模块；不需要旧 bank、R3结果、真墙或数字模板。全新依赖安装是否测试以最终 `REPRODUCE.md` 为准。

## 独立新部署：一条命令

在本目录执行：

```bash
python tools/deploy.py \
  --measurement inputs/online/c001/observation.npz \
  --config inputs/online/c001/config.json \
  --options configs/main.json \
  --method B0 \
  --deployment-budget-seconds 900 \
  --output my_new_reconstruction
```

`B0` 为完整 Maxwell 联合 FWI+TV；`B1` 为原生固定GS的 contrast-source CSI；`B2` 增加精确吸收及定时重新开放；`B3` 在B2上增加内部模态 Schur 预条件；`B4` 采用内部模态 Schur 预条件和经验吸收/重新开放。`--deployment-budget-seconds` 明确创建与已结束研究独立的新本地计算额度，不修改原 `RUN_CLOCK.json`。输出目录必须为空。

输入 NPZ **仅含复数组 y**，排序 `[frequency,source,receiver]`，主实验为3×4×16。每个实虚坐标标准差0.01。在线程序读取公共墙ROI、目标ROI和材料box；不读取 private truth、字形、真实墙/目标参数或评分目录。

输出包括 `target_material.npy`（f0时chi+i loss）、`wall_material.npy`（epsilon+i loss）、`predicted_measurement.npz`、二维材料图、`reconstruction.json`、`optimization_trace.csv`、`runtime_cost.json`、各起点与各阶段停止记录。物理重建图显示7×7材料格，不把插值当新增分辨率。

注意：上面的新部署输出是复数二维Npy；原研究B0/Full-GN的 `complete_outputs` 别名则为 `(ny,nx,2)` 的 `[epsilon_inf,ell]`。不要混读；跨方法统一读取 `endpoint.npz` 的 `target_contrast`、`target_loss`、`wall_eps`、`wall_loss`。所有ell均为色散系数，任意频率虚介电率是ell/f。细节见 `DATA_FORMATS.json`。

## 原研究批量入口

```bash
python tools/run_digits.py \
 --online-root inputs/online \
 --cases c001 c002 c003 c004 c005 c006 c007 c008 c009 \
 --options configs/main.json \
 --output results/B0_FULL_FWI
```

这组通用编号对应九份独立测量。在线调用不含数字ID；ID与真值的对应关系只在独立生成/评分文件中。该入口遵守原研究时钟，截止后拒绝新增物理调用；复现新的实际求解应使用上面的显式新部署入口。

## 独立评分

```bash
python tools/score_all.py --results results
```

评分只接受带 `ONLINE_SEALED.json` 的最终材料输出，读取真值的时间在封存之后。所有数字使用相同阈值chi>=1；同时报告连续材料误差、质心、背景伪材料、4/8连通性及4连通背景孔洞。主范数不除通道数平方根，描述性RMS另列。留出频率0.925的预测由独立物理调用产生并单独计费。
