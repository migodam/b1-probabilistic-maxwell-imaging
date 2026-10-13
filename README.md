# B1 — Probabilistic Multi-Hypothesis Maxwell Imaging

公开仓库：[https://github.com/migodam/b1-probabilistic-maxwell-imaging](https://github.com/migodam/b1-probabilistic-maxwell-imaging)。重建图网页：[Top-K imaging dashboard](https://migodam.github.io/b1-probabilistic-maxwell-imaging/)。网页版GPT阅读入口：[START_HERE.md](START_HERE.md)。

中文完整报告：[RESEARCH_REPORT_ZH.md](RESEARCH_REPORT_ZH.md)。图像与交互入口：[dashboard.html](dashboard.html)。结论：[DECISION.md](DECISION.md)。理论：[THEORY.md](THEORY.md)。数字：[RESULTS.md](RESULTS.md)。逐轮记录：[RESEARCH_LOG.md](RESEARCH_LOG.md)。

本项目已实现“16个复散射通道 → 连续材料/拓扑后验 → Top-K事件质量与代表性复材料切面”。物理是3D vector-Mie；材料模型只覆盖三个固定同心模板。竞争案例：单中心空腔50.64%，中心+球壳空隙40.43%，剩余实心8.93%。单张精确图像质量为0；百分比对应整个连续拓扑分支。

**GO：受限三维Top-K原型。HOLD：一般几何、全域误差认证与全面统计校准。未建立强TAP理论/算法贡献。** 原生full-grid部分reference失败保留；第二轮经验refinement不是严格证书。

## 使用现有已付费物理结果

依赖Python、NumPy、SciPy；成图另需Matplotlib。从项目根目录安装依赖后运行：

```bash
python3 -m pip install -r requirements.txt
b1_python=python3
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 "$b1_python" -m src.b1 infer --observation results/cases/competitive/posterior.npz --output results/new_imaging_demo
```

从项目根目录执行。输出目录必须未存在，以保留之前预测。该入口已经实际测试，输出probabilities.json、posterior.npz、实/虚部top_k_material.png与occupancy/SD图；新Maxwell调用0。`--model full-grid`使用原生已付费7阶网格，并对5阶比较；默认`refined`使用已核验的polynomial forward和17/25积分，并明确近似target。

外部观测可以用JSON：`data`为32个实数（16复通道全部实部后全部虚部），`noise_real_sd`为正且固定的σ，可选`topology_prior`为3个正数并归一化，可选`gain_sd`为g~N(1,s_g²)的标准差。必须与protocol.json的角度、偏振、k、单位及物理模型相同；该原型不能直接解释任意实验仪器的数据。程序拒绝改动既付费bank的材料prior/几何/频率配置后重新解读缓存。NPZ同样保留gain_sd，修复前的错误重放已归档。

## 检查与独立复现

```bash
"$b1_python" -m unittest discover -s tests -v
"$b1_python" -m src.b1 summarize
```

当前有24项代数/概率/缓存/预算/元数据测试、19项物理forward审计、4次已有treams独立层球检查。测试不等于科研门槛或网页用户验收。HTML现场浏览器验收未执行：工具拒绝本地file URL；科学图像已由原生Codex subagent逐张视觉核对。

累计新物理预算已用384/384，且账本有不可降低的high-water记录。**不要在当前目录再次运行完整物理campaign，也不要删除/改名账本重置预算。** 在明确获得独立复现预算的新checkout里，源实现的入口依次是：

```bash
"$b1_python" -m src.b1 audit
"$b1_python" -m src.b1 pilot
"$b1_python" -m src.b1 refine
```

实际primary solver不需要treams；refine阶段的4次独立参考复用本机A6/vendor ABI，在独立进程运行。可用B1_TREAMS_VENDOR显式指定已有兼容NumPy/SciPy/treams目录；不修改主环境。A13角向/均匀系数函数已经以未修改摘录随vendor保留，来源见SOURCE_MANIFEST.json。无需GPU/NN/外部数据。

## 文件与边界

configs/protocol.json冻结实验分布；results/PHYSICS_LEDGER.jsonl对每次物理求解先BEGIN、后END，失败/中断亦付费；JOB_LEDGER记录实际host job的wall/CPU；RESOURCE_SUMMARY为最终384而非第一轮356。results/history保留修复前文件；IMPLEMENTATION_FAILURES及元数据修复说明保留所有已知实现错误。

input source包来自用户A32附件，命名仅作为来源；项目统一称B1。A17/A30冻结代码、论文与旧负结果保持不变；没有读取旧locked/truth数组、训练NN或开启GPU。旧包74次2D解析调用与新384次3D调用并列追踪，历史与B1合计458次；A31微型物理另列，未在本项目重跑。

初期机械source/绘图工作曾按旧AGENTS规则使用DeepSeek；收到用户禁止指令后全部终止，后续三个native Codex subagents完成review、文献与视觉交付。后续不再调用DeepSeek，已记录此偏好。

项目现按用户授权公开发布。HTML可在线或离线查看；之前本地浏览器现场运行受工具限制，原始NOT_RUN记录保留。大型真实三维成像、未知位置/边界和可变维多物体后验均未实现，不用模板图冒充它们。
