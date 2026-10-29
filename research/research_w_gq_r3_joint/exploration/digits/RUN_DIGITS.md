# Digit8 新非参数材料格探索：运行记录

本轮实际重建已完成。18:15以后不启动新重job；不得覆盖既有私有truth、端点或文件访问日志。本文件记录真实执行过程，不授权越过原18:30实验截止重跑。

## 结果身份

- `digit_patterns.py`、`digit_1..9.npz/png`：九个公共字形定义，**不是九个成像结果**。
- `private/truth.npz`：只为数字8生成的独立材料参数、同离散全波clean数据；online禁止读取。
- `online/observation.npz`：只有 `y`，形状3频×4TX×16RX。
- `results/digit8/endpoint.npz`：实际110实变量自由材料反演端点，不是ellipse或标签预测。
- `results/digit8/ONLINE_SEALED.json`：全部在线输出封存后才允许评分器打开truth。
- `SUMMARY.json` 的 reconstructed IDs 表示存在真实逆问题输出，**不表示成功辨认/恢复**；本例没有成功恢复8的声明。

## 原始命令（从R3研究目录）

```sh
DIGIT_PYTHON="/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/A25/a25-adaptive-om-feedback/.venv_a25/bin/python"
"$DIGIT_PYTHON" tools/run_suite.py digit8_generate PYTHON tools/digit_generate.py
"$DIGIT_PYTHON" tools/run_suite.py digit8_online PYTHON tools/guarded_online.py --audit-log exploration/digits/results/digit8/FILE_ACCESS.jsonl tools/digit_stress.py
"$DIGIT_PYTHON" tools/run_suite.py digit8_score PYTHON tools/digit_score.py
```

`run_suite.py` 在数值库导入前设单BLAS线程、监控4GiB与原活动截止，并记录包含imports的子进程CPU/实际wall。模型构建及每次全波/梯度动作共用R3 `.fem.lock`。新数字模块另有累计450FE/3600RHS与18:30硬停，先到者为准。需要未来新实验时须另立独立目录、合同与时钟；现有生成器会拒绝覆盖truth，读guard会拒绝覆盖日志。

## 模型与范围

墙为原3×2全矩形公共格，6个epsilon与6个loss未知。目标为7×7、单元0.2λ₀，全部49个contrast和49个loss独立未知，共110实变量。5列×7行的数字8置于中间5列，另外两列也在online未知量内；不给笔画support、digit ID、孔洞标签、模板或ellipse。笔画宽0.2λ₀，总字形1.0×1.4λ₀，真孔各约0.6×0.4λ₀。该对象明显大于原小ellipse案例，不能称同尺度分辨率扩展。

P3 FEM使用公共targetstep0.1λ₀、h0.17、原PML/源/接收器，真实材料相关incident弱载荷完整保留。每次完整objective+weighted-adjoint gradient为3FE/24RHS，只有一个合法mixed wall/target方向差分核对。生成与online使用相同离散，不是独立细网格确认。

目标为独立uniform box，无空间正则、无TV、无形状prior。墙使用原proper Gaussian约束prior。两个公开均匀墙初始化、全目标contrast0.1/loss0.03；L-BFGS-B各40iter、65真实calls上限，累计包括生成与核对不超过140calls。实际两start均到40iter限制，不能称充分收敛、算法最优上限或物理无解。

## 展示与评分

读取新truth前已保存 `DISPLAY_SCORING_RULE.json`：连续contrast同色轴0..8；loss同色轴0..0.25。support主例固定χ≥1，同时报告χ≥0.5/1/2全部IoU，不挑最好阈值、不新增PASS门。报告材料相对L2、contrast质心、空白外溢、孔洞连通、完整residual和停止原因。`truth_vs_reconstruction.png/pdf` 已实际打开目视检查，重建不清楚构成数字8，两个真孔未可靠保留。

## 只依赖测量与公开配置的独立online入口

以下入口不读取private、字形模块或评分结果；生成/评分不属于online依赖。当前实验时钟已经冻结，命令仅供同一授权窗口记录或未来另立合同时复用，禁止绕过截止强行重跑。

```sh
"$DIGIT_PYTHON" tools/run_suite.py digit8_online_new PYTHON tools/guarded_online.py --audit-log NEW_OUTPUT/FILE_ACCESS.jsonl tools/digit_stress.py --measurement exploration/digits/online/observation.npz --config exploration/digits/online/config.json --output NEW_OUTPUT
```

私有generator是 `tools/digit_generate.py`，只有离线scorer `tools/digit_score.py` 会在seal后读真值。字形模块只能由generator/scorer读取，online打开日志确认未读取。

## 标准输出与审计

实际已保存 `results/digit8/{predicted_measurement.npz,target_material.npy,reconstructed_target.png,reconstructed_wall.png,runtime_cost.json,endpoint.npz,TRACE.jsonl,STARTS.json,ONLINE_SEALED.json}`。`target_material.npy`为7×7×2：实部ε=1+contrast、loss系数ℓ。`tools/digit_artifacts.py` 是纯封存数组/公开输入的别名与审计步骤，不读truth/glyph、不运行FE。

`AUDIT.json` 的 `pipeline_dataflow_audit_passed` 只表示数据流/成本/时序通过，**绝不表示数字8成像成功**。两start达到40iter限制的原始停止原因完整保留。探索子账本与R3主PHYSICS_LEDGER记录的是同一物理动作，不能重复相加。
