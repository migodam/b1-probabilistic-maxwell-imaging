# 冻结 provider 的风险评价入口

`tools/r4_risk_panel.py` 是可复用的四阶段工具；此次只执行已暴露的DEV128，没有执行LOCKED生成/真值读取。未来LOCKED必须由主线程另行授权，方法/父样本由调用方事先固定，不自动扩展。

1. `freeze --spec SPEC --out FROZEN`：冻结方法ID、模型checkpoint路径/字节数/mtime、provider来源与K、比较和父样本。旧Spline保留K32；DiagK32、FieldK0。已知六个坏header只通过加法说明、原checkpoint身份、封存TRAIN/validation gate恢复，未知坏header拒绝。新版冻结同时封存兼容说明身份。
2. `features --freeze FROZEN --observations ORIGINAL_ONLINE --out FEATURES`：子进程调用既有`r4_fusion_moments.py`，只用原始在线观测。Spline/Diag固定K32；Field确定性K0且禁止伪造x_samples。每provider独立封存来源/每父文件/资源。原观测身份在开始和结束核对。
3. `infer --freeze FROZEN --contexts ONLINE_CONTEXT --native SAVED_DETAILS --online PREDICTIONS --out FEATURES`：只读observation/q0/chi0/J/residual；native_saved只复用已付费原输出，不调用Target。模型参数、provider、K、parent、源gate身份必须匹配。失败明确回退原anchor并保留分母；原anchor缺失则停止，绝不删除父样本。每方法保存raw/原矩形截断动作；所有模型与预测封存后写ONLINE_COMPLETE。最新版入口同时封存全部moment数组身份。
4. `score --freeze FROZEN --online PREDICTIONS --truth OFFLINE_ONLY --out SCORES`：核对ONLINE_COMPLETE、独立PID、冻结/完成先后关系、模型/来源/输出身份，之后才启动独立评分进程打开真值。输出全576复MSE/MAE、ROI、原anchor与truth的host/transition/pure分层。逐父样本seed损失平均后4000配对bootstrap，seed20261025099、区间98.333333%。不把三seed当384独立父样本，不把平均损失误写成图像集成。

共同DEV面板的封存输出为`results/development_risk_strong_panel_v1`。该次完成后才增强未来代码的moment数组/兼容说明身份封存；加法FINAL_AUDIT核对本次全部384个provider数组及所有模型/预测仍与原封存身份一致，没有修改已完成ONLINE_COMPLETE或FROZEN。

`tools/run_r4_dev_risk_panel.py`只能运行固定R4DEV00000..00127；它不是LOCKED启动器。`tools/summarize_r4_risk_panel.py`汇总既有逐父损失及选中validation，不读真值。验证集1024/256与旧512/128单独列示。报告不形成后验校准、Flow独特优势或Maxwell兼容资格判断。

已完成六个provider训练的metadata修复仅改未来trainer的变量命名（moment_count与parameter count分开），不改网络/损失/优化器/数据/任何已关闭checkpoint。六个恢复模型的参数逐项不变，单样本预测与封存64批量预测的float32差异记录在summary_v2；一次过紧2e-6验证失败保留原FAILED账本，重汇总容差1e-5，未重跑训练、特征或DEV推理。
