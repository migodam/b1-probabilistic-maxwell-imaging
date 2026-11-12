# Fresnel 实测留出预测：后验独立通道验证

实际2001 Fresnel单柱和双柱2GHz数据；1053训练、711留出，source-specific receiver maps及空目标二阶多极照明标定。没有逐像素GT，所有材料NRMSE/IoU均保持缺失。

对56个raw/corrected/primary/secondary视图验证，27个不同材料；使用27次LU与972源RHS。不做候选优化，不按留出结果改变primary。

corrected primary的留出relative L2约9.15%（单柱）/9.86%（双柱），各方法接近。DR raw约10.12%/11.17%，本例corrector改善实际预测。15%为模型discrepancy，不是已证明CN99噪声界；接收horn与完整噪声协方差仍有资格限制。

逐项数据：[MEASURED_HOLDOUT.csv](MEASURED_HOLDOUT.csv)、[CORRECTOR_COMPARISON.csv](CORRECTOR_COMPARISON.csv)、[独立物理成本](UNIQUE_VALIDATION_COSTS.csv)、[运行回执](RECEIPT.json)。完整复数测量预测NPZ在Release分析附件中。
