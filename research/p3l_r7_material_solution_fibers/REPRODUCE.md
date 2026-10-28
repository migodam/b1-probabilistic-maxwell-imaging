# 可复现入口

本轮已封账，`CLOSURE.json`阻止在同一账本继续求解。**不要删除closure、清空ledger或覆盖raw来复跑。**

## 零物理复核

使用原Gaussian研究venv（NumPy/SciPy/Matplotlib/threadpoolctl），从独立R7 worktree执行：

```bash
PYTHONDONTWRITEBYTECODE=1 ../../.venv_nn/bin/python research/p3l_r7_material_solution_fibers/audit/results/audit_results.py --output /tmp/r7_independent_array_review
PYTHONDONTWRITEBYTECODE=1 ../../.venv_nn/bin/python research/p3l_r7_material_solution_fibers/src/theory_controls.py
```

第一条读取无损证据中的固定数组与中央账本，独立重算128状态、五节点、成本及封存身份，不读取GT。第二条仅是解析toy。数学/decoder QA脚本位于audit及src/validate_decoder_audit.py，历史原始输出不要覆盖，复核请写独立输出目录或使用包内复核收据。

## 原始执行顺序

E1冻结256点→decoder QA→E2物理probe→4DEV→保留错误并纠正→LOCKED_FINAL_FREEZE→generate8→parent逐个单worker→seal→单独offline score→注册DEV profile/gauge控制→close。真实命令和日志已保留。生成与online共用原forward，但在线call graph没有offline_truth读取；score是唯一正式GT评价入口。

R7的src/physics.py使用原R4/P1与A30 frozen vendor；原目标配置、geometry、background及可重建源码作为必要复用输入单独无损归档，见ASSET_MANIFEST。新物理复跑必须在独立目录与明确新预算下进行，不能误改本轮封账或历史R6。单物理worker、BLAS1、最多一个live LU、complex128/float64与原packing必须保持。

原始数组不放普通Git blob；Release完整证据恢复后才可用np.load。归档保留原文件bytes与相对路径，BLAKE2b/逐字节readback收据与代码身份随包提供；不做SHA-256附加检查。网页和Markdown只是阅读层，不能替代原始证据。
