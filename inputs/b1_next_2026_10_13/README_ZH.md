# B1 下一阶段研究包

日期：2026-10-13。

## 先读什么

`B1_NEXT_RESEARCH_REPORT_ZH.md` 是完整中文研究报告，包含B1审计、统一材料后验定义、定理/证明/反例、Maxwell与nuisance联系、文献差距、算法比较、阶段实验、预算和论文定位。

`CODEX_PROMPT_ZH.md` 是独立可复制的执行Prompt。它先接入现有B1，默认不额外授权新的物理计算；已有明确授权阶段可按其原预算继续。建议的新P1预算是4,000个材料-频率状态，需要显式采纳，不重置旧384账本。

`verify_theory.py` 和 `checks.json` 是本轮实际执行的代数、反例及已保存数字核验。运行方式：

```bash
python verify_theory.py --output checks.json
```

需要NumPy与SciPy；不需要GPU，不导入B1物理solver，不会调用Maxwell前向。`EXECUTION_RECEIPT.json` 记录执行范围和未完成验证。

## 核心裁决

继续“统一材料空间中的nuisance-aware、event-aware Maxwell posterior inference”。拓扑是生成后计算的属性和输出事件，不是必须预设的孔洞数量专家。未知半径是低成本验证台，不是最终架构。

当前最重要的新增审计项：在同一强信号case中，原生与经验细化的拓扑概率几乎相同，但材料方差相差约227倍。必须分别验证事件质量和连续材料不确定性。

## 证据边界

本包本轮新增Maxwell调用0、GPU训练0、用户仓库写入0。报告中原B1物理结果来自固定commit，未对原仓库所有二进制bank独立重算。小矩阵测试不能替代真实Maxwell验证；现有基础恒等式不被申报为新发明。完整源清单和未解决问题在报告末尾。
补充：原生强信号case的均值与5%/95%分位数，还要求方差至少3.1851e-4，已报方差仅3.2291e-6，违反必要下界约98.64倍。这暴露的是两套近似summary不相容，不是Bayes evidence公式归一化错误；E0已明确加入coherence gate。
