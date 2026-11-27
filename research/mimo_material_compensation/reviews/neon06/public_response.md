# 单次相干采集的独立审查

**核对结论（口径无额外 8 倍增益漏洞，模拟器未误算单接收链）：**

1. **√8 依据与收发功率**：[coherent_modes.py](../../environment/coherent_modes.py) 乘 $\sqrt{8}$ 是抵消八列基线每列 $1/\sqrt{8}$ 的总能分配（见 [fluence_loop.py](../../environment/fluence_loop.py)），把总能集中至单发（[proposal.json](../../environment/fluence_modes/n96/proposal.json) 发射能 $\le 0.140625$），接收端 $\|\text{rx}\|_2\le 1$ 属无源合成；[independent_check.json](../../environment/fluence_modes/n96/independent_check.json) 残差 $0.8000r$ 与注入噪声吻合。
2. **原创性界限**：[classes_96.json](../../environment/fluence_modes/classes_96.json) 显示 `class_1` 下界 `3.034` 仅胜单对 Bouchet 首波前 `center_1`（`3.014`）`0.67%`，封存幅值同为 `0.375`，收益几乎全来自单发相干聚焦。
3. **最严重实测缺口**：仅 2 例无色散同求解器仿真、零硬件实测；“开路探针+理想无损合成器”忽略了物理负载失配/插损与 8 通道幅相校准漂移，且将 64 元矩阵同一绝对噪声球 $r$ 直接套用于单路合成读数，未检验真实多链热噪声与器件非理想性。
