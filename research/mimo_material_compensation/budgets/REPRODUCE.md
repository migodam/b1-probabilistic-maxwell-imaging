# 复算有限面积外区预算

在 Inv_SLAM 根目录运行。需要 Python 3.13、NumPy、SciPy、Matplotlib、python-flint 和 threadpoolctl。所需旧项目包含读取允许输入的 Fourier 算子、加载网络、噪声测量和公开参考候选。

freeze.json 保存粗扫合同，task_amplitudes 和 --focused 是已观察阈值后的同对象开发测试。所有脚本只在本目录写入。

    cd research/mimo_material_compensation/budgets
    python3 -B quick.py
    python3 -B bound.py joint32
    python3 -B bound.py air32 air96
    python3 -B bound.py native96
    python3 -B bound.py --focused
    python3 -B shape_mean.py
    python3 -B shape_witness.py
    python3 -B conic_reference.py
    python3 -B replay.py
    python3 -B figure.py
    python3 -B bridge/bridge.py
    python3 -B bridge/replay_bridge.py

raw_*.npz 保存每个锚点的材料、复数电流、接收与发射总场、前向响应，以及实数化方向候选。前向场和导数的实数化顺序为全部实部、再全部虚部，矩阵内接收端为第一维、发射端为第二维。shape_* 保存平均值的支撑方向和形状反例。two_endpoint_raw.npz 保存小型 SOCP 的实化系数和乘子。independent_replay.npz 保存复核后的材料、裸 VIE 电流和加载端口复数矩阵。JSON 索引指向实际数组键。

budget_phase.pdf 可直接用 LaTeX 的 \includegraphics 引入，table.tex 给出关键行。所有下界小于阈值和无内见证的条目保留为未决。名义固定仪器下的点态证据与宽先验结论使用不同量词。
