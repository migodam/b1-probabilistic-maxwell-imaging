# 解析与历史物理对照

[exact_controls.py](exact_controls.py) 用精确有理数、根式和符号矩阵独立复算 Pro §§1–4。[结果](exact_controls.json) 中已比较以下机制。

| 对照 | 核验结果 |
|---|---|
| 源方向信息 | 对 (H=\operatorname{diag}(10,1)\)，只保留一阶秩得 1，检查 (0\preceq Q\preceq\operatorname{diag}(0,100)) 得到同为 10 的原与对偶端点，单位均为距离平方。 |
| 零负惯性 | (G=i)、(S=B=1) 时，两个合法电流的中点要求材料 (2/5+i/5)，落在原实材料合同外。 |
| 无损两单元补偿 | (u=1)、真峰值 2、同权重旧下界为零。保留中心关联，得到距离平方下界 (127/48)，即距离达到信号范数的 35.0802%。同组权重只使用幅度帽及粗中心扣分时，下界为零。 |
| 源矩阵乘子 | 逐区域标量乘子距离平方至多 (373/200)。矩阵乘子 (t=3) 给出 (88/47)，严格多 (69/9400)。真最小值为 2。 |
| 原惯性保持逐候选的低秩 | 三单元对照为 (d^2=1) 对矩最优值 0。五单元对照在整个实材料盒正演可逆，仍是 (d^2=2317/2400) 对矩最优值 0。点态负项秩至多 1，其平均秩为 2。 |
| 合法有损接收椭球 | 零空间无需分解的接收椭球距离原、对偶端点均为 (1/2)。另一个具有 (\Im\chi/\lvert\chi\rvert^2=1) 的材料产生椭球边界点。 |

这些代数模型保持固定 Green 与材料合同。其矩阵系数独立构造，数值幅度没有外推到 96² TM、任意权重或实测场景。五单元逆向对照只用一片共同区域，与历史四区域证书是不同的证据。

[历史对照回执](historical_controls.json) 另核对两件历史 10×20、2% 四区域见证的零目标及能量帽，以及历史 24² 环形案例 `r14_new_0018` 的严格盒内见证。该件峰值为 4.28189000958，严格端口残差上端为 0.006556961951，小于 (u=5) 合同的原半径 0.006634690037。

复算从 Inv_SLAM 根开始。符号对照需要 SymPy，历史标量核对仅需 Python。四区域回放需 NumPy、SciPy、python-flint 和 threadpoolctl，复用八个已保存候选方向及其均值。独立重放 [区间回执](historical_moment_replay.json) 恢复 8,499.9874 与 6,999.99993 的统一相干能量上端。两个接收精确均值的四区域扇区均严格满足，矩松弛目标为零。原有联合惯性距离下端分别为 0.160699、0.110608，高于配对噪声半径。

```sh
python research/z_alg_optical_inertia/ip_unified/numerics/controls/exact_controls.py
python research/z_alg_optical_inertia/ip_unified/numerics/controls/read_historical.py
base=research/z_alg_optical_inertia
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  python -B "$base/ip_submission_v2/energy_sdp/verify.py" \
  --proposals "$base/ip_submission_v2/energy_sdp/alt025.npz" \
  --output "$base/ip_unified/numerics/controls/historical_moment_replay.json"
```
