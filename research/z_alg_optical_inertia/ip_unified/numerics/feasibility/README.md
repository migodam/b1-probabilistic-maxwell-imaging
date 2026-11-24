# 原生 96² 全盒物理残差上界

**物理主张。** 对公开的六个 400 MHz、2 m 见方全域 TM 对象，是否存在一个所有 32 束照明共享的实非负材料，峰值不超过指定 (u\in\{1,1.25,1.5,2,3\}\)，并使 32×32 复端口 Frobenius 残差落入原 2% 误差球？任何已核验盒内材料都为全盒最小残差提供严格上端。只有该上端小于或等于原半径，才得到全盒可行见证。某条路径失手仍留下全盒可行性未定。

## 冻结输入和选择

使用 `physical_revision/array/observations/manifest.json` 内的 `optfresh_01` 至 `optfresh_06`，完整观测 ID 后缀 `__pw360_32x32__pair2pct`。使用同目录 `arrays.npz` 的 `B32` 与 `S32`，Green 算子来自 `synthetic/public/fine_geometry.npz` 的保存 Fourier 符号。32 发乘 `sqrt(10/32)` 保持总 fluence，32 收乘 `sqrt(20/32)` 保持共同接收权重。`padded_radius` 取冻结 manifest，2% 分母取原干净复端口 Frobenius 范数。搜索函数只访问这些冻结公开文件。没有读入生成材料、支撑、干净端口或另一阵列的历史拟合。

`search_native.py` 用两个固定初值，零材料与宽泛高斯，做 12² 起步、24² 细化的有界解析 Gauss–Newton，再开放全部 9216 像素做原生有界 L-BFGS-B。所有预测、接受与选择均用原 96² 模型。`projected_starts.py` 也试低成本投影建议，包含三个起点或继续方向，再用原生 VIE 检查。投影场的残差只用于该分支自身的优化；原生场决定候选排序。`cie_starts.py` 复用历史 Dual-CIE 生成种子，对每个受测上限逐像素盒投影后才计算原生响应，并经全盒校正。CIE 的辅助材料范围见该臂回执。预算内最小残差与未成功的步骤保留在逐案 JSON。

最终材料和 32 源电流保存为 `.npz`。每个选中候选另一次浮点正演用相对容差 `1e-10`。这给出 `numerical_residual`，状态量舍入误差另计。`enclose_native.py` 用 128 位 Arb 的独立 DFT 在原网格重做状态缺陷，并用候选特定的实零真空行简化活动块。密集 LU 的近似逆经过 (\|I-RA\|<1) 的外包检查，给出严格逆范数。预测误差不超过 \(\|S\|_F\|A^{-1}\|_2\|D\|_F\)，其中 (D) 是原物理状态缺陷。搜索当前文件若在标为真空的像素留有浮点小电流，证明时把这些电流精确设零并重新评价缺陷。材料保持不动。初次旧核对拒绝这些浮点残余的费用也保留。

## 输出


在 18 个物体与上限配对里，最低数值候选均通过了盒及独立正演核验，残差均超过原半径。下表只列该 (u) 的直接候选。全盒上端的嵌套版本可沿用更小盒的更好候选。

| 物体 | (u=1) | (u=1.25) | (u=1.5) | (u=2) | (u=3) |
| --- | ---: | ---: | ---: | ---: | ---: |
| Slotted block | 28.32954 | 27.77924 | 27.48953 | 27.19787 | 27.01649 |
| Four lobes | 33.36212 | — | — | — | 32.47177 |
| Double chevron | 43.23073 | — | — | — | 35.23524 |
| Crescent and plate | 43.15381 | — | — | — | 37.17501 |
| Staggered grooves | 44.09348 | 43.92339 | 41.47159 | 43.68954 | 13.40767 |
| Graded pinwheel | 35.75112 | — | — | — | 26.08989 |

单元格是物理残差严格上端除以原半径。空档使用原任务允许的代表性覆盖。已完成 10 条原生多初值路径、12 条投影起步路径、4 条 Dual-CIE 起步路径。三臂最终阶段合计 426 次物理检查、591,441 个原生浮点 FFT 列。各臂计费分别为 676.44、203.91、99.88 CPU 秒，另有 CIE 种子 8.58 CPU 秒、投影共享建模 0.41 CPU 秒。投影进程结束时正在切换到额外预算，进程总 CPU 约 214.48 秒。18 个严格候选回放合计 145.73 CPU 秒，调用 576 个区间 FFT 列，运行峰值约 1.26 GB。每次共享区间几何初始化、初次真空行拒绝与独立低峰生成对照见总账。

- [全盒上界汇总](upper_endpoints.csv) · [全部审计和费用](upper_endpoints.json)。对未完成或未区间核验的行，表内可采用同对象已核验更小盒的材料或精确真空上界。上端随盒嵌套不会因该层局部优化较弱而升高。`certified_candidate_outside_radius` 仅说明这份材料在球外。
- [独立盒内原生真值对照](strict_endpoints/control_disk075__pw360_32x32__pair2pct.json)。固定 (u=.75)，盘体峰值 .625，复端口严格残差上端为 0.03888818611223809，小于 0.038890131339625975。盘及其独立固定随机扰动与冻结六件目标无关。使用原 Green、B32、S32、噪声规范和原半径护栏。
- [历史生成材料核验](../../../physical_revision/array/forward_audit.json) 提供六个原数据球各自的非空性。六个生成峰值都超过 3，因此不能用作本轮最高 3 的盒内替代材料。

原波形数据与相对阈值衡量拟合，搜索没有使用 GT 评分选材。已验证的高残差候选不支持材料恢复或排除。未使用的照明是独立复测的要求，仅对另一个采集合同才有意义，本轮对冻结 32×32 用的是统一原始观测球。

## 复现

从 Inv_SLAM 根执行。依赖 NumPy、SciPy、python-flint 和 threadpoolctl。独立运行保持单线程。若要保留原始回执，用新的输出目录。逐案均保存实际墙钟、CPU、全波检查、状态及伴随右端、FFT 工作量和迭代信息。

```sh
base=research/z_alg_optical_inertia/ip_unified/numerics/feasibility
export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
python -B "$base/search_native.py" --caps 1 3 2 1.5 1.25 --coarse-steps 4 --middle-steps 5 --full-calls 16 --full-cpu 30
python -B "$base/projected_starts.py" --caps 3 1 --side 16 --coarse-calls 35 --full-calls 3
python -B "$base/cie_starts.py" --caps 1 3 --seed-wall 20 --full-calls 4
python -B "$base/enclose_native.py"
python -B "$base/summarize.py"
```
