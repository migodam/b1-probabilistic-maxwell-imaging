# 盲成像、区域材料界与互耦 MIMO 的复现

[研究网页](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/z_alg_optical_inertia/blind_uq_mimo_theory/) · [原仓库研究目录](https://github.com/migodam/b1-probabilistic-maxwell-imaging/tree/main/research/z_alg_optical_inertia/blind_uq_mimo_theory) · [本轮 Release](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/z-alg-blind-uq-mimo-20261127)

| 附件 | 内容 |
|---|---|
| [研究与复算包](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/z-alg-blind-uq-mimo-20261127/blind-uq-mimo-study.zip) | 公开 32² 和 96² 复数观测、盲图、冻结区域、实际候选、原始统计与全部科学费用和失败轨迹、源码及图 |
| [正演算子包](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/z-alg-blind-uq-mimo-20261127/blind-uq-mimo-operators.zip) | 分别以 n32 和 n96 标识的原网格算子、低秩端口加载与原始物理 E0/S |
| [独立材料评分与生成包](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/z-alg-blind-uq-mimo-20261127/blind-uq-mimo-scoring-truths.zip) | 十件数值材料的参考标签、清洁参考场及原始生成器。只在封存预测后展开 |

三个附件都展开到同一个 blind_uq_mimo_theory/ 根目录，文件名保持原样。研究包内的 archive-manifest-study.json 附逐文件清单，另两个包各有独立清单。研究包包含公开的原发射与接收阵列合同、数值证明和复算成本。Manitoba 第三方原始 ZIP 由采集方独立提供。

## 从冻结结果重算汇总

在可写的复算副本运行以下命令。Python 3.13 的 requirements.txt 记录主研究者原生复算环境。成像配对另用 NumPy 2.5.3、SciPy 1.18.1、python-flint 0.9.0、Matplotlib 3.11.1、threadpoolctl 3.7.0 和 Pillow 12.3.0。严格复算须可导入 python-flint，图形程序可使用 Agg 后端。复算脚本会在各自的 results 目录写入结果。

```sh
unzip blind-uq-mimo-study.zip
unzip blind-uq-mimo-operators.zip
cd blind_uq_mimo_theory
python3.13 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
export PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg
export VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1

python -B imaging/summarize.py
python -B bounds/summarize.py
python -B bounds/roi_mimo_summary.py
python -B diagnostics/aperture_joint.py --aggregate
```

前四条汇总从存档的数组和回执读取。成像比较见 imaging/results/summary.csv，材料界与阵列汇总见 bounds/results/summary.json 和 bounds/roi_mimo_summary.json，交错布局与原布局联合拟合见 results/aperture_joint/summary.csv。results/aperture_control/summary.csv 为单独布局的固定候选和有限预算重拟合结果，cross_audit.csv 收录同一候选对原 A 与留出数据的重新审计。

评分和有真值的图在单独的评分阶段生成。从本目录的上一级展开评分包，回到本目录后运行下列命令。

```sh
cd ..
unzip blind-uq-mimo-scoring-truths.zip
cd blind_uq_mimo_theory
python -B imaging/score.py --sealed-manifest imaging/results/sealed_32.json \
  --private-truth-directory ports/private --output-stem replay_material_scores32
python -B imaging/score.py --sealed-manifest imaging/results/sealed_96.json \
  --private-truth-directory ports/private --output-stem replay_material_scores96
python -B diagnostics/region_summary.py
python -B diagnostics/region_means.py
python -B diagnostics/refresh_report.py
python -B diagnostics/refresh_report.py
python -B imaging/plot_witness.py
python -B diagnostics/report_figures.py
python -B diagnostics/roi_atlas.py
```

两张材料评分表分别按 32² 与 96² 病例封存白名单生成。ROI 的外置信界、已取得的严格材料端点、连续路径及离线覆盖见 results/region_intervals.csv。原始汇总、逐区域观测 ID、核验入口、假阳性定义与成本分别见 [成像结果](imaging/RESULTS_ZH.md)、[区域界结果](bounds/RESULTS_ZH.md)、[阵列结果](ports/PORTS_ZH.md) 和 [数据登记](OBSERVATIONS.md)。

## 审计冻结候选的原全波结果

这些命令选用已有候选和原始离散算子，不重新寻找参数。准确的再计算时间依设备而定。

```sh
python -B imaging/certify.py --stationary n32_c02 n32_c04
python -B imaging/certify_native.py n96_c01 3
python -B imaging/certify_native.py n96_c01 3 --source matched
python -B imaging/certify_native.py n96_c02 1 --source matched --stationary
python -B bounds/verify.py --case n32_c01 --suffix stationary --dark --dark-key roi_01
python -B ports/verify_witness.py --bands 400 --counts 6 8 12 16
python -B diagnostics/aperture_control.py --verify-only
python -B diagnostics/aperture_joint.py --verify-only
```

原生 96² 的分块逆复查需数 GB 内存。随包的 ports/native_verify.py --items 1 2 复查原生真值正演，需先展开评分包。另一路原生 ROI 路径复查用 ports/native_roi_path_verify.py，编译选项按本轮 macOS/Clang 执行合同设置。不同平台可先运行只依赖 Python 数值库的严格复算入口。Flint/Arb 复核针对发布的离散算子与观测，适用的校准和范围见 [理论报告](THEORY_REFINEMENT.md)。

交错 B 的带噪读数可在解封后由 diagnostics/aperture_control.py 重新生成并重拟合，随后的 --cross-audit 从已保存材料复查原 A 和原留出。联合 A+B 的前向、更新和停机可用 diagnostics/aperture_joint.py 重新运行，或用 --cases 分批运行后 --aggregate。比较保留两次采集各自的噪声半径。

较大阵列的区域内端点回放可运行 bounds/roi_mimo.py --q 6 --q 8 --q 12 --q 16 --replay-only。带耦合强度及标定误差的逐单元固定候选复查见 ports/roi_coupling.py，可运行 python -B ports/roi_coupling.py --threads 4。其显式区间乘法也使用 macOS/Clang 执行合同。两类回放对原区域和候选作检验；对应费用和未定单元都保存在结果表内。

## 重新生成合成采集

以下入口会重写复算副本的部分公开观测、算子和代价表，也会在评分目录写出生成用的中间记录。用评分包中的冻结随机生成器依次生成主观测和配对矩阵。阵列生成者、相位与终端来自 [采集设计](ports/DESIGN.md)。

```sh
python -B ports/private/generate.py --side 32 --items 1 2 3 4 5 6 7 8
python -B ports/private/generate.py --side 96 --items 1 2
python -B ports/sweep.py --side 32 --items 1 2 3 4 5 6 7 8
python -B ports/sweep.py --side 96 --items 1 2
```

原生盲成像阶段先使用公开 A 的训练观测。以下示例从头复算一个 32² 对象，再做原系数继续、稳定性控制、物理见证、等更新及等墙钟暗模对照。用于材料评分的预测和选根必须先冻结。从头批跑可传入全部十个病例号。遇到稳定性迭代上限，保留首轮文件并按原协议使用 stationary.py --extend。留出评估为最后一步。

```sh
python -B imaging/run.py baseline n32_c01
python -B imaging/run.py continuation n32_c01
python -B imaging/stationary.py n32_c01
python -B imaging/run.py b1 --stationary n32_c01
python -B imaging/run.py b2 --stationary n32_c01
python -B imaging/run.py held --stationary n32_c01
```

## 原始 Manitoba 数据

本轮发布已派生的空室重复量、对称性残差和复算脚本。原始 ZIP 从 University of Manitoba EIGOR 下载。将路径传给 --archive 或 MANITOBA_TOMOGRAPHIC_ZIP。两条命令读取原始微波记录，分别重算不同阵列下的重复空室波动和原读数互易投影。 [采集方数据页](https://home.cc.umanitoba.ca/~lovetrij/EMILab/mwt_data.html) · [原 Dataset 1 ZIP](https://home.cc.umanitoba.ca/~lovetrij/EMILab/research/EIGOR/Tomographic_Dataset_1.zip) · [采集论文](https://home.cc.umanitoba.ca/~lovetrij/Research/PDFs/Journals/JP50.pdf)

```sh
MANITOBA_TOMOGRAPHIC_ZIP=/path/to/Tomographic_Dataset_1.zip python -B diagnostics/chamber.py
python -B ports/chamber_reciprocity.py --archive /path/to/Tomographic_Dataset_1.zip
```

从完整研究工作区重新打包时，运行 python -B diagnostics/package_research.py。仅刷新报告、图和研究数组时，传入 --part study，算子和评分附件继续沿用。

互耦三档的同一候选池重放使用 ports/roi_coupling.py；对 144 Ω 重新拟合的独立研究用 ports/roi_coupling_refit.py。报告中的两个费用分别计入，详见 [逐阵列与互耦读数](ports/ROI_COUPLING_ZH.md)。
