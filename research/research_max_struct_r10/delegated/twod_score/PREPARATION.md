# Fresh二维评分脚本已准备，尚未执行

入口：`experiments/score_twod_comparison.py`。当前只完成源码阅读、实现与语法检查；本评分任务没有打开任何fresh GT数组或sealed score metadata，没有生成评分结果或结果图。

运行必须同时具备root明确postseal授权、两队列完全终态、120+48任务完整且各24对象均有七个预定方法。脚本先核验公开任务、queue child结果和每份SEALED预测，再写公开seal；其后才第一次读取`sealed_truth/twod/score_metadata.json`、`generator_metadata.json`和GT数组。不存在自动等待后自行授权或按部分已完成对象先评分的路径。

已实现的冻结口径：

- 主NRMSE沿用sealed amplitude=`max(true chi)−.1`；主IoU为预测/GT均严格`chi>1`。主success只由NRMSE≤.15、fixed IoU≥.7及合法材料组成。
- `legacy_Otsu_IoU`对应历史半幅度GT support；`Otsu_full_support_IoU`对应GT `chi>1`。两者独立命名，不替代主support门。
- low2.5、generator-declared high12/18、background.1各区域RMSE和像素数；空区域null。576维、实数、有限及原box合法性；真实float32预算数组使用单ULP容差。
- CIE-single、Dual-CIE、A+global、A+source、matrixTwofold、matrixGreen19单列。矩阵候选由`PUBLIC_TASK.material_readout.label`辨认，不会并入canonical a_plus。GlobalDual只标记内部matched-current消融。
- 全168行终点、显式失败/缺失、逐对象配对、wide表、low/high/mixed各8例及总体汇总。失败在material success的固定分母中保留；连续指标的完整配对子集与排除数量明确报告。
- 对象paired bootstrap固定seed=20261107、10000次、percentile95% CI，均值与中位差另列；只作描述性证据，没有总体推广或多重比较校正保证。
- 队列GetProcessTimes完整child寿命为可加外部成本；R10/R8/aux/corrector为嵌套视图，不相加。原物理576维LU/RHS与19/28维分母小LU/RHS分开。
- CPU50/100/200/400与physical LU50/100/300/600仅选择最近已付原Phi-incumbent，不看GT选时点。内部真实保存clock、physical LU、外部完整child成本三分；外部预算仅给已完成终点或在clock嵌套条件可证时使用`internal+(external_total−internal_last)`保守上界，未知/未覆盖记null及原因。
- 材料图覆盖全部24对象，按三regime分三页，包含失败标记和共享物理色阶；另有成本/连续误差图。未执行前不生成这些图。

默认结果目录为`delegated/twod_score/scored/`，图为`FIGURES/twod_comparison/`，均write-once。正式执行时应传最终queue receipt（可使用fetch得到的终态receipt覆盖默认参数）；旧RUNNING快照会被拒绝。

root明确授权后才可使用下列调用形式，授权说明必须填写真实收到的指示：

```bash
"/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python" \
  research/research_max_struct_r10/experiments/score_twod_comparison.py \
  --authorize-postseal ROOT_EXPLICIT_POSTSEAL_AUTHORIZATION \
  --authorization-note "填写root实际postseal授权"
```

历史成本说明已轻量更正：六条R9 cie_dual/global/A+与六条新DEV候选都记录同Mac arm64主机及Python3.13.13，可以作历史描述性cost比较；不是同期同负载重跑。两条R6 CIE参考缺完整平台字段，另列unknown。fresh同队列比较仍为主要成本依据。40ac此前也已有另一R9 profile full-current路线恢复，不能将当前DEV称为项目首次解出或已胜所有物理路线。
