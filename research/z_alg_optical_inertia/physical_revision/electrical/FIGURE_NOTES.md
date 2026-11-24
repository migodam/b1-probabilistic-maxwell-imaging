# Electrical size and angular sampling figures

## English captions

**Angular interpolation and Manitoba receiver-dark states.** (a) The synthetic full-wave complex responses of six materials at 10, 16 and 24 uniform incident angles are interpolated to the 32-angle responses by a trigonometric interpolant. For an even source count s the fit uses orders −s/2 + 1 through s/2 − 1. The incident amplitude is proportional to 1/√s so all arrays sum to the same incident fluence at equal dwell. The plotted relative Frobenius error compares equal-fluence responses on the same receive ring. Filled markers use 20 receivers and open outlines use 32. The 32-angle references are finite simulations. (b) The curves show the numerical count of negative directions restricted to the receiver null space of one Manitoba antenna color group, read from the `hidden` column. The same real-sector rule is used at every frequency. The normalization is μ = ‖S‖₂²/100 for the corresponding geometry, with zero optical weight. The horizontal line is the 12 incident fields in that group. The upper axis gives frequency and the lower axis gives 2ka for the circumscribed circle of the 0.14 m square. [PDF](electrical_scaling_en.pdf) · [PNG](electrical_scaling_en.png)

**Manitoba large nylon object.** The separate chamber panel shows the relative full-wave distance lower bound at u = 1 and two empirical calibration scales. An ordinate of zero on the lower-bound curve means that this bound supplies no positive distance. The scales come from the two-metal comparison and the empty-scene reference with a separate metal object as a check. On this object, 3 GHz is the only frequency where the lower bound at u = 1 exceeds both empirical scales. [PDF](manitoba_large_nylon_cap_en.pdf) · [PNG](manitoba_large_nylon_cap_en.png)

## 中文图注

**角向插值与 Manitoba 接收盲区的负方向。** （a）对六件仿真材料，取 10、16、24 个等间隔入射角的全波复响应，以三角多项式插值到同一布阵的 32 个入射角，画出相对 Frobenius 误差。偶数源数 s 对应角向阶数从 −s/2 + 1 到 s/2 − 1，单位入射振幅随 1/√s 缩放，等驻留时间下总入射能量通量相同。实心为 20 个接收点，空心轮廓为 32 个。两者分别与本布阵的 32 个入射角响应比较。（b）Manitoba 暗室同一分色组在接收算子零空间中的数值负方向数，读取 `hidden` 字段。各频率使用相同的实扇区权重规则 μ = ‖S‖₂²/100，光学权重为零。虚线标出每组 12 束入射。底轴按 0.14 m 方域的外接圆计算电尺寸 2ka，顶轴是对应频率。 [PDF](electrical_scaling_zh.pdf) · [PNG](electrical_scaling_zh.png)

**Manitoba 大尼龙目标。** 单图展示检验上限 u = 1 时的相对全波距离下界和两个经验标定尺度。一条尺度来自两件金属参考，另一条用空场作参考及独立金属目标作检查。下界为零表示这一次没有得到正距离。对该目标，只有 3 GHz 的下界同时高过两个经验尺度。 [PDF](manitoba_large_nylon_cap_zh.pdf) · [PNG](manitoba_large_nylon_cap_zh.png)

## Generation and inputs

Run the following command from Inv_SLAM with Python, NumPy, and Matplotlib. The plotter reads the existing evidence tables. It does not solve the wave equation or interpolate missing data. A CJK font is needed for the Chinese figures. The PDF retains vector curves and labels. The PNG is rendered at 300 dpi.

```sh
python -B research/z_alg_optical_inertia/physical_revision/electrical/plot_electrical.py
```

- Panel (a) reads all 36 observations in `angular_interpolation.csv` using the `clean_ports_rel_error` column. The table generator cancels the paired 1% and 2% increments by taking 2Y₁ − Y₂ before the angular calculation. [Angular table](angular_interpolation.csv) · [Table generator](angular.py)
- Panel (b) joins the seven Manitoba geometry rows to the receiver-null count at u = 0.50, 0.75, and 1.00. It obtains s = 12 from both tables and checks their agreement. [Geometry](geometry.csv) · [Counts](state_counts.csv) · [Table generator](scale.py)
- The chamber panel uses the seven `nylon_large` rows and checks that the signed margin equals the lower bound minus the larger of the two empirical scales. [Chamber table](chamber_by_frequency.csv) · [Plot source](plot_electrical.py)
